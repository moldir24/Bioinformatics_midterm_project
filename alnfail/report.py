"""Assemble figures and a results digest from the collected tables.

`summary.md` restates the numbers; it does not interpret them. Interpretation
(why a tool fails where it fails, and what that means for a user) belongs in
the written report and is the authors' work.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from . import plots
from .profile import load_profile
from .stats import wilson_interval


def _read(tables: str, name: str) -> pd.DataFrame:
    path = os.path.join(tables, f"{name}.tsv")
    if not os.path.exists(path) or os.path.getsize(path) < 2:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype={"stratum": str, "key": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _md_table(frame: pd.DataFrame, formats: dict | None = None) -> str:
    """A DataFrame as a GitHub-flavoured Markdown table."""
    if frame.empty:
        return "_no rows_\n"
    formats = formats or {}
    lines = ["| " + " | ".join(frame.columns) + " |", "|" + "|".join("---" for _ in frame.columns) + "|"]
    for _, row in frame.iterrows():
        cells = []
        for col in frame.columns:
            value = row[col]
            if col in formats and pd.notna(value):
                cells.append(formats[col].format(value))
            elif isinstance(value, float):
                cells.append("" if np.isnan(value) else f"{value:.4g}")
            else:
                cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _per_thousand(k, n) -> str:
    """'1.8 (1.2-2.6)' errors per 1,000 with a 95 % Wilson interval."""
    if n == 0:
        return "no reads"
    lo, hi = wilson_interval(k, n)
    return f"{1000 * k / n:.2f} ({1000 * float(lo):.2f}–{1000 * float(hi):.2f})"


def build_report(tables: str, outdir: str, meta_path: str | None = None, benchmarks: str | None = None) -> None:
    os.makedirs(outdir, exist_ok=True)
    meta = json.load(open(meta_path)) if meta_path else {}
    labels = meta.get("labels", {})
    acc, strata = _read(tables, "accuracy"), _read(tables, "strata")
    cal = _read(tables, "mapq_calibration_pooled")
    head = _read(tables, "headline")
    figures: list[tuple[str, str]] = []

    def figure(name: str, caption: str, func, *args):
        path = os.path.join(outdir, name)
        func(*args, path)
        if os.path.exists(path):
            figures.append((name, caption))

    if not head.empty:  # denominators for "reads kept": all simulated reads, mapped or not
        labels["reads_per_platform_aligner"] = {f"{r.platform}|{r.aligner}": int(r.n) for r in head.itertuples()}
    if not acc.empty:
        figure("fig_accuracy.png", "Misplaced and unmapped reads per aligner, genome and platform (task 4-6).",
               plots.accuracy_overview, acc, labels)
        figure("fig_divergence.png", "Accuracy as the sequenced genome diverges from the reference (task 7).",
               plots.divergence_curves, acc, labels)
    if not cal.empty:
        figure("fig_mapq_calibration.png", "Reported MAPQ against the error rate actually observed (task 6).",
               plots.mapq_calibration, cal, labels)
        figure("fig_mapq_threshold.png", "Reads kept and errors kept as the MAPQ threshold is lowered (task 6).",
               plots.mapq_threshold, cal, labels)
    if not strata.empty:
        for (genome, platform), _ in strata.groupby(["genome", "platform"]):
            figure(f"fig_strata_{genome}_{platform}.png",
                   f"Failure map for {labels.get('genomes', {}).get(genome, genome)}, {platform} (task 7).",
                   lambda path, g=genome, p=platform: plots.strata_heatmap(strata, g, p, labels, path))

    for item in meta.get("profiles", []):  # real against simulated error profile (task 3)
        if not (os.path.exists(item["real"]) and os.path.exists(item["sim"])):
            continue
        real, sim = load_profile(item["real"]), load_profile(item["sim"])
        func = plots.calibration_short if real["meta"]["kind"] == "short" else plots.calibration_long
        title = labels.get("platforms", {}).get(item["platform"], item["platform"])
        figure(f"fig_calibration_{item['dataset']}.png", f"Simulator against real reads: {title} (task 3).",
               lambda path, r=real, s=sim, t=title, f=func: f(r, s, t, path))

    points = _read(tables, "transfer_strata_points")
    if not points.empty:
        figure("fig_transfer.png", "Per-context disagreement with the aligner consensus, simulation against real reads (task 8).",
               plots.transfer_scatter, points)
    bench = pd.read_csv(benchmarks, sep="\t") if benchmarks and os.path.exists(benchmarks) and os.path.getsize(benchmarks) > 2 else pd.DataFrame()
    if not bench.empty:
        figure("fig_resources.png", "Wall time and peak memory per aligner.", plots.resources, bench, labels)
    ideal_path = os.path.join(tables, "mapq_ideal_calibration.tsv")
    if os.path.exists(ideal_path):
        figure("fig_mapq_theory.png", "An aligner that reports the exact posterior is perfectly calibrated (task 1).",
               plots.mapq_theory, pd.read_csv(ideal_path, sep="\t"))

    # ---------------------------------------------------------------- digest
    out = ["# Results digest", "",
           "Generated by `python -m alnfail report`. Every number comes from `results/tables/*.tsv`; "
           "rates in brackets are 95 % Wilson intervals. This file states what was measured. "
           "Explaining it is the job of the written report.", ""]

    if not head.empty:
        out += ["## 1. Headline: placement accuracy and the MAPQ 30 promise", "",
                "Simulated reads with known origin, all genomes pooled, no added divergence. "
                "A MAPQ of 30 claims at most 1 wrong alignment per 1,000.", ""]
        rows = []
        for r in head.sort_values(["platform", "aligner"]).itertuples():
            rows.append({
                "platform": r.platform, "aligner": r.aligner, "reads": f"{r.n:,}",
                "correct %": f"{100 * r.rate_correct:.2f}", "misplaced %": f"{100 * r.rate_misplaced:.3f}",
                "unmapped %": f"{100 * r.rate_unmapped:.3f}",
                "reads at MAPQ≥30 %": f"{100 * r.n_confident / max(r.n, 1):.1f}",
                "wrong per 1,000 at MAPQ≥30": _per_thousand(r.n_confident_wrong, r.n_confident),
                "× the promised 1 per 1,000": "" if pd.isna(r.excess_over_claim) else f"{r.excess_over_claim:.2f}",
            })
        out += [_md_table(pd.DataFrame(rows)), ""]

    if not cal.empty:
        out += ["## 2. Calibration by reported MAPQ", "",
                "For each band of reported MAPQ: how many alignments, how many were wrong, and the error rate "
                "the band is supposed to have at most (its lower edge).", ""]
        bands = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 59), (60, 60)]
        rows = []
        for (platform, aligner), g in cal.groupby(["platform", "aligner"]):
            for lo, hi in bands:
                b = g[(g["mapq"] >= lo) & (g["mapq"] <= hi)]
                n, k = int(b["n_mapped"].sum()), int(b["n_wrong"].sum())
                if n == 0:
                    continue
                rows.append({"platform": platform, "aligner": aligner, "MAPQ band": f"{lo}" if lo == hi else f"{lo}–{hi}",
                             "alignments": f"{n:,}", "wrong": f"{k:,}",
                             "wrong per 1,000": _per_thousand(k, n),
                             "claimed at most per 1,000": f"{1000 * 10 ** (-lo / 10):.3g}"})
        out += [_md_table(pd.DataFrame(rows)), ""]

    if not strata.empty:
        out += ["## 3. Hardest genomic contexts", "",
                "The ten contexts with the highest misplacement rate per platform and genome "
                "(at least 200 simulated reads; no added divergence), with the rate of every aligner.", ""]
        base = strata[(strata["divergence"].astype(float) == 0) & (strata["n"] >= 200)]
        for (genome, platform), g in base.groupby(["genome", "platform"]):
            wide = g.pivot_table(index=["stratum_set", "stratum"], columns="aligner", values="rate_misplaced", aggfunc="first")
            wide = (wide * 100).round(3)
            wide["worst"] = wide.max(axis=1)
            top = wide.sort_values("worst", ascending=False).head(10).drop(columns="worst").reset_index()
            out += [f"**{labels.get('genomes', {}).get(genome, genome)}, {platform}** (misplaced, % of reads)", "",
                    _md_table(top), ""]

    check = _read(tables, "calibration_check")
    if not check.empty:
        out += ["## 4. Is the simulator realistic? (task 3)", "",
                "Each row compares a quantity measured on real reads with the same quantity measured, by the same "
                "procedure, on simulated reads. `ok = False` rows are limitations of the simulation.", ""]
        out += [_md_table(check[["dataset", "metric", "unit", "real", "simulated", "difference", "kind", "ok"]]), ""]
        bad = check[~check["ok"].astype(bool)]
        out += [f"{len(check) - len(bad)} of {len(check)} checks inside tolerance."
                + ("" if bad.empty else " Outside tolerance: " + "; ".join(f"{r.dataset}: {r.metric}" for r in bad.itertuples()) + "."), ""]

    t_obs, t_str, prov = _read(tables, "transfer_observables"), _read(tables, "transfer_strata"), _read(tables, "provider_agreement")
    if not t_obs.empty or not t_str.empty:
        out += ["## 5. Does the simulation transfer to real reads? (task 8)", ""]
        if not t_obs.empty:
            bad = t_obs[~t_obs["transfers"].astype(bool)]
            out += [f"Truth-free quantities (mapped fraction, MAPQ distribution, edit distance, clipping, split reads), "
                    f"simulated against real: {len(t_obs) - len(bad)} of {len(t_obs)} comparisons agree within tolerance. "
                    "All rows are in `transfer_observables.tsv`; the ones the simulation does **not** reproduce are:", "",
                    _md_table(bad) if len(bad) else "_none_\n", ""]
        if not t_str.empty:
            out += ["Rank agreement across genomic contexts. `proxy_validity` asks whether disagreement with the "
                    "consensus tracks true error in simulation; `transfer` asks whether the contexts that disagree "
                    "in simulation also disagree on real reads.", "", _md_table(t_str), ""]
        if not prov.empty:
            prov = prov.assign(agree_pct=100 * prov["n_agree_provider"] / prov["n_compared"].clip(lower=1),
                               confident_disagree_per_1000=1000 * prov["n_confident_disagree_provider"] / prov["n_confident"].clip(lower=1))
            out += ["Agreement with the data provider's own whole-genome alignment of the same reads.", "",
                    _md_table(prov[["dataset", "aligner", "n_compared", "agree_pct", "confident_disagree_per_1000"]]), ""]

    qc = _read(tables, "qc_decisions")
    if not qc.empty:
        out += ["## 6. Read QC decisions", "", _md_table(qc[qc["step"].str.startswith("filter")][
            ["dataset", "step", "metric", "value", "threshold", "why"]]), ""]
    if not bench.empty:
        b = bench[bench["step"].isin(["map", "index"])]
        if not b.empty:
            agg = b.groupby(["step", "group"]).agg(runs=("s", "size"), median_seconds=("s", "median"),
                                                   total_cpu_seconds=("cpu_time", "sum"), peak_rss_mb=("max_rss", "max")).reset_index()
            out += ["## 7. Runtime and memory", "", _md_table(agg.rename(columns={"group": "aligner"})), ""]
    if figures:
        out += ["## Figures", ""] + [f"- `{name}`: {caption}" for name, caption in figures] + [""]
    with open(os.path.join(outdir, "summary.md"), "w") as fh:
        fh.write("\n".join(out))
    print(f"{len(figures)} figures and summary.md written to {outdir}")
