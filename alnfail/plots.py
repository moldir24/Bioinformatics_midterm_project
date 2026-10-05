"""Figures for the report.

One visual language throughout:
  * an aligner always has the same colour, in every figure
  * colour is never the only carrier of identity: every multi-series panel has
    a legend, and the underlying numbers are in results/tables/*.tsv
  * magnitude (heat maps) uses one hue from light to dark, never a rainbow
  * one y-axis per panel; different measures get different panels
  * error rates are drawn on a log axis, because the interesting ones span
    1e-5 to 1e-1 and a linear axis would hide everything below 1 %
"""
from __future__ import annotations

import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# colour follows the tool, not its position in a particular figure
ALIGNER_COLOUR = {
    "minimap2": "#2a78d6",  # blue
    "bwamem2": "#eb6834",  # orange
    "bowtie2": "#1baf7a",  # aqua
    "strobealign": "#eda100",  # yellow
    "winnowmap": "#e87ba4",  # magenta
    "bwa": "#008300",  # green
    "ngmlr": "#4a3aa7",  # violet
}
ALIGNER_MARKER = {"minimap2": "o", "bwamem2": "s", "bowtie2": "^", "strobealign": "D", "winnowmap": "v",
                  "bwa": "P", "ngmlr": "X"}
FALLBACK = "#e34948"
REAL_COLOUR, SIM_COLOUR = "#0b0b0b", "#2a78d6"
SEQUENTIAL = LinearSegmentedColormap.from_list("blue_ramp", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 8.5, "lines.linewidth": 2.0,
    "lines.markersize": 6, "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def colour(aligner: str) -> str:
    return ALIGNER_COLOUR.get(aligner, FALLBACK)


def marker(aligner: str) -> str:
    return ALIGNER_MARKER.get(aligner, "o")


def _legend(fig, aligners, **kwargs):
    handles = [Line2D([0], [0], color=colour(a), marker=marker(a), markeredgecolor=SURFACE, markeredgewidth=1.2,
                      linewidth=2, label=a) for a in aligners]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.02), **kwargs)


def _save(fig, path: str):
    fig.savefig(path)
    plt.close(fig)


def _floor(values, default=1e-5):
    positive = np.asarray(values, float)
    positive = positive[np.isfinite(positive) & (positive > 0)]
    return max(10 ** np.floor(np.log10(positive.min())), 1e-7) if positive.size else default


def _pct(x: float) -> str:
    if x >= 0.1:
        return f"{100 * x:.0f}"
    if x >= 0.01:
        return f"{100 * x:.1f}"
    return f"{100 * x:.2f}" if x >= 0.0001 else ("0" if x == 0 else "<0.01")


# --------------------------------------------------------------------------- #
def accuracy_overview(acc: pd.DataFrame, labels: dict, path: str):
    """Misplaced and unmapped rates per aligner, genome and platform."""
    base = acc[acc["divergence"].astype(float) == 0]
    platforms = [p for p in labels.get("platform_order", sorted(base["platform"].unique())) if p in set(base["platform"])]
    if base.empty or not platforms:
        return
    genomes = [g for g in labels.get("genome_order", sorted(base["genome"].unique())) if g in set(base["genome"])]
    measures = [("rate_misplaced", "n_misplaced", "Misplaced reads (%)"), ("rate_unmapped", "n_unmapped", "Unmapped reads (%)")]
    fig, axes = plt.subplots(2, len(platforms), figsize=(max(1.1 * len(genomes), 3.2) * len(platforms) + 0.8, 6.4),
                             squeeze=False, sharey="row", sharex=True)
    from .stats import wilson_interval

    aligners_all = []
    for col, platform in enumerate(platforms):
        sub = base[base["platform"] == platform]
        aligners = [a for a in ALIGNER_COLOUR if a in set(sub["aligner"])] + sorted(set(sub["aligner"]) - set(ALIGNER_COLOUR))
        aligners_all += [a for a in aligners if a not in aligners_all]
        for row, (rate, count, ylabel) in enumerate(measures):
            ax = axes[row][col]
            for i, aligner in enumerate(aligners):
                s = sub[sub["aligner"] == aligner].set_index("genome").reindex(genomes)
                x = np.arange(len(genomes)) + (i - (len(aligners) - 1) / 2) * 0.16
                lo, hi = wilson_interval(s[count].fillna(0), s["n"].fillna(0))
                y = s[rate].to_numpy(float) * 100
                ax.errorbar(x, np.where(y > 0, y, np.nan), yerr=[np.clip(y - lo * 100, 0, None), np.clip(hi * 100 - y, 0, None)],
                            fmt=marker(aligner), color=colour(aligner), ecolor=colour(aligner), elinewidth=1.0,
                            capsize=0, markeredgecolor=SURFACE, markeredgewidth=1.0, markersize=7)
            ax.set_yscale("log")
            ax.set_xlim(-0.5, len(genomes) - 0.5)
            ax.set_xticks(np.arange(len(genomes)))
            ax.set_xticklabels([textwrap.fill(labels.get("genomes", {}).get(g, g), 16) for g in genomes], fontsize=8)
            ax.grid(axis="x", visible=False)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(labels.get("platforms", {}).get(platform, platform))
    fig.suptitle("Where reads end up: wrong place versus no place", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.925, "Simulated reads with known origin, no added divergence. Bars are 95 % Wilson intervals; "
             "a missing point means no such read was seen.", fontsize=8.5, color=INK2)
    _legend(fig, aligners_all)
    fig.subplots_adjust(top=0.86, bottom=0.2, hspace=0.12, wspace=0.08)
    _save(fig, path)


MAPQ_BANDS = ((0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 60))


def _banded(s: pd.DataFrame, min_reads: int) -> pd.DataFrame:
    """Pool exact MAPQ values into bands; x is the read-weighted mean MAPQ of the band."""
    from .stats import phred, wilson_interval

    rows = []
    for lo, hi in MAPQ_BANDS:
        b = s[(s["mapq"] >= lo) & (s["mapq"] <= hi)]
        n, k = int(b["n_mapped"].sum()), int(b["n_wrong"].sum())
        if n < min_reads:
            continue
        low, high = wilson_interval(k, n)
        rows.append({"x": float((b["mapq"] * b["n_mapped"]).sum() / n), "seen": k > 0,
                     "y": float(phred(k / n)) if k else float(phred(high)),
                     "y_lo": float(phred(high)), "y_hi": float(phred(max(float(low), 1e-7)))})
    return pd.DataFrame(rows)


def mapq_calibration(cal: pd.DataFrame, labels: dict, path: str, min_reads: int = 100):
    """Reported MAPQ against the Phred-scaled error rate actually observed."""
    if cal.empty:
        return
    platforms = [p for p in labels.get("platform_order", sorted(cal["platform"].unique())) if p in set(cal["platform"])]
    fig, axes = plt.subplots(1, len(platforms), figsize=(3.9 * len(platforms) + 0.4, 4.6), squeeze=False, sharey=True)
    aligners_all = []
    for ax, platform in zip(axes[0], platforms):
        sub = cal[cal["platform"] == platform]
        ax.plot([0, 60], [0, 60], color=MUTED, linewidth=1.0)
        ax.text(2, 61, "perfectly calibrated ↗", color=MUTED, fontsize=8, ha="left", va="bottom")
        for aligner in [a for a in ALIGNER_COLOUR if a in set(sub["aligner"])]:
            b = _banded(sub[sub["aligner"] == aligner], min_reads)
            if b.empty:
                continue
            aligners_all += [aligner] if aligner not in aligners_all else []
            seen, unseen = b[b["seen"]], b[~b["seen"]]
            ax.errorbar(seen["x"], seen["y"], yerr=[seen["y"] - seen["y_lo"], seen["y_hi"] - seen["y"]],
                        color=colour(aligner), marker=marker(aligner), markeredgecolor=SURFACE, markeredgewidth=1.0,
                        linewidth=1.4, elinewidth=1.0, capsize=0)
            # no error observed: the true value is at least this high (a lower bound)
            ax.scatter(unseen["x"], unseen["y"], marker=marker(aligner), facecolors=SURFACE,
                       edgecolors=colour(aligner), linewidths=1.4, s=40, zorder=3)
        ax.set_title(labels.get("platforms", {}).get(platform, platform))
        ax.set_xlabel("MAPQ reported by the aligner")
        ax.set_xlim(-2, 63)
        ax.set_ylim(-2, 66)
    axes[0][0].set_ylabel("MAPQ earned: −10·log10(observed error rate)")
    fig.suptitle("Is the reported mapping quality honest?", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.9, f"Simulated reads pooled into MAPQ bands (0, 1–9, 10–19, … 60; at least {min_reads} reads per band), with 95 % "
             "intervals.\nBelow the line = over-confident. Hollow points: no error seen, so the true value lies at or above the point.",
             fontsize=8.5, color=INK2, va="top")
    _legend(fig, aligners_all)
    fig.subplots_adjust(top=0.78, bottom=0.24, wspace=0.08)
    _save(fig, path)


def mapq_threshold(cal: pd.DataFrame, labels: dict, path: str):
    """Trade-off when filtering by MAPQ: reads kept versus errors kept."""
    if cal.empty:
        return
    totals = labels.get("reads_per_platform_aligner", {})
    platforms = [p for p in labels.get("platform_order", sorted(cal["platform"].unique())) if p in set(cal["platform"])]
    fig, axes = plt.subplots(1, len(platforms), figsize=(3.9 * len(platforms) + 0.4, 4.4), squeeze=False, sharey=True)
    aligners_all = []
    for ax, platform in zip(axes[0], platforms):
        sub = cal[cal["platform"] == platform]
        for aligner in [a for a in ALIGNER_COLOUR if a in set(sub["aligner"])]:
            s = sub[sub["aligner"] == aligner].sort_values("mapq", ascending=False)
            aligners_all += [aligner] if aligner not in aligners_all else []
            total = totals.get(f"{platform}|{aligner}", s["n_mapped"].sum())
            y = s["cum_error"].replace(0, np.nan)
            ax.plot(s["cum_mapped"] / total * 100, y, color=colour(aligner), linewidth=1.8)
            at30 = s[s["mapq"] >= 30].tail(1)
            if len(at30) and at30["cum_error"].iloc[0] > 0:
                ax.plot(at30["cum_mapped"] / total * 100, at30["cum_error"], marker=marker(aligner), color=colour(aligner),
                        markeredgecolor=SURFACE, markeredgewidth=1.2, markersize=8)
        ax.axhline(1e-3, color=MUTED, linewidth=1.0)
        ax.text(0.98, 0.02, "grey line: 1 in 1,000,\nthe MAPQ 30 promise", transform=ax.transAxes, color=MUTED,
                fontsize=8, va="bottom", ha="right")
        ax.set_yscale("log")
        ax.set_title(labels.get("platforms", {}).get(platform, platform))
        ax.set_xlabel("Reads kept (% of all reads)")
    axes[0][0].set_ylabel("Share of kept reads that are misplaced")
    fig.suptitle("What a MAPQ filter buys", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.9, "Each curve lowers the MAPQ threshold from 60 (left) to 0 (right). Markers: threshold MAPQ ≥ 30. "
             "A curve starts where its first error appears\n(zero cannot be drawn on a log axis), so an aligner without a "
             "curve made no error at any threshold.", fontsize=8.5, color=INK2, va="top")
    _legend(fig, aligners_all)
    fig.subplots_adjust(top=0.78, bottom=0.25, wspace=0.08)
    _save(fig, path)


ORDERED_SETS = {
    "gc_content": ["<30%", "30-40%", "40-50%", "50-60%", "60-70%", ">=70%"],
    "variant_density": ["0", "0-2/kb", "2-5/kb", "5-10/kb", ">=10/kb"],
    "read_length": ["<500", "0.5-1k", "1-2k", "2-5k", "5-10k", "10-20k", "20-50k", ">=50k"],
    "repeat_age": ["young_<5%", "middle_5-15%", "old_>15%"],
    "segdup": ["SD_<95%", "SD_95-99%", "SD_>=99%", "SD"],
}
SET_TITLES = {"gc_content": "GC content", "variant_density": "variants", "read_length": "read length",
              "repeat_class": "repeat class", "repeat_age": "repeat age (divergence from consensus)",
              "segdup": "segmental duplication", "giab_confident": "benchmark regions", "feature": "annotated feature"}


def _strata_rows(sub: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Choose and order the heat-map rows.

    Kept: every label of every stratum set, the 'none' baseline of the repeat
    and duplication sets (reads in unique sequence), and for those sets the
    reads lying entirely inside one element. Numeric sets keep their natural
    order; the others are sorted by name.
    """
    cover = sub["stratum_set"].str.endswith("_cover")
    base_set = sub["stratum_set"].str.replace("_cover", "", regex=False)
    keep = (~cover & (sub["stratum"] != "none")) | (~cover & base_set.isin(["repeat_class", "segdup"])) \
        | (cover & (sub["stratum"] == "entire_>=95%") & base_set.isin(["repeat_class", "segdup", "feature"]))
    sub = sub[keep].copy()
    if sub.empty:
        return sub, []
    base_set = sub["stratum_set"].str.replace("_cover", "", regex=False)
    title = base_set.map(lambda s: SET_TITLES.get(s, s.replace("_", " ")))
    label = np.where(sub["stratum_set"].str.endswith("_cover"), "read entirely inside one element",
                     np.where(sub["stratum"] == "none", "none (unique sequence)", sub["stratum"].str.replace("_", " ")))
    sub["row"] = title + ": " + label

    def rank(r):
        base = r["stratum_set"].replace("_cover", "")
        if r["stratum_set"].endswith("_cover"):
            position = 999
        elif r["stratum"] == "none":
            position = -1
        elif base in ORDERED_SETS and r["stratum"] in ORDERED_SETS[base]:
            position = ORDERED_SETS[base].index(r["stratum"])
        else:
            position = 500
        return (base, position, r["stratum"])

    firsts = sub.drop_duplicates("row")
    order = [row for _, row in sorted(((rank(r), r["row"]) for _, r in firsts.iterrows()))]
    return sub, order


def strata_heatmap(strata: pd.DataFrame, genome: str, platform: str, labels: dict, path: str, min_reads: int = 100):
    """The failure map: misplacement rate per genomic context and aligner."""
    sub = strata[(strata["genome"] == genome) & (strata["platform"] == platform) & (strata["divergence"].astype(float) == 0)]
    sub = sub[sub["n"] >= min_reads]
    sub, order = _strata_rows(sub)
    if sub.empty:
        return
    aligners = [a for a in ALIGNER_COLOUR if a in set(sub["aligner"])]
    panels = [("rate_misplaced", "Misplaced (% of reads)"), ("rate_confident_wrong", "Misplaced among MAPQ ≥ 30 (%)"),
              ("rate_unmapped", "Unmapped (% of reads)")]
    fig, axes = plt.subplots(1, 3, figsize=(3.1 * 3 + 3.2, 0.27 * len(order) + 2.2), squeeze=False, sharey=True)
    for ax, (value, title) in zip(axes[0], panels):
        grid = sub.pivot_table(index="row", columns="aligner", values=value, aggfunc="first").reindex(index=order, columns=aligners)
        data = grid.to_numpy(float)
        shown = np.where(data > 0, data, np.nan)
        norm = LogNorm(vmin=1e-4, vmax=1.0)
        ax.imshow(np.where(np.isnan(shown), 1e-4, np.clip(shown, 1e-4, 1)), cmap=SEQUENTIAL, norm=norm, aspect="auto")
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                v = data[i, j]
                if np.isnan(v):
                    ax.text(j, i, "–", ha="center", va="center", color=MUTED, fontsize=8)
                    continue
                ax.text(j, i, _pct(v), ha="center", va="center", fontsize=8,
                        color="#ffffff" if v > 0.02 else INK)
        ax.set_xticks(range(len(aligners)))
        ax.set_xticklabels(aligners, rotation=30, ha="right")
        ax.set_title(title, fontsize=9.5)
        ax.grid(False)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        # white gaps between cells instead of borders
        ax.set_xticks(np.arange(-0.5, len(aligners)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(order)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        ax.tick_params(which="minor", length=0)
    axes[0][0].set_yticks(range(len(order)))
    axes[0][0].set_yticklabels(order)
    name = f"{labels.get('genomes', {}).get(genome, genome)}, {labels.get('platforms', {}).get(platform, platform)}"
    fig.suptitle(f"Where each aligner breaks: {name}", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 1 - 0.75 / fig.get_figheight(), f"Percent of simulated reads per context (contexts with at least {min_reads} reads). "
             "Darker = worse. – = not enough reads.", fontsize=8.5, color=INK2)
    fig.subplots_adjust(top=1 - 1.25 / fig.get_figheight(), wspace=0.05)
    _save(fig, path)


def divergence_curves(acc: pd.DataFrame, labels: dict, path: str):
    """Accuracy as the sample diverges from the reference."""
    acc = acc.assign(d=acc["divergence"].astype(float))
    series = acc.groupby(["genome", "platform"])["d"].nunique()
    keep = series[series > 1].index.tolist()
    if not keep:
        return
    order_g, order_p = labels.get("genome_order", []), labels.get("platform_order", [])
    keep.sort(key=lambda gp: (order_g.index(gp[0]) if gp[0] in order_g else 99,
                              order_p.index(gp[1]) if gp[1] in order_p else 99))
    measures = [("rate_correct", "Correctly placed (%)"), ("rate_misplaced", "Misplaced (%)"), ("rate_unmapped", "Unmapped (%)")]
    fig, axes = plt.subplots(len(measures), len(keep), figsize=(3.7 * len(keep) + 0.5, 2.5 * len(measures) + 1.2),
                             squeeze=False, sharex=True, sharey="row")
    aligners_all = []
    for col, (genome, platform) in enumerate(keep):
        sub = acc[(acc["genome"] == genome) & (acc["platform"] == platform)]
        for row, (value, ylabel) in enumerate(measures):
            ax = axes[row][col]
            for aligner in [a for a in ALIGNER_COLOUR if a in set(sub["aligner"])]:
                s = sub[sub["aligner"] == aligner].sort_values("d")
                aligners_all += [aligner] if aligner not in aligners_all else []
                ax.plot(s["d"] * 100, s[value] * 100, color=colour(aligner), marker=marker(aligner),
                        markeredgecolor=SURFACE, markeredgewidth=1.0)
            if row == 0:
                ax.set_title(f"{labels.get('genomes', {}).get(genome, genome)}\n{labels.get('platforms', {}).get(platform, platform)}")
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == len(measures) - 1:
                ax.set_xlabel("Divergence from the reference (%)")
    fig.suptitle("Mapping reads from a genome that is not the reference", x=0.02, ha="left", fontsize=12, fontweight="bold")
    _legend(fig, aligners_all)
    fig.subplots_adjust(top=0.88, bottom=0.14, hspace=0.18, wspace=0.08)
    _save(fig, path)


def calibration_short(real: dict, sim: dict, title: str, path: str):
    """Per-cycle error and quality: real reads against simulated reads."""
    mates = min(real["cycle_err"].shape[0], sim["cycle_err"].shape[0])
    fig, axes = plt.subplots(2, mates, figsize=(4.6 * mates + 0.4, 5.6), squeeze=False, sharex=True, sharey="row")
    for m in range(mates):
        n = int(min(np.flatnonzero(real["length_probs"][m]).max(), np.flatnonzero(sim["length_probs"][m]).max()))
        x = np.arange(1, n + 1)
        for row, (key, ylabel, scale) in enumerate((("cycle_err", "Mismatch rate (%)", 100), ("cycle_meanq", "Mean base quality (Phred)", 1))):
            ax = axes[row][m]
            ax.plot(x, real[key][m][:n] * scale, color=REAL_COLOUR, linewidth=1.6, label="real reads")
            ax.plot(x, sim[key][m][:n] * scale, color=SIM_COLOUR, linewidth=1.6, label="simulated reads")
            if m == 0:
                ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(f"Read {m + 1}")
            else:
                ax.set_xlabel("Sequencing cycle (position in read)")
    axes[0][0].legend(loc="upper left")
    fig.suptitle(f"Does the simulator reproduce the instrument? {title}", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.subplots_adjust(top=0.88, hspace=0.12, wspace=0.06)
    _save(fig, path)


def calibration_long(real: dict, sim: dict, title: str, path: str):
    """Read length, per-read error and error composition: real against simulated."""
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    for ax, key, xlabel, logx in ((axes[0], "lengths", "Read length (bp)", True), (axes[1], "errors", "Per-read error rate (%)", True)):
        scale = 100 if key == "errors" else 1
        r, s = real[key] * scale, sim[key] * scale
        lo, hi = max(min(r.min(), s.min()), 1e-3), max(r.max(), s.max())
        bins = np.logspace(np.log10(lo), np.log10(hi), 45)
        ax.hist(r, bins=bins, histtype="step", color=REAL_COLOUR, linewidth=1.6, density=True, label="real reads")
        ax.hist(s, bins=bins, histtype="step", color=SIM_COLOUR, linewidth=1.6, density=True, label="simulated reads")
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_yticklabels([])
    fig.legend(handles=[Line2D([0], [0], color=REAL_COLOUR, linewidth=1.6, label="real reads"),
                        Line2D([0], [0], color=SIM_COLOUR, linewidth=1.6, label="simulated reads")],
               loc="upper right", ncol=2, bbox_to_anchor=(0.9, 0.99))
    kinds = ["mismatch_rate", "ins_rate", "del_rate"]
    x = np.arange(3)
    axes[2].bar(x - 0.17, [real["meta"][k] * 100 for k in kinds], width=0.3, color=REAL_COLOUR)
    axes[2].bar(x + 0.17, [sim["meta"][k] * 100 for k in kinds], width=0.3, color=SIM_COLOUR)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(["substitutions", "insertions", "deletions"])
    axes[2].set_ylabel("Errors per 100 aligned bases")
    axes[2].grid(axis="x", visible=False)
    hp = np.arange(1, real["hp_del_factor"].size)
    for key, style, name in (("hp_del_factor", "-", "deletions"), ("hp_ins_factor", ":", "insertions")):
        axes[3].plot(hp, real[key][1:], color=REAL_COLOUR, linestyle=style, linewidth=1.6)
        axes[3].plot(hp, sim[key][1:], color=SIM_COLOUR, linestyle=style, linewidth=1.6)
    axes[3].set_xlabel("Homopolymer length (bp, last bin pooled)")
    axes[3].set_ylabel("Indel rate relative to genome average")
    axes[3].legend(handles=[Line2D([0], [0], color=INK2, linestyle="-", label="deletions"),
                            Line2D([0], [0], color=INK2, linestyle=":", label="insertions")], loc="upper left")
    fig.suptitle(f"Does the simulator reproduce the instrument? {title}", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.subplots_adjust(top=0.84, wspace=0.28)
    _save(fig, path)


def transfer_scatter(points: pd.DataFrame, path: str):
    """Per-stratum disagreement: simulation (x) against real reads (y)."""
    if points.empty:
        return
    pairs = points[["simulated", "real"]].drop_duplicates().values.tolist()
    per_pair = [[a for a in ALIGNER_COLOUR if a in set(points.loc[points["real"] == real_id, "aligner"])]
                for _, real_id in pairs]
    n_cols = max(len(a) for a in per_pair)
    fig, axes = plt.subplots(len(pairs), n_cols, figsize=(2.9 * n_cols + 1.0, 2.9 * len(pairs) + 1.1),
                             squeeze=False, sharex=True, sharey=True)
    floor = _floor(np.concatenate([points["sim_disagreement"], points["real_disagreement"]]), 1e-4) / 3
    for r, ((sim_id, real_id), aligners) in enumerate(zip(pairs, per_pair)):
        for c in range(n_cols):
            ax = axes[r][c]
            if c >= len(aligners):
                ax.set_visible(False)
                continue
            aligner = aligners[c]
            s = points[(points["simulated"] == sim_id) & (points["real"] == real_id) & (points["aligner"] == aligner)]
            ax.plot([floor, 1], [floor, 1], color=MUTED, linewidth=1.0)
            ax.scatter(np.maximum(s["sim_disagreement"], floor), np.maximum(s["real_disagreement"], floor), s=26,
                       color=colour(aligner), edgecolors=SURFACE, linewidths=0.8, zorder=3)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(aligner, fontsize=9.5)
            if c == 0:
                ax.set_ylabel(f"{real_id}\nreal reads", fontsize=8)
            ax.set_xlabel("simulated reads")
            ax.tick_params(labelbottom=True)
    fig.suptitle("Do the same genomic contexts fail on real reads?", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 1 - 0.62 / fig.get_figheight(), "Share of reads placed away from the aligner consensus. One dot per genomic context; "
             "on the line = same rate in simulation and reality.\nDots on the lower or left edge are contexts with no disagreement at all.",
             fontsize=8.5, color=INK2, va="top")
    fig.subplots_adjust(top=1 - 1.35 / fig.get_figheight(), wspace=0.08, hspace=0.45)
    _save(fig, path)


def resources(bench: pd.DataFrame, labels: dict, path: str):
    """Wall-clock time and peak memory per aligner (separate panels, one axis each)."""
    maps = bench[bench["step"] == "map"].copy()
    if maps.empty:
        return
    maps["platform"] = maps["target"].map(labels.get("dataset_platform", {}))
    maps = maps.dropna(subset=["platform"])
    agg = maps.groupby(["platform", "group"]).agg(seconds=("s", "median"), rss=("max_rss", "max")).reset_index()
    platforms = [p for p in labels.get("platform_order", sorted(agg["platform"].unique())) if p in set(agg["platform"])]
    fig, axes = plt.subplots(2, len(platforms), figsize=(3.4 * len(platforms) + 0.6, 5.4), squeeze=False, sharey="row")
    for col, platform in enumerate(platforms):
        s = agg[agg["platform"] == platform]
        s = s.set_index("group").reindex([a for a in ALIGNER_COLOUR if a in set(s["group"])]).reset_index()
        for row, (value, ylabel) in enumerate((("seconds", "Median wall time per dataset (s)"), ("rss", "Peak memory (MB)"))):
            ax = axes[row][col]
            ax.bar(s["group"], s[value], width=0.5, color=[colour(a) for a in s["group"]])
            for x, v in enumerate(s[value]):
                ax.text(x, v, f"{v:,.0f}", ha="center", va="bottom", fontsize=8, color=INK2)
            ax.grid(axis="x", visible=False)
            ax.tick_params(axis="x", rotation=30)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(labels.get("platforms", {}).get(platform, platform))
    fig.suptitle("What each aligner costs", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.subplots_adjust(top=0.88, hspace=0.3, wspace=0.08)
    _save(fig, path)


def mapq_theory(ideal: pd.DataFrame, path: str):
    """The exact posterior is calibrated: claimed against observed error."""
    s = ideal.groupby("mapq")[["n_reads", "n_wrong"]].sum().reset_index()
    s = s[(s["n_wrong"] > 0) & (s["n_reads"] >= 100)]
    fig, ax = plt.subplots(figsize=(4.6, 4.3))
    ax.plot([0, 60], [0, 60], color=MUTED, linewidth=1.0)
    ax.plot(s["mapq"], -10 * np.log10(s["n_wrong"] / s["n_reads"]), marker="o", color="#2a78d6", markeredgecolor=SURFACE,
            markeredgewidth=1.2, linewidth=0)
    ax.set_xlabel("MAPQ from the exact posterior")
    ax.set_ylabel("MAPQ earned: −10·log10(observed error rate)")
    ax.set_title("An ideal aligner sits on the diagonal")
    ax.set_xlim(-2, 62)
    ax.set_ylim(-2, 62)
    _save(fig, path)
