"""From per-read outcomes to the three result tables (tasks 6 and 7).

    accuracy          how often each aligner is correct / misplaced / unmapped
    mapq_calibration  for every reported MAPQ value, how often the alignment
                      was actually wrong
    strata            the same outcomes broken down by genomic context

All tables are tidy (one observation per row) and carry raw counts, so that
any rate can be recomputed and given a confidence interval later.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pysam

from .annot import StratumSet, count_in_intervals, gc_bin, length_bin, variant_density_bin
from .seqio import gc_fraction, to_array
from .stats import phred, wilson_interval

MAPQ_CONFIDENT = 30  # the threshold most variant-calling pipelines use to trust an alignment


def cover_bin(frac: np.ndarray) -> np.ndarray:
    """How much of the read lies inside the annotated element."""
    labels = np.array(["none", "partial_<50%", "mostly_50-95%", "entire_>=95%"], dtype=object)
    idx = np.where(frac <= 0, 0, 1 + np.searchsorted(np.array([0.5, 0.95]), frac, side="right"))
    return labels[idx]


def annotate_intervals(df: pd.DataFrame, strata: dict[str, StratumSet], fasta_path: str | None = None,
                       variant_positions: dict[str, np.ndarray] | None = None, long_reads: bool = False,
                       divergence: float | None = None) -> pd.DataFrame:
    """Attach stratum labels to a frame of intervals (contig, start, end).

    Used for the true intervals of simulated reads and, in task 8, for the
    consensus placement of real reads, so both are stratified by one piece of
    code.
    """
    out = df.copy()
    contigs = out["contig"].astype(str).to_numpy(dtype=object)
    starts = out["start"].to_numpy(np.int64)
    ends = out["end"].to_numpy(np.int64)
    span = np.maximum(ends - starts, 1)
    for name, sset in strata.items():
        label, frac = sset.assign(contigs, starts, ends)
        if len(sset.labels) == 1:  # a yes/no set: say what "no" means instead of 'none'
            label = np.where(label == "none", f"outside_{sset.labels[0]}", label)
        out[name] = label
        out[f"{name}_cover"] = cover_bin(frac)
    if "gc" not in out.columns:
        if fasta_path is None:
            raise ValueError("need the reference FASTA to compute GC content")
        fasta = pysam.FastaFile(fasta_path)
        out["gc"] = [gc_fraction(to_array(fasta.fetch(c, int(s), int(e)))) if e > s else np.nan
                     for c, s, e in zip(contigs, starts, ends)]
    out["gc_content"] = gc_bin(out["gc"].to_numpy(float))
    if variant_positions is not None:
        out["var_count"] = count_in_intervals(variant_positions, contigs, starts, ends)
        out["variant_density"] = variant_density_bin(out["var_count"].to_numpy(), span)
    if long_reads:
        out["read_length"] = length_bin(span)
    if divergence is not None:
        out["divergence"] = f"{divergence:g}"
    return out


def stratum_columns(df: pd.DataFrame, strata_names: list[str]) -> list[str]:
    """The columns of an annotated frame that define strata."""
    cols = []
    for name in strata_names:
        cols += [name, f"{name}_cover"]
    cols += [c for c in ("gc_content", "variant_density", "read_length") if c in df.columns]
    return [c for c in cols if c in df.columns]


def _counts(frame: pd.DataFrame) -> dict:
    outcome = frame["outcome"].to_numpy()
    mapq = frame["mapq"].to_numpy()
    mapped = outcome != "unmapped"
    wrong = outcome == "misplaced"
    confident = mapped & (mapq >= MAPQ_CONFIDENT)
    return {
        "n": int(outcome.size),
        "n_correct": int((outcome == "correct").sum()),
        "n_misplaced": int(wrong.sum()),
        "n_unmapped": int((~mapped).sum()),
        "n_confident": int(confident.sum()),
        "n_confident_wrong": int((confident & wrong).sum()),
    }


def summarise_dataset(truth: pd.DataFrame, evals: dict[str, pd.DataFrame], meta: dict,
                      strata_cols: list[str]) -> dict[str, pd.DataFrame]:
    """Summarise one simulated dataset across all aligners.

    truth   annotated truth table (one row per read)
    evals   {aligner: per-read outcome table from evaluate.evaluate}
    meta    constant columns: dataset, genome, platform, divergence
    """
    accuracy, mapq_rows, strata_rows, kinds = [], [], [], []
    keyed = truth.set_index(["qname", "mate"])
    for aligner, ev in evals.items():
        frame = ev.set_index(["qname", "mate"]).join(keyed[strata_cols], how="left")
        base = {**meta, "aligner": aligner}

        row = {**base, **_counts(frame)}
        for threshold in (0.5, 0.9):  # sensitivity of "correct" to the overlap criterion
            row[f"n_correct_iou{int(threshold * 100)}"] = int((frame["iou"] > threshold).sum())
        row["n_wrong_strand"] = int(((frame["outcome"] == "correct") & ~frame["strand_ok"]).sum())
        accuracy.append(row)

        mapped = frame[frame["outcome"] != "unmapped"]
        per_q = mapped.groupby("mapq", observed=True)["outcome"].agg(
            n_mapped="size", n_wrong=lambda s: int((s == "misplaced").sum())).reset_index()
        for _, r in per_q.iterrows():
            mapq_rows.append({**base, "mapq": int(r["mapq"]), "n_mapped": int(r["n_mapped"]), "n_wrong": int(r["n_wrong"])})

        for kind, count in frame.loc[frame["outcome"] == "misplaced", "misplaced_kind"].value_counts().items():
            kinds.append({**base, "misplaced_kind": kind, "n": int(count)})

        for col in strata_cols:
            for label, group in frame.groupby(col, observed=True):
                strata_rows.append({**base, "stratum_set": col, "stratum": str(label), **_counts(group)})

    return {"accuracy": pd.DataFrame(accuracy), "mapq": pd.DataFrame(mapq_rows),
            "strata": pd.DataFrame(strata_rows), "misplaced_kinds": pd.DataFrame(kinds)}


# --------------------------------------------------------------------------- #
# derived tables used by the report
# --------------------------------------------------------------------------- #
def add_rates(df: pd.DataFrame) -> pd.DataFrame:
    """Counts -> rates with 95 % Wilson intervals."""
    out = df.copy()
    n = out["n"].to_numpy(float)
    for name in ("correct", "misplaced", "unmapped"):
        out[f"rate_{name}"] = out[f"n_{name}"] / np.maximum(n, 1)
    lo, hi = wilson_interval(out["n_misplaced"], out["n"])
    out["rate_misplaced_lo"], out["rate_misplaced_hi"] = lo, hi
    conf = out["n_confident"].to_numpy(float)
    out["rate_confident_wrong"] = np.where(conf > 0, out["n_confident_wrong"] / np.maximum(conf, 1), np.nan)
    lo, hi = wilson_interval(out["n_confident_wrong"], out["n_confident"])
    out["rate_confident_wrong_lo"], out["rate_confident_wrong_hi"] = lo, hi
    return out


def calibration_table(mapq: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Reported MAPQ against the error rate actually observed.

    For each MAPQ value q: the aligner claims P(wrong) = 10^(-q/10); we
    observed n_wrong out of n_mapped. `empirical_mapq` is the observed error
    rate on the Phred scale. With zero observed errors the true rate is not
    zero, only too small to see with n reads, so the interval's upper bound is
    reported and the row is flagged `no_errors_seen`.

    The cumulative columns answer the practical question "if I keep reads with
    MAPQ >= q, what fraction of my reads is that and how many are wrong?".
    """
    grouped = mapq.groupby(group_cols + ["mapq"], observed=True)[["n_mapped", "n_wrong"]].sum().reset_index()
    parts = []
    for _, sub in grouped.groupby(group_cols, observed=True):
        sub = sub.sort_values("mapq", ascending=False).copy()
        sub["cum_mapped"] = sub["n_mapped"].cumsum()
        sub["cum_wrong"] = sub["n_wrong"].cumsum()
        parts.append(sub)
    out = pd.concat(parts, ignore_index=True) if parts else grouped
    out["claimed_error"] = np.power(10.0, -out["mapq"] / 10.0)
    out["observed_error"] = out["n_wrong"] / out["n_mapped"]
    lo, hi = wilson_interval(out["n_wrong"], out["n_mapped"])
    out["observed_error_lo"], out["observed_error_hi"] = lo, hi
    out["no_errors_seen"] = out["n_wrong"] == 0
    out["empirical_mapq"] = phred(np.where(out["no_errors_seen"], hi, out["observed_error"]))
    out["cum_error"] = out["cum_wrong"] / out["cum_mapped"]
    lo, hi = wilson_interval(out["cum_wrong"], out["cum_mapped"])
    out["cum_error_lo"], out["cum_error_hi"] = lo, hi
    return out.sort_values(group_cols + ["mapq"]).reset_index(drop=True)


def headline(accuracy: pd.DataFrame, calibration: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """One row per aligner and platform: the numbers the report opens with."""
    acc = accuracy.groupby(group_cols, observed=True)[
        ["n", "n_correct", "n_misplaced", "n_unmapped", "n_confident", "n_confident_wrong"]].sum().reset_index()
    acc = add_rates(acc)
    acc["claimed_error_at_mapq30"] = 1e-3
    # how many times more errors than promised among alignments reported at MAPQ >= 30
    acc["excess_over_claim"] = acc["rate_confident_wrong"] / 1e-3
    return acc
