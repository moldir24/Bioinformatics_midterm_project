"""Score every simulated read: correct, misplaced or unmapped (tasks 4-6).

Three outcomes, because an aligner can fail in two different ways:

    unmapped    the aligner declined to answer. Costs sensitivity, but the
                failure is visible and harmless downstream.
    misplaced   the aligner answered and was wrong. This is the dangerous
                one: the read becomes false evidence somewhere else.
    correct     the reported position overlaps the true origin.

A read counts as correct when it is on the right contig and its aligned
reference interval overlaps the true interval by more than `min_overlap` of
their union (intersection over union). The default 0.1 follows the criterion
of minimap2's own evaluation script (paftools.js mapeval), so our numbers
are comparable with published ones. It is deliberately tolerant: an aligner
that clips a noisy read end still found the right locus. The stricter
thresholds are kept as extra columns so the report can show how much the
conclusion depends on that choice.

Reads that are absent from an aligner's output are unmapped: some tools
(NGMLR) simply do not print reads they cannot place.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

OUTCOMES = ("correct", "misplaced", "unmapped")


def interval_iou(start_a, end_a, start_b, end_b) -> np.ndarray:
    """Intersection over union of half-open intervals, element-wise."""
    inter = np.clip(np.minimum(end_a, end_b) - np.maximum(start_a, start_b), 0, None)
    union = np.maximum(end_a, end_b) - np.minimum(start_a, start_b)
    return np.where(union > 0, inter / np.maximum(union, 1), 0.0)


def evaluate(truth: pd.DataFrame, alignments: pd.DataFrame, min_overlap: float = 0.1) -> pd.DataFrame:
    """Join truth with one aligner's primary alignments and score each read.

    truth        qname, mate, contig, start, end, strand, ...
    alignments   the table written by samtable.bam_to_table
    """
    aln = alignments.rename(columns={"contig": "aln_contig", "start": "aln_start", "end": "aln_end",
                                     "strand": "aln_strand"})
    merged = truth[["qname", "mate", "contig", "start", "end", "strand"]].merge(
        aln, on=["qname", "mate"], how="left", validate="one_to_one")
    mapped = merged["mapped"].eq(True).to_numpy()  # read absent from the output -> not mapped
    same_contig = mapped & (merged["contig"].astype(str) == merged["aln_contig"].astype(str)).to_numpy()
    aln_start = merged["aln_start"].fillna(-1).to_numpy(np.int64)
    aln_end = merged["aln_end"].fillna(-1).to_numpy(np.int64)
    true_start = merged["start"].to_numpy(np.int64)
    true_end = merged["end"].to_numpy(np.int64)
    iou = np.where(same_contig, interval_iou(true_start, true_end, aln_start, aln_end), 0.0)

    outcome = np.where(~mapped, "unmapped", np.where(iou > min_overlap, "correct", "misplaced"))
    # where did a misplaced read go? another chromosome, far away on the same one, or nearly right
    kind = np.full(mapped.size, "", dtype=object)
    wrong = outcome == "misplaced"
    kind[wrong & ~same_contig] = "other_contig"
    kind[wrong & same_contig & (iou > 0)] = "same_locus_low_overlap"
    kind[wrong & same_contig & (iou == 0)] = "same_contig_elsewhere"

    out = pd.DataFrame({
        "qname": merged["qname"].to_numpy(),
        "mate": merged["mate"].to_numpy(),
        "outcome": pd.Categorical(outcome, categories=OUTCOMES),
        "mapq": merged["mapq"].fillna(0).to_numpy(np.int16),
        "iou": iou.astype(np.float32),
        "strand_ok": mapped & (merged["strand"].astype(str) == merged["aln_strand"].astype(str)).to_numpy(),
        "misplaced_kind": kind,
        "nm": merged["nm"].fillna(-1).to_numpy(np.int32),
        "clipped": merged["clipped"].fillna(0).to_numpy(np.int32),
        "split": merged["split"].eq(True).to_numpy(),
        "aln_contig": merged["aln_contig"].fillna("").astype(str).to_numpy(),
        "aln_start": aln_start,
        "aln_end": aln_end,
    })
    return out
