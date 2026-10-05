"""Scoring, calibration arithmetic and the consensus used for real data."""
import numpy as np
import pandas as pd

from alnfail.consensus import consensus
from alnfail.evaluate import evaluate, interval_iou
from alnfail.mapq_theory import ideal_mapq, two_copy_experiment
from alnfail.profile import walk_cigar
from alnfail.qc import mean_quality, n50
from alnfail.stats import ks_distance, phred, spearman, total_variation, wilson_interval
from alnfail.summarise import calibration_table


def _alignment(qname, contig, start, end, mapq, mapped=True, strand="+", mate=0):
    return {"qname": qname, "mate": mate, "mapped": mapped, "contig": contig, "start": start, "end": end,
            "strand": strand, "mapq": mapq, "nm": 0, "read_len": end - start, "clipped": 0, "proper": False,
            "split": False}


def test_three_outcomes():
    truth = pd.DataFrame({"qname": list("abcde"), "mate": 0, "contig": "chr1", "start": 1000, "end": 1150,
                          "strand": "+"})
    aln = pd.DataFrame([
        _alignment("a", "chr1", 1002, 1152, 60),  # right place
        _alignment("b", "chr1", 5000, 5150, 60),  # same contig, elsewhere
        _alignment("c", "chr2", 1000, 1150, 3),  # same coordinates on another contig is still wrong
        _alignment("d", "", -1, -1, 0, mapped=False),  # reported as unmapped
    ])  # "e" is missing from the output altogether, as NGMLR does for unplaced reads
    out = evaluate(truth, aln).set_index("qname")
    assert list(out["outcome"]) == ["correct", "misplaced", "misplaced", "unmapped", "unmapped"]
    assert out.loc["b", "misplaced_kind"] == "same_contig_elsewhere"
    assert out.loc["c", "misplaced_kind"] == "other_contig"


def test_overlap_threshold_is_intersection_over_union():
    assert interval_iou(np.array([0]), np.array([100]), np.array([0]), np.array([100]))[0] == 1.0
    assert interval_iou(np.array([0]), np.array([100]), np.array([100]), np.array([200]))[0] == 0.0
    assert abs(interval_iou(np.array([0]), np.array([100]), np.array([50]), np.array([150]))[0] - 1 / 3) < 1e-9


def test_wilson_interval_behaves_at_the_edges():
    lo, hi = wilson_interval(0, 1000)
    assert lo == 0 and 0.003 < hi < 0.005  # "no errors in 1,000" still allows about 4 per 1,000
    lo, hi = wilson_interval(10, 1000)
    assert lo < 0.01 < hi
    lo, hi = wilson_interval(0, 0)
    assert np.isnan(lo) and np.isnan(hi)


def test_phred_scale():
    assert abs(phred(0.001) - 30) < 1e-9
    assert abs(phred(0.5) - 3.0103) < 1e-3  # two equally good places: MAPQ 3


def test_calibration_table_accumulates_from_high_mapq_down():
    mapq = pd.DataFrame({"aligner": "x", "mapq": [60, 30, 0], "n_mapped": [1000, 100, 50], "n_wrong": [0, 1, 25]})
    table = calibration_table(mapq, ["aligner"]).set_index("mapq")
    assert table.loc[60, "cum_mapped"] == 1000 and table.loc[30, "cum_mapped"] == 1100 and table.loc[0, "cum_mapped"] == 1150
    assert table.loc[0, "cum_wrong"] == 26
    assert abs(table.loc[30, "empirical_mapq"] - 20) < 1e-6  # 1 wrong in 100 is Phred 20, not the claimed 30
    assert table.loc[60, "no_errors_seen"] and table.loc[60, "empirical_mapq"] < 30  # only a lower bound


def test_ideal_mapq_values():
    assert abs(ideal_mapq(0, 0.01) - 3.0103) < 1e-3  # identical copies
    assert 24 < ideal_mapq(1, 0.01) < 25.5  # one distinguishing mismatch at 1 % error
    assert ideal_mapq(2, 0.01) > 49


def test_exact_posterior_is_calibrated():
    table = two_copy_experiment(error_rate=0.05, reads_per_setting=200_000, differences=(0, 1, 2, 3))
    pooled = table.groupby("mapq")[["n_reads", "n_wrong"]].sum()
    pooled = pooled[pooled["n_wrong"] >= 30]
    observed = -10 * np.log10(pooled["n_wrong"] / pooled["n_reads"])
    assert len(pooled) >= 2
    assert np.all(np.abs(observed - pooled.index) < 2.0)  # earned MAPQ equals claimed MAPQ


def test_cigar_walk_matches_a_hand_alignment():
    # 2S 3M 1I 2M 2D 1M : read offsets 0-1 clipped, then aligned / inserted / aligned
    walk = walk_cigar([(4, 2), (0, 3), (1, 1), (0, 2), (2, 2), (0, 1)])
    assert list(walk["q_idx"]) == [2, 3, 4, 6, 7, 8]
    assert list(walk["r_idx"]) == [0, 1, 2, 3, 4, 7]
    assert list(walk["ins_q"]) == [5] and list(walk["ins_r"]) == [3] and list(walk["ins_len"]) == [1]
    assert list(walk["del_r"]) == [5] and list(walk["del_len"]) == [2]


def test_consensus_needs_a_strict_majority():
    def table(rows):
        return pd.DataFrame([_alignment(*r) for r in rows])

    tables = {
        "a": table([("r1", "chr1", 100, 250, 60), ("r2", "chr1", 100, 250, 60), ("r3", "chr1", 100, 250, 60)]),
        "b": table([("r1", "chr1", 102, 252, 60), ("r2", "chr1", 100, 250, 60), ("r3", "chr2", 100, 250, 60)]),
        "c": table([("r1", "chr1", 100, 250, 40), ("r2", "chr1", 9000, 9150, 55), ("r3", "chr3", 100, 250, 60)]),
        "d": table([("r1", "chr2", 100, 250, 50), ("r2", "chr1", 9000, 9150, 0), ("r3", "chr4", 100, 250, 60)]),
    }
    cons = consensus(tables).set_index("qname")
    assert cons.loc["r1", "has_consensus"] and cons.loc["r1", "status_d"] == "disagree" and cons.loc["r1", "status_a"] == "agree"
    assert not cons.loc["r2", "has_consensus"]  # two against two
    assert not cons.loc["r3", "has_consensus"]  # everyone somewhere else


def test_small_statistics():
    assert ks_distance([1, 2, 3], [1, 2, 3]) == 0.0
    assert ks_distance([1, 2, 3], [10, 11, 12]) == 1.0
    assert total_variation([1, 0], [0, 1]) == 1.0
    assert abs(spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1) < 1e-9
    assert n50(np.array([2, 2, 2, 3, 3, 4, 8, 8])) == 8
    # one terrible base drags the mean quality down far more than the mean of the Phred numbers suggests
    assert mean_quality("I" * 99 + "!") < 21


def test_reads_without_stored_qualities_are_not_filtered_as_low_quality(tmp_path):
    import gzip

    from alnfail.qc import filter_long

    fastq = tmp_path / "r.fq.gz"
    with gzip.open(fastq, "wt") as fh:
        fh.write("@no_quality\n" + "ACGT" * 500 + "\n+\n" + "!" * 2000 + "\n")  # qualities unknown
        fh.write("@bad\n" + "ACGT" * 500 + "\n+\n" + "$" * 2000 + "\n")  # Q3 throughout
        fh.write("@good\n" + "ACGT" * 500 + "\n+\n" + "I" * 2000 + "\n")
        fh.write("@short\n" + "ACGT" * 50 + "\n+\n" + "I" * 200 + "\n")
    rows = filter_long(str(fastq), str(tmp_path / "c.fq.gz"), str(tmp_path / "h.fq.gz"), 1000, 7.0, "t")
    kept = gzip.open(tmp_path / "c.fq.gz", "rt").read() + gzip.open(tmp_path / "h.fq.gz", "rt").read()
    assert "@no_quality" in kept and "@good" in kept
    assert "@bad" not in kept and "@short" not in kept
    assert any(r["metric"] == "reads_without_quality_fraction" for r in rows)
