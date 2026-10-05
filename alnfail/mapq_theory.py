"""What a mapping quality is supposed to mean (task 1).

MAPQ is a probability written on the Phred scale:

    MAPQ = -10 * log10( P(the reported position is wrong) )

so MAPQ 10 claims a 1-in-10 chance of being wrong, 20 claims 1 in 100,
30 claims 1 in 1,000 and 60 claims 1 in a million.

The probability is a Bayesian posterior. If a read r could have come from
candidate positions u1, u2, ... with likelihoods P(r | u), and every position
is equally likely beforehand, then

    P(u_best is right | r) = P(r | u_best) / sum over all u of P(r | u)

A read is hard to place not when it has errors, but when a *second* position
explains it almost as well as the best one.

This module makes that concrete with the smallest possible genome: one locus
and one paralogous copy that differs from it at d positions. Reads are drawn
from the true locus with a known per-base error rate, placed at whichever
copy explains them better, and given the exact posterior as their MAPQ. By
construction this ideal aligner is perfectly calibrated, which is the
reference line real aligners are measured against in task 6.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .stats import phred


def ideal_mapq(delta_mismatches, error_rate: float, n_equal_copies: int = 1):
    """Exact MAPQ for a read whose best hit has `delta_mismatches` fewer
    mismatches than each of `n_equal_copies` competing positions.

    One extra mismatch multiplies the likelihood by (e/3) / (1 - e): the base
    must have been misread (probability e) into that particular wrong base
    (one of three). The best hit therefore wins the posterior by that factor
    to the power delta.
    """
    delta = np.asarray(delta_mismatches, dtype=float)
    ratio = (error_rate / 3.0) / (1.0 - error_rate)  # likelihood penalty per extra mismatch
    p_wrong = n_equal_copies * ratio ** delta / (1.0 + n_equal_copies * ratio ** delta)
    return phred(p_wrong)


def two_copy_experiment(read_length: int = 100, error_rate: float = 0.01, differences=(0, 1, 2, 3, 4, 6),
                        reads_per_setting: int = 200_000, seed: int = 26) -> pd.DataFrame:
    """Monte Carlo check that the exact posterior is calibrated.

    For each number d of differences between the two copies, simulate reads
    from copy A. A read has an error at each of the d distinguishing sites
    with probability e, and with probability 1/3 that error happens to
    produce exactly copy B's base. Mismatches elsewhere hit both copies
    equally and cancel out of the posterior.

    The read is assigned to the copy with fewer mismatches (ties are broken
    at random) and labelled with the exact posterior MAPQ. We then count how
    often reads carrying each MAPQ value were actually wrong.
    """
    rng = np.random.default_rng(seed)
    ratio = (error_rate / 3.0) / (1.0 - error_rate)
    rows = []
    for d in differences:
        n = reads_per_setting
        if d == 0:
            mism_a = np.zeros(n, int)
            mism_b = np.zeros(n, int)
        else:
            u = rng.random((n, d))
            to_b = u < error_rate / 3.0  # misread into copy B's base: matches B, mismatches A
            to_other = (u >= error_rate / 3.0) & (u < error_rate)  # misread into a third base: mismatches both
            mism_a = (to_b | to_other).sum(1)
            mism_b = (~to_b).sum(1)  # every distinguishing site not converted to B's base
        delta = mism_b - mism_a  # > 0: copy A (the truth) fits better
        choose_a = (delta > 0) | ((delta == 0) & (rng.random(n) < 0.5))
        wrong = ~choose_a
        p_wrong = np.where(delta == 0, 0.5, ratio ** np.abs(delta) / (1.0 + ratio ** np.abs(delta)))
        mapq = np.round(phred(p_wrong)).astype(int)
        for q in np.unique(mapq):
            sel = mapq == q
            rows.append({"differences_between_copies": d, "read_length": read_length, "error_rate": error_rate,
                         "mapq": int(q), "n_reads": int(sel.sum()), "n_wrong": int(wrong[sel].sum())})
    out = pd.DataFrame(rows)
    out["observed_error"] = out["n_wrong"] / out["n_reads"]
    out["claimed_error"] = np.power(10.0, -out["mapq"] / 10.0)
    return out


def reference_table() -> pd.DataFrame:
    """The numbers worth knowing by heart at the defence."""
    rows = []
    for copies, label in ((1, "one equally good alternative"), (2, "two"), (4, "four"), (9, "nine")):
        p = copies / (copies + 1.0)
        rows.append({"situation": f"identical repeat, {label}", "p_wrong": p, "ideal_mapq": float(phred(p))})
    for error in (0.001, 0.01, 0.1):
        for delta in (1, 2, 3):
            q = float(ideal_mapq(delta, error))
            rows.append({"situation": f"best hit ahead by {delta} mismatch(es), base error {error:g}",
                         "p_wrong": float(10 ** (-q / 10)), "ideal_mapq": q})
    return pd.DataFrame(rows)
