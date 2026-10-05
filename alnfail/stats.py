"""Small statistical helpers used across the benchmark.

Everything here works on numpy arrays so it can be applied to whole
columns at once.
"""
from __future__ import annotations

import numpy as np

Z95 = 1.959963984540054  # two-sided 95 % normal quantile


def wilson_interval(k, n, z: float = Z95):
    """Wilson score interval for a binomial proportion k/n.

    Preferred over the normal approximation because our interesting
    proportions are tiny (error rates near 1e-3 and below), where the normal
    interval collapses to zero width or goes negative.

    Returns (low, high) arrays. Where n == 0 both bounds are NaN.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n
        denom = 1.0 + z * z / n
        centre = (p + z * z / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
        low = np.clip(centre - half, 0.0, 1.0)
        high = np.clip(centre + half, 0.0, 1.0)
    # exact end points: with k == 0 the lower bound is 0, with k == n the upper is 1
    low = np.where(k == 0, 0.0, low)
    high = np.where(k == n, 1.0, high)
    low = np.where(n == 0, np.nan, low)
    high = np.where(n == 0, np.nan, high)
    return low, high


def phred(p, floor: float = 1e-10):
    """Probability of error -> Phred scale, Q = -10 * log10(p)."""
    p = np.asarray(p, dtype=float)
    return -10.0 * np.log10(np.clip(p, floor, 1.0))


def unphred(q):
    """Phred scale -> probability of error, p = 10 ** (-Q / 10)."""
    return np.power(10.0, -np.asarray(q, dtype=float) / 10.0)


def ks_distance(a, b) -> float:
    """Two-sample Kolmogorov-Smirnov distance: the largest gap between the two
    empirical cumulative distributions. 0 means identical, 1 means disjoint."""
    a = np.sort(np.asarray(a, dtype=float))
    b = np.sort(np.asarray(b, dtype=float))
    if a.size == 0 or b.size == 0:
        return float("nan")
    grid = np.concatenate([a, b])
    cdf_a = np.searchsorted(a, grid, side="right") / a.size
    cdf_b = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def total_variation(p, q) -> float:
    """Total variation distance between two discrete distributions given as
    count (or probability) vectors over the same categories. 0 = identical,
    1 = no overlap. Read it as "the share of probability mass you would have to
    move to turn one histogram into the other"."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.sum() == 0 or q.sum() == 0:
        return float("nan")
    return float(0.5 * np.abs(p / p.sum() - q / q.sum()).sum())


def spearman(x, y) -> float:
    """Spearman rank correlation without scipy (average ranks for ties)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return float("nan")

    def rank(v):
        order = np.argsort(v, kind="mergesort")
        ranks = np.empty(v.size, dtype=float)
        ranks[order] = np.arange(v.size, dtype=float)
        # average the ranks of tied values
        _, inv, counts = np.unique(v, return_inverse=True, return_counts=True)
        sums = np.bincount(inv, weights=ranks)
        return (sums / counts)[inv]

    rx, ry = rank(x), rank(y)
    sx, sy = rx.std(), ry.std()
    if sx == 0 or sy == 0:
        return float("nan")
    return float(((rx - rx.mean()) * (ry - ry.mean())).mean() / (sx * sy))


def normalise(counts):
    """Counts -> probabilities; an all-zero vector becomes uniform."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.full(counts.shape, 1.0 / max(counts.size, 1))
    return counts / total
