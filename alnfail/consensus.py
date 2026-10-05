"""Does the simulation tell the truth about real data? (task 8)

Real reads have no known origin, so accuracy cannot be measured on them
directly. Two things can be measured, and both are also measured on the
simulated reads so that like is compared with like:

  1. Observables. What an aligner reports regardless of truth: the share of
     reads it maps, the distribution of MAPQ, edit distance, soft clipping.
     If simulated reads were realistic, an aligner would treat them the way
     it treats real ones and these distributions would coincide.

  2. Disagreement with the consensus. For each read, the placement supported
     by a strict majority of aligners is taken as a stand-in for truth; an
     aligner that puts the read elsewhere "disagrees". Disagreement is not
     error: aligners sharing an algorithm can be wrong together. So the
     stand-in is first checked on simulated reads, where truth is known,
     before it is trusted on real ones.

The Genome in a Bottle truth set enters twice: its variants define the
variant-density strata, and its high-confidence regions define the
`giab_confident` stratum, on real and simulated reads alike. GIAB certifies
variants, not read positions, so it cannot give per-read truth; the report
must say so.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .evaluate import interval_iou
from .stats import spearman, total_variation
from .summarise import MAPQ_CONFIDENT

MAPQ_BINS = np.arange(0, 62)


def consensus(tables: dict[str, pd.DataFrame], min_overlap: float = 0.1) -> pd.DataFrame:
    """Per-read consensus placement across aligners.

    tables   {aligner: table from samtable.bam_to_table}

    Returns one row per read with the consensus interval (when a strict
    majority of aligners agree on one) and each aligner's status:
    'agree', 'disagree' or 'unmapped'.
    """
    aligners = list(tables)
    keys = pd.concat([t[["qname", "mate"]] for t in tables.values()]).drop_duplicates().reset_index(drop=True)
    n = len(keys)
    k = len(aligners)
    contig_names = sorted({c for t in tables.values() for c in t["contig"].astype(str).unique()})
    contig_code = {c: i for i, c in enumerate(contig_names)}

    mapped = np.zeros((k, n), bool)
    contig = np.full((k, n), -1, np.int64)
    start = np.zeros((k, n), np.int64)
    end = np.zeros((k, n), np.int64)
    mapq = np.zeros((k, n), np.int16)
    for i, aligner in enumerate(aligners):
        m = keys.merge(tables[aligner], on=["qname", "mate"], how="left", validate="one_to_one")
        mapped[i] = m["mapped"].eq(True).to_numpy()  # missing read -> not mapped
        contig[i] = m["contig"].astype(str).map(contig_code).fillna(-1).to_numpy(np.int64)
        start[i] = m["start"].fillna(0).to_numpy(np.int64)
        end[i] = m["end"].fillna(0).to_numpy(np.int64)
        mapq[i] = m["mapq"].fillna(0).to_numpy(np.int16)

    # agree[i, j]: aligners i and j put the read in the same place
    agree = np.zeros((k, k, n), bool)
    for i in range(k):
        for j in range(i, k):
            same = mapped[i] & mapped[j] & (contig[i] == contig[j])
            both = same & (interval_iou(start[i], end[i], start[j], end[j]) > min_overlap)
            agree[i, j] = agree[j, i] = both
    support = agree.sum(1)  # for each aligner's placement: how many aligners back it (itself included)
    best = support.argmax(0)
    rows = np.arange(n)
    has_consensus = support[best, rows] >= (k // 2 + 1)

    out = keys.copy()
    out["has_consensus"] = has_consensus
    out["support"] = support[best, rows]
    out["contig"] = np.where(has_consensus, np.array(contig_names + [""], dtype=object)[contig[best, rows]], "")
    out["start"] = np.where(has_consensus, start[best, rows], 0)
    out["end"] = np.where(has_consensus, end[best, rows], 0)
    for i, aligner in enumerate(aligners):
        status = np.where(~mapped[i], "unmapped", np.where(agree[i, best, rows], "agree", "disagree"))
        out[f"status_{aligner}"] = status
        out[f"mapq_{aligner}"] = mapq[i]
    return out


def consensus_summary(cons: pd.DataFrame, aligners: list[str], meta: dict, strata_cols: list[str],
                      truth: pd.DataFrame | None = None, min_overlap: float = 0.1) -> pd.DataFrame:
    """Disagreement with the consensus, per aligner and stratum.

    cons    consensus() output, annotated with stratum columns (only reads
            with a consensus have a position and hence a stratum)
    truth   for simulated data: the truth table, to measure how often the
            consensus itself is right
    """
    with_cons = cons[cons["has_consensus"]].copy()
    if truth is not None:
        t = truth[["qname", "mate", "contig", "start", "end"]].rename(
            columns={"contig": "t_contig", "start": "t_start", "end": "t_end"})
        with_cons = with_cons.merge(t, on=["qname", "mate"], how="left")
        same = (with_cons["contig"].astype(str) == with_cons["t_contig"].astype(str)).to_numpy()
        iou = interval_iou(with_cons["start"].to_numpy(), with_cons["end"].to_numpy(),
                           with_cons["t_start"].to_numpy(), with_cons["t_end"].to_numpy())
        with_cons["consensus_correct"] = same & (iou > min_overlap)

    rows = []
    for aligner in aligners:
        status_all = cons[f"status_{aligner}"].to_numpy()
        base = {**meta, "aligner": aligner}
        rows.append({**base, "stratum_set": "all", "stratum": "all", "n_reads": int(len(cons)),
                     "n_no_consensus": int((~cons["has_consensus"]).sum()),
                     "n_unmapped_total": int((status_all == "unmapped").sum()),
                     **_disagreement(with_cons, aligner)})
        for col in strata_cols:
            if col not in with_cons.columns:
                continue
            for label, group in with_cons.groupby(col, observed=True):
                rows.append({**base, "stratum_set": col, "stratum": str(label), "n_reads": int(len(group)),
                             "n_no_consensus": 0, "n_unmapped_total": int((group[f"status_{aligner}"] == "unmapped").sum()),
                             **_disagreement(group, aligner)})
    return pd.DataFrame(rows)


def _disagreement(frame: pd.DataFrame, aligner: str) -> dict:
    status = frame[f"status_{aligner}"].to_numpy()
    mapq = frame[f"mapq_{aligner}"].to_numpy()
    disagree = status == "disagree"
    out = {
        "n_consensus": int(len(frame)),
        "n_agree": int((status == "agree").sum()),
        "n_disagree": int(disagree.sum()),
        "n_disagree_confident": int((disagree & (mapq >= MAPQ_CONFIDENT)).sum()),
        "n_unmapped": int((status == "unmapped").sum()),
    }
    if "consensus_correct" in frame.columns:
        out["n_consensus_correct"] = int(frame["consensus_correct"].sum())
    return out


def observables(tables: dict[str, pd.DataFrame], meta: dict) -> pd.DataFrame:
    """Truth-free quantities per aligner, in tidy form (metric, key, value)."""
    rows = []
    for aligner, t in tables.items():
        base = {**meta, "aligner": aligner}
        n = len(t)
        mapped = t[t["mapped"]]
        span = np.maximum((mapped["end"] - mapped["start"]).to_numpy(float), 1.0)
        has_nm = (mapped["nm"] >= 0).to_numpy()
        clip = mapped["clipped"].to_numpy(float) / np.maximum(mapped["read_len"].to_numpy(float), 1.0)
        scalars = {
            "reads": n,
            "mapped_fraction": len(mapped) / max(n, 1),
            "confident_fraction": float((mapped["mapq"] >= MAPQ_CONFIDENT).sum()) / max(n, 1),
            "mapq0_fraction": float((mapped["mapq"] == 0).sum()) / max(n, 1),
            "edit_distance_per_base": float((mapped["nm"].to_numpy(float)[has_nm] / span[has_nm]).mean()) if has_nm.any() else np.nan,
            "clipped_fraction_mean": float(clip.mean()) if len(mapped) else np.nan,
            "clipped_reads_fraction": float((clip > 0.05).mean()) if len(mapped) else np.nan,
            "split_fraction": float(mapped["split"].mean()) if len(mapped) else np.nan,
        }
        for metric, value in scalars.items():
            rows.append({**base, "metric": metric, "key": "", "value": value})
        hist = np.bincount(np.clip(mapped["mapq"].to_numpy(np.int64), 0, 60), minlength=61)
        for q in np.flatnonzero(hist):
            rows.append({**base, "metric": "mapq_count", "key": str(int(q)), "value": int(hist[q])})
    return pd.DataFrame(rows)


def provider_agreement(tables: dict[str, pd.DataFrame], provider: pd.DataFrame, meta: dict,
                       min_overlap: float = 0.1) -> pd.DataFrame:
    """How often each aligner reproduces the data provider's own placement.

    The provider (Genome in a Bottle) aligned the same reads against the whole
    genome with a different pipeline. Agreement with it is a second,
    independent stand-in for truth on real data.
    """
    rows = []
    prov = provider.rename(columns={"contig": "p_contig", "start": "p_start", "end": "p_end", "mapq": "p_mapq"})
    for aligner, t in tables.items():
        m = t.merge(prov, on=["qname", "mate"], how="inner")
        mapped = m["mapped"].to_numpy()
        same = mapped & (m["contig"].astype(str) == m["p_contig"].astype(str)).to_numpy()
        iou = interval_iou(m["start"].to_numpy(), m["end"].to_numpy(), m["p_start"].to_numpy(), m["p_end"].to_numpy())
        agree = same & (iou > min_overlap)
        confident = mapped & (m["mapq"].to_numpy() >= MAPQ_CONFIDENT)
        rows.append({**meta, "aligner": aligner, "n_compared": int(len(m)), "n_mapped": int(mapped.sum()),
                     "n_agree_provider": int(agree.sum()),
                     "n_confident": int(confident.sum()), "n_confident_disagree_provider": int((confident & ~agree).sum())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# simulation versus real
# --------------------------------------------------------------------------- #
def transfer_observables(obs: pd.DataFrame, pairs: list[tuple[str, str]], tolerance: float = 0.05) -> pd.DataFrame:
    """Compare truth-free observables between matched simulated and real datasets.

    pairs       [(simulated dataset, real dataset)], same genome and platform
    tolerance   absolute difference below which a share of reads "transfers"
                (0.05 = five percentage points). The MAPQ distribution may
                differ by a total variation distance of twice that; the edit
                distance per base by 25 % of its simulated value.
    """
    rows = []
    for sim_id, real_id in pairs:
        sim, real = obs[obs["dataset"] == sim_id], obs[obs["dataset"] == real_id]
        for aligner in sorted(set(sim["aligner"]) & set(real["aligner"])):
            s, r = sim[sim["aligner"] == aligner], real[real["aligner"] == aligner]
            for metric in ("mapped_fraction", "confident_fraction", "mapq0_fraction", "edit_distance_per_base",
                           "clipped_reads_fraction", "split_fraction"):
                sv = s.loc[s["metric"] == metric, "value"]
                rv = r.loc[r["metric"] == metric, "value"]
                if sv.empty or rv.empty:
                    continue
                sim_value, real_value = float(sv.iloc[0]), float(rv.iloc[0])
                diff = real_value - sim_value
                if metric == "edit_distance_per_base":
                    # a rate, not a share of reads: judged relative to its own size
                    ok = abs(diff) <= 0.25 * max(sim_value, 1e-9)
                else:
                    ok = abs(diff) <= tolerance
                rows.append({"simulated": sim_id, "real": real_id, "aligner": aligner, "metric": metric,
                             "simulated_value": sim_value, "real_value": real_value,
                             "difference": diff, "transfers": bool(ok)})
            hist = []
            for part in (s, r):
                h = np.zeros(61)
                counts = part[part["metric"] == "mapq_count"]
                h[counts["key"].astype(int).to_numpy()] = counts["value"].to_numpy(float)
                hist.append(h)
            tv = total_variation(hist[0], hist[1])
            rows.append({"simulated": sim_id, "real": real_id, "aligner": aligner,
                         "metric": "mapq_distribution_distance", "simulated_value": 0.0, "real_value": tv,
                         "difference": tv, "transfers": bool(tv <= 2 * tolerance)})
    return pd.DataFrame(rows)


def transfer_strata(sim_truth: pd.DataFrame, cons: pd.DataFrame, pairs: list[tuple[str, str]],
                    min_reads: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Is the map of hard regions the same in simulation and in real data?

    sim_truth   stratified truth-based table for simulated data (summarise.strata)
    cons        consensus_summary rows for simulated and real datasets

    For each matched pair and each (aligner, stratum) with at least
    `min_reads` reads on both sides, three rates are lined up:
        sim_error          truth-based misplacement rate in simulation
        sim_disagreement   consensus-based disagreement rate in simulation
        real_disagreement  consensus-based disagreement rate on real reads
    Spearman correlations then answer two questions:
        sim_error  vs sim_disagreement    is disagreement a usable proxy for error?
        sim_disagreement vs real_disagreement   do the same strata fail on real data?
    """
    joined, summary = [], []
    for sim_id, real_id in pairs:
        t = sim_truth[sim_truth["dataset"] == sim_id]
        s = cons[cons["dataset"] == sim_id]
        r = cons[cons["dataset"] == real_id]
        key = ["aligner", "stratum_set", "stratum"]
        m = (t[key + ["n", "n_misplaced"]].rename(columns={"n": "sim_n"})
             .merge(s[key + ["n_consensus", "n_disagree"]].rename(columns={"n_consensus": "sim_n_consensus",
                                                                           "n_disagree": "sim_n_disagree"}), on=key)
             .merge(r[key + ["n_consensus", "n_disagree"]].rename(columns={"n_consensus": "real_n_consensus",
                                                                           "n_disagree": "real_n_disagree"}), on=key))
        m = m[(m["sim_n"] >= min_reads) & (m["sim_n_consensus"] >= min_reads) & (m["real_n_consensus"] >= min_reads)].copy()
        if m.empty:
            continue
        m["sim_error"] = m["n_misplaced"] / m["sim_n"]
        m["sim_disagreement"] = m["sim_n_disagree"] / m["sim_n_consensus"]
        m["real_disagreement"] = m["real_n_disagree"] / m["real_n_consensus"]
        m.insert(0, "real", real_id)
        m.insert(0, "simulated", sim_id)
        joined.append(m)
        for aligner, g in m.groupby("aligner"):
            summary.append({
                "simulated": sim_id, "real": real_id, "aligner": aligner, "strata_compared": int(len(g)),
                "proxy_validity_spearman": spearman(g["sim_error"], g["sim_disagreement"]),
                "transfer_spearman": spearman(g["sim_disagreement"], g["real_disagreement"]),
                "median_ratio_real_over_sim": float(np.nanmedian(
                    (g["real_disagreement"] + 1e-6) / (g["sim_disagreement"] + 1e-6))),
            })
    return (pd.concat(joined, ignore_index=True) if joined else pd.DataFrame(),
            pd.DataFrame(summary))
