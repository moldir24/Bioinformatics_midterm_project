"""Learn a platform error profile from real reads (task 3).

A simulator is only as trustworthy as its error model. Rather than typing in
error rates from a vendor brochure, we measure them: real reads are aligned
to their reference and every aligned base is compared with the reference.

Two model families, because the platforms fail in different ways:

  short reads (Illumina)
      Errors depend on the sequencing cycle and are flagged by the base
      quality. We learn, per mate:
        * the read-length distribution (after adapter trimming)
        * a first-order Markov chain of quality values along the read
          (quality at cycle c given quality at cycle c-1), which reproduces
          both the per-cycle decay and the within-read correlation
        * the mismatch probability given (quality, position in read)
        * which base is called when a base is wrong (substitution matrix)
        * rare indels and the fragment-length distribution

  long reads (ONT, PacBio HiFi)
      Errors are indel-dominated, vary strongly from read to read, and pile
      up in homopolymers. We learn:
        * the joint sample of (read length, per-read error rate)
        * how errors split into substitutions / insertions / deletions
        * indel length distributions
        * how much more likely an indel is inside a homopolymer of length L

Known variant sites of the sequenced individual can be masked so that true
variants are not counted as sequencing errors.

The same code profiles simulated reads, so "does the simulation match?" is
answered by comparing two profiles produced by one identical procedure.
"""
from __future__ import annotations

import json

import numpy as np
import pysam

from .seqio import CODE, homopolymer_run_lengths
from .stats import ks_distance, normalise, total_variation
from .variants import VariantSet, mask_positions

QMAX = 94  # Phred values 0..93 fit the printable FASTQ range
MAX_CYCLES = 400  # longest short read we model
CYCLE_BINS = 4  # read split into quarters for the error-given-quality table
HPMAX = 10  # homopolymer lengths >= HPMAX are pooled
INDEL_LEN_MAX = 50  # longer gaps are structural variants, not sequencing error
PSEUDO = 500.0  # pseudo-count (in bases) when smoothing sparse error cells

_Q_CONSUME = np.zeros(10, bool)
_Q_CONSUME[[0, 1, 4, 7, 8]] = True  # M I S = X advance along the read
_R_CONSUME = np.zeros(10, bool)
_R_CONSUME[[0, 2, 3, 7, 8]] = True  # M D N = X advance along the reference
_ALIGNED = np.zeros(10, bool)
_ALIGNED[[0, 7, 8]] = True  # M = X put a read base opposite a reference base


def walk_cigar(cigartuples):
    """Expand a CIGAR into index arrays, without a Python loop over bases.

    Returns a dict with
        q_idx, r_idx   read offset and reference offset of every aligned pair
        ins_q, ins_r, ins_len   insertions (read offset, reference offset, length)
        del_r, del_len          deletions  (reference offset of first deleted base, length)
    Reference offsets are relative to the alignment start.
    """
    ops = np.fromiter((c[0] for c in cigartuples), dtype=np.int64, count=len(cigartuples))
    lens = np.fromiter((c[1] for c in cigartuples), dtype=np.int64, count=len(cigartuples))
    q_off = np.cumsum(lens * _Q_CONSUME[ops]) - lens * _Q_CONSUME[ops]
    r_off = np.cumsum(lens * _R_CONSUME[ops]) - lens * _R_CONSUME[ops]
    m = _ALIGNED[ops]
    m_len = lens[m]
    total = int(m_len.sum())
    within = np.arange(total) - np.repeat(np.cumsum(m_len) - m_len, m_len)
    ins = ops == 1
    dele = ops == 2
    return {
        "q_idx": np.repeat(q_off[m], m_len) + within,
        "r_idx": np.repeat(r_off[m], m_len) + within,
        "ins_q": q_off[ins], "ins_r": r_off[ins], "ins_len": lens[ins],
        "del_r": r_off[dele], "del_len": lens[dele],
    }


def _usable(read, min_mapq: int) -> bool:
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.mapping_quality < min_mapq or read.query_sequence is None)


def _masked(mask: np.ndarray, ref_start: int, ref_end: int) -> np.ndarray:
    """Masked reference positions that fall inside [ref_start, ref_end)."""
    if mask.size == 0:
        return mask
    return mask[np.searchsorted(mask, ref_start):np.searchsorted(mask, ref_end)]


# --------------------------------------------------------------------------- #
# short reads
# --------------------------------------------------------------------------- #
def learn_short(bam_path: str, fasta_path: str, variants: VariantSet | None = None,
                max_reads: int = 100_000, min_mapq: int = 20, source: str = "") -> dict:
    """Learn the short-read model from a BAM of real (or simulated) reads."""
    fasta = pysam.FastaFile(fasta_path)
    masks: dict[str, np.ndarray] = {}

    qual = [np.full((max_reads, MAX_CYCLES), 255, np.uint8) for _ in (0, 1)]
    lengths = [np.zeros(max_reads, np.int32) for _ in (0, 1)]
    n_kept = [0, 0]
    err_num = np.zeros((2, CYCLE_BINS, QMAX))
    err_den = np.zeros((2, CYCLE_BINS, QMAX))
    cyc_num = np.zeros((2, MAX_CYCLES))
    cyc_den = np.zeros((2, MAX_CYCLES))
    subst = np.zeros((4, 4))
    ins_len = np.zeros(INDEL_LEN_MAX + 1)
    del_len = np.zeros(INDEL_LEN_MAX + 1)
    aligned_bases = 0
    fragments: list[int] = []
    paired = False

    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if not _usable(read, min_mapq):
                continue
            mate = 1 if read.is_read2 else 0
            paired = paired or read.is_paired
            if n_kept[mate] >= max_reads:
                if n_kept[0] >= max_reads and (not paired or n_kept[1] >= max_reads):
                    break
                continue
            length = read.query_length
            if length == 0 or length > MAX_CYCLES or read.query_qualities is None:
                continue
            contig = read.reference_name
            if contig not in masks:
                masks[contig] = mask_positions(variants, contig)
            walk = walk_cigar(read.cigartuples)
            rs = read.reference_start
            ref = CODE[np.frombuffer(fasta.fetch(contig, rs, read.reference_end).encode(), np.uint8)]
            seq = CODE[np.frombuffer(read.query_sequence.encode(), np.uint8)]
            q = np.frombuffer(read.query_qualities.tobytes(), np.uint8)

            q_idx, r_idx = walk["q_idx"], walk["r_idx"]
            ref_base, read_base = ref[r_idx], seq[q_idx]
            keep = (ref_base < 4) & (read_base < 4)
            hidden = _masked(masks[contig], rs, read.reference_end)
            if hidden.size:
                keep &= ~np.isin(r_idx + rs, hidden)
            q_idx, ref_base, read_base = q_idx[keep], ref_base[keep], read_base[keep]
            mism = ref_base != read_base

            # the sequencer's view: cycle 0 is the first base it read
            cycle = (length - 1 - q_idx) if read.is_reverse else q_idx
            cyc_den[mate] += np.bincount(cycle, minlength=MAX_CYCLES)
            cyc_num[mate] += np.bincount(cycle[mism], minlength=MAX_CYCLES)
            cell = (cycle * CYCLE_BINS // length) * QMAX + np.minimum(q[q_idx], QMAX - 1)
            err_den[mate] += np.bincount(cell, minlength=CYCLE_BINS * QMAX).reshape(CYCLE_BINS, QMAX)
            err_num[mate] += np.bincount(cell[mism], minlength=CYCLE_BINS * QMAX).reshape(CYCLE_BINS, QMAX)
            if mism.any():
                true_b, called_b = ref_base[mism], read_base[mism]
                if read.is_reverse:  # complement back to the strand that was sequenced
                    true_b, called_b = 3 - true_b, 3 - called_b
                subst += np.bincount(true_b * 4 + called_b, minlength=16).reshape(4, 4)

            aligned_bases += int(keep.sum())
            for arr, lens, pos in ((ins_len, walk["ins_len"], walk["ins_r"]), (del_len, walk["del_len"], walk["del_r"])):
                if lens.size:
                    ok = lens <= INDEL_LEN_MAX
                    if hidden.size:
                        ok &= ~np.isin(pos + rs, hidden)
                    arr += np.bincount(lens[ok], minlength=INDEL_LEN_MAX + 1)

            q_orig = q[::-1] if read.is_reverse else q
            i = n_kept[mate]
            qual[mate][i, :length] = np.minimum(q_orig, QMAX - 1)
            lengths[mate][i] = length
            n_kept[mate] += 1
            if read.is_proper_pair and read.is_read1 and read.template_length != 0:
                fragments.append(abs(read.template_length))

    if n_kept[0] == 0:
        raise ValueError(f"{bam_path}: no usable alignments (primary, MAPQ >= {min_mapq})")

    prof: dict = {"q_first": [], "q_trans": [], "length_probs": [], "err_by_q": [], "cycle_err": [], "cycle_meanq": []}
    n_mates = 2 if paired and n_kept[1] > 0 else 1
    for mate in range(n_mates):
        n = n_kept[mate]
        qm = qual[mate][:n]
        lm = lengths[mate][:n]
        lmax = int(lm.max())
        prof["length_probs"].append(normalise(np.bincount(lm, minlength=MAX_CYCLES + 1)))
        prof["q_first"].append(normalise(np.bincount(qm[:, 0], minlength=256)[:QMAX]))
        trans = np.zeros((MAX_CYCLES, QMAX, QMAX), np.float32)
        meanq = np.zeros(MAX_CYCLES)
        for c in range(lmax):
            col = qm[:, c]
            present = col != 255
            if present.any():
                meanq[c] = col[present].mean()
            if c == 0:
                continue
            both = present & (qm[:, c - 1] != 255)
            pairs = qm[both, c - 1].astype(np.int64) * QMAX + col[both]
            trans[c] = np.bincount(pairs, minlength=QMAX * QMAX).reshape(QMAX, QMAX)
        prof["q_trans"].append(trans)
        prof["cycle_meanq"].append(meanq)
        # mismatch probability given quality and quarter of the read, shrunk
        # towards the per-quality average where a cell has few observations
        by_q = err_num[mate].sum(0) / np.maximum(err_den[mate].sum(0), 1)
        smoothed = (err_num[mate] + PSEUDO * by_q) / (err_den[mate] + PSEUDO)
        prof["err_by_q"].append(smoothed)
        prof["cycle_err"].append(cyc_num[mate] / np.maximum(cyc_den[mate], 1))

    frag = np.array(fragments, np.int64)
    if frag.size:
        frag = frag[frag <= np.quantile(frag, 0.999)]  # drop chimeric outliers
    n_ins, n_del = ins_len.sum(), del_len.sum()
    meta = {
        "kind": "short", "paired": n_mates == 2, "source": source,
        "reads_used": [int(x) for x in n_kept[:n_mates]], "aligned_bases": int(aligned_bases),
        "mismatch_rate": float(cyc_num.sum() / max(cyc_den.sum(), 1)),
        "ins_rate": float(n_ins / max(aligned_bases, 1)), "del_rate": float(n_del / max(aligned_bases, 1)),
        "mean_read_length": float(np.concatenate([lengths[m][:n_kept[m]] for m in range(n_mates)]).mean()),
        "fragment_median": float(np.median(frag)) if frag.size else None,
        "fragment_max": int(frag.max()) if frag.size else None,
        "variant_sites_masked": variants is not None,
        "min_mapq": min_mapq,
    }
    return {
        "meta": meta,
        "q_first": np.stack(prof["q_first"]), "q_trans": np.stack(prof["q_trans"]),
        "length_probs": np.stack(prof["length_probs"]), "err_by_q": np.stack(prof["err_by_q"]),
        "cycle_err": np.stack(prof["cycle_err"]), "cycle_meanq": np.stack(prof["cycle_meanq"]),
        "subst": subst, "ins_len": ins_len, "del_len": del_len,
        "fragments": frag if frag.size else np.zeros(0, np.int64),
    }


# --------------------------------------------------------------------------- #
# long reads
# --------------------------------------------------------------------------- #
def learn_long(bam_path: str, fasta_path: str, variants: VariantSet | None = None,
               max_reads: int = 20_000, min_mapq: int = 20, source: str = "") -> dict:
    """Learn the long-read model from a BAM of real (or simulated) reads."""
    fasta = pysam.FastaFile(fasta_path)
    masks: dict[str, np.ndarray] = {}
    read_len, read_err, read_q = [], [], []
    sub_b = ins_b = del_b = aligned = 0
    ins_len = np.zeros(INDEL_LEN_MAX + 1)
    del_len = np.zeros(INDEL_LEN_MAX + 1)
    hp_den = np.zeros(HPMAX + 1)
    hp_ins = np.zeros(HPMAX + 1)
    hp_del = np.zeros(HPMAX + 1)
    subst = np.zeros((4, 4))
    n = 0

    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if not _usable(read, min_mapq):
                continue
            contig = read.reference_name
            if contig not in masks:
                masks[contig] = mask_positions(variants, contig)
            rs, re_ = read.reference_start, read.reference_end
            walk = walk_cigar(read.cigartuples)
            ref_ascii = np.frombuffer(fasta.fetch(contig, rs, re_).encode(), np.uint8)
            ref = CODE[ref_ascii]
            seq = CODE[np.frombuffer(read.query_sequence.encode(), np.uint8)]
            hidden = _masked(masks[contig], rs, re_)

            q_idx, r_idx = walk["q_idx"], walk["r_idx"]
            ref_base, read_base = ref[r_idx], seq[q_idx]
            keep = (ref_base < 4) & (read_base < 4)
            if hidden.size:
                keep &= ~np.isin(r_idx + rs, hidden)
            mism = keep & (ref_base != read_base)

            def events(lens, pos):
                ok = lens <= INDEL_LEN_MAX
                if hidden.size:
                    ok &= ~np.isin(pos + rs, hidden)
                return lens[ok], pos[ok]

            i_len, i_pos = events(walk["ins_len"], walk["ins_r"])
            d_len, d_pos = events(walk["del_len"], walk["del_r"])
            n_m, n_x = int(keep.sum()), int(mism.sum())
            n_i, n_d = int(i_len.sum()), int(d_len.sum())
            columns = n_m + n_i + n_d
            if columns < 200:
                continue

            read_len.append(read.infer_read_length())
            read_err.append((n_x + n_i + n_d) / columns)
            quals = read.query_qualities
            read_q.append(float(np.mean(quals)) if quals is not None and len(quals) else np.nan)
            sub_b += n_x
            ins_b += n_i
            del_b += n_d
            aligned += n_m
            ins_len += np.bincount(i_len, minlength=INDEL_LEN_MAX + 1)
            del_len += np.bincount(d_len, minlength=INDEL_LEN_MAX + 1)
            if n_x:
                tb, cb = ref_base[mism], read_base[mism]
                if read.is_reverse:
                    tb, cb = 3 - tb, 3 - cb
                subst += np.bincount(tb * 4 + cb, minlength=16).reshape(4, 4)

            # homopolymer context of every indel, against the context of all bases
            hp = np.minimum(homopolymer_run_lengths(ref_ascii), HPMAX)
            hp_den += np.bincount(hp, minlength=HPMAX + 1)
            last = hp.size - 1
            if d_pos.size:
                hp_del += np.bincount(hp[np.minimum(d_pos, last)], minlength=HPMAX + 1)
            if i_pos.size:
                left = hp[np.clip(i_pos - 1, 0, last)]
                right = hp[np.minimum(i_pos, last)]
                hp_ins += np.bincount(np.maximum(left, right), minlength=HPMAX + 1)

            n += 1
            if n >= max_reads:
                break

    if n == 0:
        raise ValueError(f"{bam_path}: no usable alignments (primary, MAPQ >= {min_mapq})")

    err_bases = max(sub_b + ins_b + del_b, 1)
    overall_ins = ins_len.sum() / max(hp_den.sum(), 1)
    overall_del = del_len.sum() / max(hp_den.sum(), 1)
    # relative indel-event rate at each homopolymer length (1.0 = genome average),
    # smoothed so that rare long homopolymers do not produce wild factors
    hp_ins_factor = ((hp_ins + 5.0 * 1.0) / (hp_den * overall_ins + 5.0)) if overall_ins > 0 else np.ones(HPMAX + 1)
    hp_del_factor = ((hp_del + 5.0 * 1.0) / (hp_den * overall_del + 5.0)) if overall_del > 0 else np.ones(HPMAX + 1)
    lengths = np.array(read_len, np.int64)
    errors = np.array(read_err, float)
    order = np.sort(lengths)[::-1]
    n50 = int(order[np.searchsorted(np.cumsum(order), order.sum() / 2)])
    meta = {
        "kind": "long", "paired": False, "source": source, "reads_used": int(n),
        "aligned_bases": int(aligned),
        "mismatch_rate": float(sub_b / max(aligned, 1)),
        "ins_rate": float(ins_b / max(aligned, 1)), "del_rate": float(del_b / max(aligned, 1)),
        "frac_sub": float(sub_b / err_bases), "frac_ins": float(ins_b / err_bases), "frac_del": float(del_b / err_bases),
        "mean_read_length": float(lengths.mean()), "read_n50": n50,
        "mean_error": float(errors.mean()), "median_error": float(np.median(errors)),
        "mean_ins_len": float((ins_len * np.arange(INDEL_LEN_MAX + 1)).sum() / max(ins_len.sum(), 1)),
        "mean_del_len": float((del_len * np.arange(INDEL_LEN_MAX + 1)).sum() / max(del_len.sum(), 1)),
        "variant_sites_masked": variants is not None, "min_mapq": min_mapq,
    }
    return {
        "meta": meta, "lengths": lengths, "errors": errors, "mean_quals": np.array(read_q, float),
        "ins_len": ins_len, "del_len": del_len, "subst": subst,
        "hp_ins_factor": hp_ins_factor, "hp_del_factor": hp_del_factor, "hp_den": hp_den,
    }


# --------------------------------------------------------------------------- #
# closing the loop
# --------------------------------------------------------------------------- #
def refine_profile(real: dict, pilot: dict, limit: float = 4.0) -> dict:
    """Correct the simulator's injected error rates with a pilot run.

    An error rate measured through an alignment is not exactly the rate at
    which errors were made: the aligner turns an adjacent insertion and
    deletion into a mismatch, slides indels inside homopolymers, clips noisy
    read ends. So injecting the measured rates does not reproduce them.

    We therefore simulate a pilot set with the measured rates, measure it by
    the same procedure, and scale each error type by real / pilot. The scale
    factors are stored under `inject` and applied by the simulator; the
    measured (target) rates stay untouched in the profile.
    """
    out = dict(real)
    meta = dict(real["meta"])
    previous = meta.get("inject", {"sub": 1.0, "ins": 1.0, "del": 1.0})
    inject = {}
    for key, name in (("mismatch_rate", "sub"), ("ins_rate", "ins"), ("del_rate", "del")):
        target, observed = real["meta"][key], pilot["meta"][key]
        factor = target / observed if target > 0 and observed > 0 else 1.0
        inject[name] = float(np.clip(previous[name] * factor, 1.0 / limit, limit))
    meta["inject"] = inject
    meta["pilot_observed"] = {k: pilot["meta"][k] for k in ("mismatch_rate", "ins_rate", "del_rate")}
    out["meta"] = meta
    return out


# --------------------------------------------------------------------------- #
# storage
# --------------------------------------------------------------------------- #
def save_profile(path: str, prof: dict) -> None:
    """Arrays go into one compressed .npz; scalar facts into a JSON string
    inside it and, for people, into a sidecar .json next to it."""
    arrays = {k: v for k, v in prof.items() if k != "meta"}
    np.savez_compressed(path, meta=json.dumps(prof["meta"]), **arrays)
    sidecar = str(path)[:-4] + ".json" if str(path).endswith(".npz") else str(path) + ".json"
    with open(sidecar, "w") as fh:
        json.dump(prof["meta"], fh, indent=2)


def load_profile(path: str) -> dict:
    with np.load(path, allow_pickle=False) as data:
        prof = {k: data[k] for k in data.files if k != "meta"}
        prof["meta"] = json.loads(str(data["meta"]))
    return prof


# --------------------------------------------------------------------------- #
# does the simulation match?  (the "calibrated" part of task 3)
# --------------------------------------------------------------------------- #
def _rel(real: float, sim: float) -> float:
    return float((sim - real) / real) if real else float("nan")


def compare_profiles(real: dict, sim: dict, tolerance: float = 0.25) -> list[dict]:
    """Compare a profile measured on real reads with one measured, by the same
    procedure, on simulated reads.

    Scalar rates are compared by relative difference; distributions by a
    distance in [0, 1] (Kolmogorov-Smirnov or total variation). `ok` marks
    rows inside the tolerance; the remaining rows are where the simulator is
    not yet realistic and must be reported as a limitation.
    """
    rows: list[dict] = []

    def scalar(name, unit):
        r, s = real["meta"].get(name), sim["meta"].get(name)
        if r is None or s is None:
            return
        rel = _rel(r, s)
        rows.append({"metric": name, "unit": unit, "real": r, "simulated": s, "difference": rel,
                     "kind": "relative difference", "ok": bool(abs(rel) <= tolerance)})

    def dist(name, value):
        rows.append({"metric": name, "unit": "distance 0-1", "real": 0.0, "simulated": value, "difference": value,
                     "kind": "distribution distance", "ok": bool(value <= tolerance)})

    for name, unit in (("mismatch_rate", "per aligned base"), ("ins_rate", "per aligned base"),
                       ("del_rate", "per aligned base"), ("mean_read_length", "bp")):
        scalar(name, unit)

    if real["meta"]["kind"] == "short":
        scalar("fragment_median", "bp")
        mates = min(real["cycle_err"].shape[0], sim["cycle_err"].shape[0])
        for mate in range(mates):
            n = int(min(np.flatnonzero(real["length_probs"][mate]).max(), np.flatnonzero(sim["length_probs"][mate]).max()))
            r_err, s_err = real["cycle_err"][mate][:n], sim["cycle_err"][mate][:n]
            r_q, s_q = real["cycle_meanq"][mate][:n], sim["cycle_meanq"][mate][:n]
            corr = float(np.corrcoef(r_err, s_err)[0, 1]) if n > 2 and r_err.std() > 0 and s_err.std() > 0 else float("nan")
            rows.append({"metric": f"per-cycle error curve, mate {mate + 1}", "unit": "Pearson r", "real": 1.0,
                         "simulated": corr, "difference": 1.0 - corr, "kind": "curve correlation",
                         "ok": bool(corr >= 1.0 - tolerance)})
            gap = float(np.abs(r_q - s_q).max())
            rows.append({"metric": f"per-cycle mean quality, mate {mate + 1}", "unit": "max Phred gap", "real": 0.0,
                         "simulated": gap, "difference": gap, "kind": "curve gap", "ok": bool(gap <= 2.0)})
            dist(f"quality distribution, mate {mate + 1}",
                 total_variation(real["q_trans"][mate].sum((0, 1)), sim["q_trans"][mate].sum((0, 1))))
        dist("substitution spectrum", total_variation(real["subst"].ravel(), sim["subst"].ravel()))
        if real["fragments"].size and sim["fragments"].size:
            dist("fragment length distribution", ks_distance(real["fragments"], sim["fragments"]))
    else:
        for name, unit in (("mean_error", "per aligned column"), ("frac_sub", "share of errors"),
                           ("frac_ins", "share of errors"), ("frac_del", "share of errors"), ("read_n50", "bp")):
            scalar(name, unit)
        dist("read length distribution", ks_distance(real["lengths"], sim["lengths"]))
        dist("per-read error distribution", ks_distance(real["errors"], sim["errors"]))
        dist("insertion length distribution", total_variation(real["ins_len"], sim["ins_len"]))
        dist("deletion length distribution", total_variation(real["del_len"], sim["del_len"]))
        for name in ("hp_ins_factor", "hp_del_factor"):
            # compare the enrichment at homopolymers of length >= 4, where it matters
            r = float(np.average(real[name][4:], weights=real["hp_den"][4:] + 1))
            s = float(np.average(sim[name][4:], weights=real["hp_den"][4:] + 1))
            rel = _rel(r, s)
            rows.append({"metric": f"{name} (homopolymer >= 4)", "unit": "x genome average", "real": r, "simulated": s,
                         "difference": rel, "kind": "relative difference", "ok": bool(abs(rel) <= tolerance)})
    return rows
