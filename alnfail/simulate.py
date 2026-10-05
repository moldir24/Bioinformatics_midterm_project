"""Simulate reads whose true origin is known by construction (tasks 4 and 5).

Every simulated read is built in four layers, and the reference coordinate of
every base is carried through all of them:

    1. pick a window of the reference                     (uniform over the genome)
    2. turn it into the sample's sequence                 (truth-set variants of one haplotype)
    3. optionally diverge it from the reference           (substitutions + small indels)
    4. sequence it with the platform's measured errors    (profile learned in profile.py)

Because layer 1 fixes where the read comes from and layers 2-4 only edit
bases, the "true" alignment interval is exact: it is the span of reference
positions of the bases that ended up in the read. That interval is written
to a truth table (Parquet) keyed by read name and mate.

Short and long reads use different error layers because the platforms fail
differently (see profile.py).
"""
from __future__ import annotations

import zlib

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .profile import CYCLE_BINS, HPMAX, INDEL_LEN_MAX, MAX_CYCLES, QMAX
from .seqio import BASES, CODE, fastq_writer, homopolymer_run_lengths, to_array
from .variants import VariantSet

BATCH = 20_000


def make_rng(seed: int, dataset_id: str) -> np.random.Generator:
    """One reproducible random stream per dataset: same seed and id, same reads."""
    return np.random.default_rng([int(seed), zlib.crc32(dataset_id.encode())])


class Genome:
    """Random access to a reference through its .fai index: only the windows we
    sample are ever read from disk, never the whole genome."""

    def __init__(self, fasta_path: str, contigs: list[str] | None = None):
        self.fasta = pysam.FastaFile(fasta_path)
        names = list(self.fasta.references)
        if contigs:
            missing = [c for c in contigs if c not in names]
            if missing:
                raise ValueError(f"contigs not in {fasta_path}: {missing}")
            names = list(contigs)
        self.names = names
        self.lengths = np.array([self.fasta.get_reference_length(c) for c in names], np.int64)
        self.cum_weights = np.cumsum(self.lengths / self.lengths.sum())

    def fetch(self, contig: str, start: int, end: int) -> np.ndarray:
        return to_array(self.fasta.fetch(contig, start, end))

    def window(self, rng: np.random.Generator, size: int, max_tries: int = 200):
        """A random window with no N bases. Windows are drawn uniformly over
        the sequence (contig chosen in proportion to its length), so assembly
        gaps are skipped but every other base is equally likely."""
        for _ in range(max_tries):
            ci = min(int(np.searchsorted(self.cum_weights, rng.random(), side="right")), self.lengths.size - 1)
            length = int(self.lengths[ci])
            size_here = min(size, length)
            start = int(rng.integers(0, length - size_here + 1))
            seq = self.fetch(self.names[ci], start, start + size_here)
            if not (CODE[seq] == 4).any():
                return self.names[ci], start, seq
        raise RuntimeError("could not find an N-free window; is the reference mostly gaps?")


# --------------------------------------------------------------------------- #
# edit layers shared by both platforms
# --------------------------------------------------------------------------- #
def _apply_indels(codes, rpos, del_start, del_len, ins_after, ins_len, ins_bases):
    """Delete runs and insert bases in one pass, keeping `rpos` in step.

    del_start[i], del_len[i]   delete del_len bases starting at index del_start
    ins_after[j], ins_len[j]   insert ins_len bases after index ins_after
    ins_bases                  concatenated inserted base codes
    Inserted bases inherit the reference position of the base they follow.
    """
    n = codes.size
    keep = np.ones(n, bool)
    if del_start.size:
        idx = np.repeat(del_start, del_len) + (np.arange(int(del_len.sum())) - np.repeat(np.cumsum(del_len) - del_len, del_len))
        keep[idx[idx < n]] = False
    out_codes, out_pos = codes[keep], rpos[keep]
    if ins_after.size:
        alive = keep[ins_after]  # an insertion anchored on a deleted base is dropped
        if alive.any():
            new_index = np.cumsum(keep)  # index *after* each surviving base
            where = np.repeat(new_index[ins_after[alive]], ins_len[alive])
            bases = ins_bases[np.repeat(alive, ins_len)]
            anchors = np.repeat(rpos[ins_after[alive]], ins_len[alive])
            out_codes = np.insert(out_codes, where, bases)
            out_pos = np.insert(out_pos, where, anchors)
    return out_codes, out_pos


def apply_divergence(codes, rpos, rate: float, rng: np.random.Generator, indel_share: float = 0.1):
    """Make the sampled sequence differ from the reference at `rate` per base.

    This models reads from a genome that is related to, but not the same as,
    the reference: another strain, another individual, another species. 90 %
    of the events are substitutions and 10 % short indels, roughly the ratio
    seen between closely related genomes. Returns (codes, rpos, n_events).
    """
    if rate <= 0:
        return codes, rpos, 0
    n = codes.size
    codes = codes.copy()
    sub = rng.random(n) < rate * (1.0 - indel_share)
    n_sub = int(sub.sum())
    if n_sub:
        codes[sub] = (codes[sub] + rng.integers(1, 4, n_sub)) % 4  # always a different base
    n_indel = int(rng.binomial(n, rate * indel_share))
    if n_indel:
        pos = np.sort(rng.choice(n, size=min(n_indel, n), replace=False))
        is_del = rng.random(pos.size) < 0.5
        size = rng.geometric(0.6, pos.size)  # mostly 1-2 bp
        ins_len = size[~is_del]
        codes, rpos = _apply_indels(codes, rpos, pos[is_del], size[is_del], pos[~is_del], ins_len,
                                    rng.integers(0, 4, int(ins_len.sum())).astype(codes.dtype))
    return codes, rpos, n_sub + n_indel


def _template(genome: Genome, rng, size: int, variants: VariantSet | None, divergence: float):
    """Layers 1-3: a window of the reference turned into sample sequence.

    Returns (contig, base codes, reference position per base, haplotype,
    variants applied, divergence events).
    """
    contig, start, seq = genome.window(rng, size)
    positions = None
    hap, n_var = 0, 0
    if variants is not None:
        hap = int(rng.integers(0, 2))
        seq, positions, n_var = variants.apply(contig, hap, start, seq)
    if positions is None:
        positions = np.arange(start, start + seq.size, dtype=np.int64)
    codes = CODE[seq]
    codes, positions, n_div = apply_divergence(codes, positions, divergence, rng)
    return contig, codes, positions, hap, n_var, n_div


def _cdf_rows(matrix: np.ndarray) -> np.ndarray:
    """Row-normalised cumulative sums, last column forced to 1."""
    total = matrix.sum(-1, keepdims=True)
    cdf = np.cumsum(matrix / np.where(total > 0, total, 1.0), -1)
    cdf[..., -1] = 1.0
    return cdf


def _draw(cdf_rows: np.ndarray, rng) -> np.ndarray:
    """One categorical draw per row of a table of cumulative probabilities."""
    return (rng.random(cdf_rows.shape[0])[:, None] > cdf_rows[:, :-1]).sum(1)


class _TruthWriter:
    """Collects per-read truth and writes one Parquet file."""

    COLUMNS = ("qname", "mate", "contig", "start", "end", "strand", "length", "gc",
               "n_var", "n_div", "n_err", "hap")

    def __init__(self):
        self.cols = {c: [] for c in self.COLUMNS}

    def add(self, **values):
        for key in self.COLUMNS:
            self.cols[key].append(values[key])

    def write(self, path: str):
        table = pa.table({
            "qname": pa.array(np.concatenate(self.cols["qname"]), pa.string()),
            "mate": pa.array(np.concatenate(self.cols["mate"]), pa.int8()),
            "contig": pa.array(np.concatenate(self.cols["contig"]), pa.string()).dictionary_encode(),
            "start": pa.array(np.concatenate(self.cols["start"]), pa.int64()),
            "end": pa.array(np.concatenate(self.cols["end"]), pa.int64()),
            "strand": pa.array(np.concatenate(self.cols["strand"]), pa.string()).dictionary_encode(),
            "length": pa.array(np.concatenate(self.cols["length"]), pa.int32()),
            "gc": pa.array(np.concatenate(self.cols["gc"]), pa.float32()),
            "n_var": pa.array(np.concatenate(self.cols["n_var"]), pa.int32()),
            "n_div": pa.array(np.concatenate(self.cols["n_div"]), pa.int32()),
            "n_err": pa.array(np.concatenate(self.cols["n_err"]), pa.int32()),
            "hap": pa.array(np.concatenate(self.cols["hap"]), pa.int8()),
        })
        pq.write_table(table, path, compression="zstd")


# --------------------------------------------------------------------------- #
# short reads
# --------------------------------------------------------------------------- #
class ShortReadModel:
    """Sampler built from a short-read profile."""

    def __init__(self, prof: dict):
        self.meta = prof["meta"]
        self.n_mates = prof["q_first"].shape[0]
        self.fragments = prof["fragments"]
        if self.n_mates == 2 and self.fragments.size == 0:
            raise ValueError("paired-end profile without any properly paired read: cannot learn fragment lengths")
        self.length_cdf = [np.cumsum(p) for p in prof["length_probs"]]
        inject = self.meta.get("inject", {"sub": 1.0, "ins": 1.0, "del": 1.0})  # pilot-run correction
        self.err = np.clip(prof["err_by_q"] * inject["sub"], 0.0, 0.75)  # [mate, quarter of read, quality]
        self.subst_cdf = _cdf_rows(prof["subst"] * (1 - np.eye(4)) + 1e-12 * (1 - np.eye(4)))
        self.ins_rate, self.del_rate = self.meta["ins_rate"] * inject["ins"], self.meta["del_rate"] * inject["del"]
        self.ins_len_cdf = np.cumsum(prof["ins_len"] / max(prof["ins_len"].sum(), 1))
        self.del_len_cdf = np.cumsum(prof["del_len"] / max(prof["del_len"].sum(), 1))
        # quality chain, compressed to the quality values the instrument actually emits
        self.qvals, self.first_cdf, self.trans_cdf = [], [], []
        for mate in range(self.n_mates):
            trans = prof["q_trans"][mate].astype(float)
            seen = np.flatnonzero((trans.sum((0, 1)) + trans.sum((0, 2)) + prof["q_first"][mate]) > 0)
            self.qvals.append(seen)
            self.first_cdf.append(np.cumsum(prof["q_first"][mate][seen] / prof["q_first"][mate][seen].sum()))
            t = trans[:, seen][:, :, seen]
            last_good = None
            cdfs = np.zeros_like(t)
            for c in range(1, MAX_CYCLES):
                if t[c].sum() == 0:  # beyond the longest observed read: reuse the last cycle
                    cdfs[c] = cdfs[last_good] if last_good is not None else np.cumsum(np.eye(seen.size), 1)
                    continue
                marginal = t[c].sum(0)
                rows = np.where(t[c].sum(1, keepdims=True) > 0, t[c], marginal[None, :])  # unseen previous value
                cdfs[c] = _cdf_rows(rows)
                last_good = c
            self.trans_cdf.append(cdfs)

    def draw_lengths(self, mate: int, n: int, rng) -> np.ndarray:
        return np.searchsorted(self.length_cdf[mate], rng.random(n), side="right").clip(1, MAX_CYCLES)

    def draw_fragments(self, n: int, rng) -> np.ndarray:
        return rng.choice(self.fragments, size=n)

    def qualities(self, mate: int, n: int, width: int, rng) -> np.ndarray:
        """Phred values [reads, cycles] from the per-cycle Markov chain."""
        idx = np.zeros((n, width), np.int64)
        idx[:, 0] = np.searchsorted(self.first_cdf[mate], rng.random(n), side="right").clip(0, self.qvals[mate].size - 1)
        for c in range(1, width):
            idx[:, c] = _draw(self.trans_cdf[mate][c][idx[:, c - 1]], rng)
        return self.qvals[mate][idx].astype(np.uint8)

    def sequence(self, mate: int, codes: np.ndarray, rpos: np.ndarray, lengths: np.ndarray, rng):
        """Apply sequencing errors to a batch of reads (padded matrices).

        Returns (codes, rpos, quals, lengths, n_errors).
        """
        n, width = codes.shape
        quals = self.qualities(mate, n, width, rng)
        cols = np.arange(width)[None, :]
        valid = cols < lengths[:, None]
        quarter = np.minimum(cols * CYCLE_BINS // lengths[:, None], CYCLE_BINS - 1)
        p_err = self.err[mate][quarter, np.minimum(quals, QMAX - 1)]
        mism = (rng.random((n, width)) < p_err) & valid
        n_err = mism.sum(1).astype(np.int32)
        if mism.any():
            true_base = codes[mism]
            codes[mism] = _draw(self.subst_cdf[true_base], rng)

        # indels are rare on this platform (about one per ten thousand bases),
        # so they are applied read by read instead of as matrices
        for rate, len_cdf, is_del in ((self.del_rate, self.del_len_cdf, True), (self.ins_rate, self.ins_len_cdf, False)):
            total = int(lengths.sum())
            k = int(rng.poisson(total * rate)) if rate > 0 else 0
            if k == 0:
                continue
            rows = rng.choice(n, size=k, p=lengths / total)  # longer reads are hit more often
            sizes = np.searchsorted(len_cdf, rng.random(k), side="right").clip(1, INDEL_LEN_MAX)
            for row, size in zip(rows, sizes):
                length = int(lengths[row])
                if length < 2 * size + 2:
                    continue
                at = int(rng.integers(1, length - size))
                c, r, q = codes[row, :length], rpos[row, :length], quals[row, :length]
                if is_del:  # the read simply becomes `size` bases shorter
                    c, r, q = np.delete(c, slice(at, at + size)), np.delete(r, slice(at, at + size)), np.delete(q, slice(at, at + size))
                else:  # inserted bases push the tail off the end of the fixed-length read
                    c = np.insert(c, at, rng.integers(0, 4, size))[:length]
                    r = np.insert(r, at, np.repeat(r[at - 1], size))[:length]
                    q = np.insert(q, at, np.repeat(q[at], size))[:length]
                new_len = c.size
                codes[row, :new_len], rpos[row, :new_len], quals[row, :new_len] = c, r, q
                lengths[row] = new_len
                n_err[row] += 1
        return codes, rpos, quals, lengths, n_err


def _mate(codes: np.ndarray, pos: np.ndarray, length: int, reverse: bool, adapter: np.ndarray | None):
    """Cut one mate out of a fragment.

    The forward mate reads the fragment from its left end; the reverse mate
    reads the reverse complement from its right end. If the fragment is
    shorter than the read and an adapter is given, the read runs through into
    the adapter and then into poly-G (what a two-colour instrument reports
    when there is nothing left to read). Such bases get reference position -1.
    """
    size = codes.size
    k = min(length, size)
    if reverse:
        c, p = 3 - codes[size - k:][::-1], pos[size - k:][::-1]
    else:
        c, p = codes[:k], pos[:k]
    if k < length and adapter is not None:
        pad = np.concatenate((adapter, np.full(length, 2, np.uint8)))[: length - k]
        c, p = np.concatenate((c, pad)), np.concatenate((p, np.full(pad.size, -1, np.int64)))
    return c, p


def simulate_short(fasta_path: str, prof: dict, n: int, out_r1: str, out_r2: str | None, truth_path: str,
                   seed: int, dataset_id: str, contigs: list[str] | None = None,
                   variants: VariantSet | None = None, divergence: float = 0.0,
                   adapters: tuple[str, str] | None = None) -> dict:
    """Simulate n short fragments (pairs if the profile is paired-end).

    `adapters` switches on adapter read-through for fragments shorter than
    the read. The benchmark reads are simulated without it, because real
    reads are adapter-trimmed before alignment; it is used to create
    realistic raw input for testing the QC step.
    """
    adapter_codes = [CODE[to_array(a)] for a in adapters] if adapters else [None, None]
    rng = make_rng(seed, dataset_id)
    genome = Genome(fasta_path, contigs)
    model = ShortReadModel(prof)
    paired = model.n_mates == 2 and out_r2 is not None
    truth = _TruthWriter()
    made = 0

    def run(handle1, handle2):
        nonlocal made
        while made < n:
            b = min(BATCH, n - made)
            len1 = model.draw_lengths(0, b, rng)
            len2 = model.draw_lengths(1, b, rng) if paired else None
            frag = model.draw_fragments(b, rng) if paired else len1.copy()
            width = int(max(len1.max(), len2.max() if paired else 0))
            # per mate: base codes and the reference position of every base, padded to `width`
            mats = [(np.zeros((b, width), np.uint8), np.zeros((b, width), np.int64)) for _ in range(2 if paired else 1)]
            contig_col = np.empty(b, dtype=object)
            forward = rng.random(b) < 0.5  # strand of mate 1
            hap = np.zeros(b, np.int8)
            n_var = np.zeros(b, np.int32)
            n_div = np.zeros(b, np.int32)
            for i in range(b):
                longest = int(max(len1[i], len2[i] if paired else 0))
                # without adapters a fragment is at least as long as its reads
                want = int(frag[i]) if adapters else int(max(frag[i], longest))
                # extra bases so that deletions (variants or divergence) cannot leave the template short
                contig, codes, pos, hap[i], n_var[i], n_div[i] = _template(
                    genome, rng, want + 20 + int(0.05 * want), variants, divergence)
                codes, pos = codes[:want], pos[:want]
                contig_col[i] = contig
                # mate 1 reads the fragment from one end, mate 2 from the other, on opposite strands
                c, p = _mate(codes, pos, int(len1[i]), not forward[i], adapter_codes[0])
                len1[i] = c.size
                mats[0][0][i, :c.size], mats[0][1][i, :c.size] = c, p
                if paired:
                    c, p = _mate(codes, pos, int(len2[i]), bool(forward[i]), adapter_codes[1])
                    len2[i] = c.size
                    mats[1][0][i, :c.size], mats[1][1][i, :c.size] = c, p

            names = np.array([f"{dataset_id}:{made + i}" for i in range(b)], dtype=object)
            for mate, (handle, lens) in enumerate(((handle1, len1), (handle2, len2))[: 2 if paired else 1]):
                codes, rpos = mats[mate]
                cols = np.arange(width)[None, :]
                valid = (cols < lens[:, None]) & (rpos >= 0)
                gc = (((codes == 1) | (codes == 2)) & valid).sum(1) / np.maximum(valid.sum(1), 1)  # before errors
                codes, rpos, quals, lens, n_err = model.sequence(mate, codes, rpos, lens, rng)
                valid = (cols < lens[:, None]) & (rpos >= 0)  # genomic bases only, not adapter
                start = np.where(valid, rpos, np.iinfo(np.int64).max).min(1)
                end = np.where(valid, rpos, -1).max(1) + 1
                ascii_bases = BASES[codes]
                ascii_quals = (quals + 33).astype(np.uint8)
                lines = []
                for i in range(b):
                    k = lens[i]
                    lines.append(f"@{names[i]}\n{ascii_bases[i, :k].tobytes().decode()}\n+\n{ascii_quals[i, :k].tobytes().decode()}\n")
                handle.write("".join(lines))
                mate_forward = forward if mate == 0 else ~forward
                truth.add(qname=names, mate=np.full(b, mate + 1 if paired else 0, np.int8), contig=contig_col,
                          start=start, end=end, strand=np.where(mate_forward, "+", "-").astype(object),
                          length=lens.astype(np.int32), gc=gc.astype(np.float32), n_var=n_var, n_div=n_div,
                          n_err=n_err, hap=hap)
            made += b

    with fastq_writer(out_r1) as h1:
        if paired:
            with fastq_writer(out_r2) as h2:
                run(h1, h2)
        else:
            run(h1, None)
    truth.write(truth_path)
    return {"dataset": dataset_id, "fragments": made, "paired": paired, "divergence": divergence,
            "variants": variants is not None}


# --------------------------------------------------------------------------- #
# long reads
# --------------------------------------------------------------------------- #
class LongReadModel:
    """Sampler built from a long-read profile."""

    def __init__(self, prof: dict):
        self.meta = m = prof["meta"]
        self.lengths, self.errors = prof["lengths"], prof["errors"]
        inject = m.get("inject", {"sub": 1.0, "ins": 1.0, "del": 1.0})  # pilot-run correction
        self.frac_sub = m["frac_sub"] * inject["sub"]
        self.frac_ins = m["frac_ins"] * inject["ins"]
        self.frac_del = m["frac_del"] * inject["del"]
        self.mean_ins, self.mean_del = max(m["mean_ins_len"], 1.0), max(m["mean_del_len"], 1.0)
        self.ins_len_cdf = np.cumsum(prof["ins_len"] / max(prof["ins_len"].sum(), 1))
        self.del_len_cdf = np.cumsum(prof["del_len"] / max(prof["del_len"].sum(), 1))
        self.hp_ins, self.hp_del = prof["hp_ins_factor"], prof["hp_del_factor"]
        self.subst_cdf = _cdf_rows(prof["subst"] * (1 - np.eye(4)) + 1e-12 * (1 - np.eye(4)))

    def draw(self, rng) -> tuple[int, float]:
        """Length and error rate are drawn together from one real read, which
        keeps their correlation (long reads are not uniformly accurate)."""
        i = int(rng.integers(0, self.lengths.size))
        return int(self.lengths[i]), float(self.errors[i])

    def sequence(self, codes: np.ndarray, rpos: np.ndarray, err: float, rng):
        """Apply indel-dominated, homopolymer-aware errors to one template."""
        n = codes.size
        codes = codes.copy()
        hp = np.minimum(homopolymer_run_lengths(codes), HPMAX)
        hp_next = np.concatenate((hp[1:], hp[-1:]))
        p_sub = err * self.frac_sub
        p_del = err * self.frac_del / self.mean_del * self.hp_del[hp]
        p_ins = err * self.frac_ins / self.mean_ins * self.hp_ins[np.maximum(hp, hp_next)]
        u = rng.random((3, n))
        sub = u[0] < p_sub
        n_sub = int(sub.sum())
        if n_sub:
            codes[sub] = _draw(self.subst_cdf[codes[sub]], rng)
        del_start = np.flatnonzero(u[1] < p_del)
        ins_after = np.flatnonzero(u[2] < p_ins)
        del_len = np.searchsorted(self.del_len_cdf, rng.random(del_start.size), side="right").clip(1, INDEL_LEN_MAX)
        ins_len = np.searchsorted(self.ins_len_cdf, rng.random(ins_after.size), side="right").clip(1, INDEL_LEN_MAX)
        # half of the inserted bases repeat the base before them (homopolymer
        # extension, the typical nanopore/polymerase slip); the rest are random
        total_ins = int(ins_len.sum())
        ins_bases = rng.integers(0, 4, total_ins).astype(codes.dtype)
        if total_ins:
            copy = rng.random(total_ins) < 0.5
            ins_bases[copy] = np.repeat(codes[ins_after], ins_len)[copy]
        out_codes, out_pos = _apply_indels(codes, rpos, del_start, del_len, ins_after, ins_len, ins_bases)
        return out_codes, out_pos, n_sub + int(del_start.size) + int(ins_after.size)


def simulate_long(fasta_path: str, prof: dict, n: int, out_fastq: str, truth_path: str, seed: int,
                  dataset_id: str, contigs: list[str] | None = None, variants: VariantSet | None = None,
                  divergence: float = 0.0) -> dict:
    """Simulate n long reads."""
    rng = make_rng(seed, dataset_id)
    genome = Genome(fasta_path, contigs)
    model = LongReadModel(prof)
    truth = _TruthWriter()
    cols = {k: [] for k in _TruthWriter.COLUMNS}
    with fastq_writer(out_fastq) as handle:
        for i in range(n):
            length, err = model.draw(rng)
            contig, codes, pos, hap, n_var, n_div = _template(
                genome, rng, int(length * 1.15) + 100, variants, divergence)
            gc = float(((codes == 1) | (codes == 2)).mean())
            codes, pos, n_err = model.sequence(codes, pos, err, rng)
            codes, pos = codes[:length], pos[:length]
            forward = rng.random() < 0.5
            if not forward:
                codes = 3 - codes[::-1]
            q = int(np.clip(round(-10 * np.log10(max(err, 1e-6))), 1, 60))
            name = f"{dataset_id}:{i}"
            handle.write(f"@{name}\n{BASES[codes].tobytes().decode()}\n+\n{chr(q + 33) * codes.size}\n")
            for key, value in (("qname", name), ("mate", 0), ("contig", contig), ("start", int(pos.min())),
                               ("end", int(pos.max()) + 1), ("strand", "+" if forward else "-"),
                               ("length", codes.size), ("gc", gc), ("n_var", n_var), ("n_div", n_div),
                               ("n_err", n_err), ("hap", hap)):
                cols[key].append(value)
    truth.add(**{k: np.array(v, dtype=object if k in ("qname", "contig", "strand") else None) for k, v in cols.items()})
    truth.write(truth_path)
    return {"dataset": dataset_id, "reads": n, "paired": False, "divergence": divergence,
            "variants": variants is not None}
