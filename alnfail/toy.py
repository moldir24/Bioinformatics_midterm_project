"""A miniature test world for verifying the pipeline in minutes.

The handbook asks for "a small test dataset so that [the reproduction]
command can be verified quickly". The real inputs are a human genome and
300x of reads, so the test world is a 350 kb toy genome that has, in
miniature, every feature the real analysis stratifies on:

  * interspersed repeats of different ages (SINE, LINE, LTR, DNA)
  * a satellite array, simple repeats, low-complexity sequence
  * segmental duplications at 99.5 %, 96 % and 92 % identity
  * GC-poor and GC-rich blocks
  * a diploid "individual" with SNVs, indels and one hypervariable region

together with annotation files in exactly the formats of the real sources
(RepeatMasker .out, UCSC genomicSuperDups, GFF3, VCF, an NCBI-style assembly
report) and deliberately inconsistent contig names (chrA / A / TOY_000001.1),
so that the coordinate and identifier handling is exercised too.

IMPORTANT: the "real" reads of the test world are synthetic stand-ins,
produced by hidden instrument models defined below. They exist only so that
the code path can be checked end to end without a network connection. No
number derived from the test world is a result.
"""
from __future__ import annotations

import gzip
import os

import numpy as np

from .profile import CYCLE_BINS, HPMAX, INDEL_LEN_MAX, MAX_CYCLES, QMAX
from .seqio import BASES

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # Illumina TruSeq read-through, read 1
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # Illumina TruSeq read-through, read 2


def _random_seq(rng, n: int, gc: float = 0.41) -> np.ndarray:
    p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
    return rng.choice(4, size=n, p=p).astype(np.uint8)


def _mutate(rng, codes: np.ndarray, divergence: float) -> np.ndarray:
    out = codes.copy()
    hit = rng.random(out.size) < divergence
    out[hit] = (out[hit] + rng.integers(1, 4, int(hit.sum()))) % 4
    return out


def build_genome(seed: int = 26) -> dict:
    """Build the toy genome and every annotation that describes it."""
    rng = np.random.default_rng(seed)
    sizes = {"chrA": 220_000, "chrB": 130_000}
    contigs = {}
    for name, size in sizes.items():  # background with 10 kb blocks of varying GC
        blocks = [_random_seq(rng, 10_000, gc) for gc in rng.choice([0.24, 0.34, 0.41, 0.41, 0.5, 0.62, 0.74], size // 10_000)]
        contigs[name] = np.concatenate(blocks)

    occupied = {name: np.zeros(size, bool) for name, size in sizes.items()}
    repeats = []  # (contig, start, end, name, class/family, perc_div)

    def place(contig, codes, start=None):
        size = codes.size
        for _ in range(5000):
            s = int(rng.integers(1000, sizes[contig] - size - 1000)) if start is None else start
            if not occupied[contig][max(0, s - 50): s + size + 50].any():
                contigs[contig][s: s + size] = codes
                occupied[contig][s: s + size] = True
                return s
            if start is not None:
                break
        raise RuntimeError("toy genome is too crowded")

    # large elements first, while there is still room for them
    monomer = _random_seq(rng, 171, 0.38)  # alpha-satellite-like tandem array
    for contig, n_mono in (("chrA", 110), ("chrB", 50)):
        array = np.concatenate([_mutate(rng, monomer, 0.02) for _ in range(n_mono)])
        s = place(contig, array)
        repeats.append((contig, s, s + array.size, "ALRToy", "Satellite/centr", 2.0))

    segdups = []  # (contig, start, end, other contig, other start, other end, fraction identical)
    for size, identity, src, dst in ((15_000, 0.995, "chrA", "chrB"), (8_000, 0.96, "chrA", "chrA"), (5_000, 0.92, "chrB", "chrA")):
        for _ in range(500):
            a = int(rng.integers(2000, sizes[src] - size - 2000))
            if not occupied[src][a: a + size].any():
                break
        occupied[src][a: a + size] = True
        copy = _mutate(rng, contigs[src][a: a + size], 1.0 - identity)
        b = place(dst, copy)
        segdups.append((src, a, a + size, dst, b, b + size, identity))
        segdups.append((dst, b, b + size, src, a, a + size, identity))

    families = [  # name, class/family, consensus length, copies, divergences to draw from
        ("AluToy", "SINE/Alu", 300, 110, (1.0, 8.0, 20.0)),
        ("L1Toy", "LINE/L1", 3000, 10, (2.0, 10.0)),
        ("ERVToy", "LTR/ERV1", 600, 18, (12.0,)),
        ("TiggerToy", "DNA/TcMar-Tigger", 400, 14, (18.0,)),
    ]
    for name, cls, size, copies, divergences in families:
        consensus = _random_seq(rng, size, 0.5)
        for _ in range(copies):
            contig = "chrA" if rng.random() < 0.63 else "chrB"
            div = float(rng.choice(divergences))
            copy = _mutate(rng, consensus, div / 100.0)
            if cls.startswith("LINE") and rng.random() < 0.5:  # 5'-truncated copies, as in real L1
                copy = copy[int(rng.integers(500, 2400)):]
            s = place(contig, copy)
            repeats.append((contig, s, s + copy.size, name, cls, div))

    for unit in ("CA", "AAAT", "GGC", "TTAGGG"):  # microsatellites
        for _ in range(5):
            contig = "chrA" if rng.random() < 0.63 else "chrB"
            n = int(rng.integers(120, 360))
            codes = np.array(["ACGT".index(c) for c in (unit * (n // len(unit) + 1))[:n]], np.uint8)
            s = place(contig, _mutate(rng, codes, 0.01))
            repeats.append((contig, s, s + n, f"({unit})n", "Simple_repeat", 1.0))
    for _ in range(12):  # A-rich low complexity
        contig = "chrA" if rng.random() < 0.63 else "chrB"
        codes = rng.choice(4, size=150, p=[0.82, 0.04, 0.04, 0.10]).astype(np.uint8)
        s = place(contig, codes)
        repeats.append((contig, s, s + 150, "A-rich", "Low_complexity", 5.0))

    # a diploid individual: ~1 SNV per kb, some indels, one hypervariable 3 kb region
    variants = []
    hyper = ("chrA", 60_000, 63_000)
    for contig, size in sizes.items():
        n_snv = int(size / 1000)
        positions = set(rng.choice(np.arange(200, size - 200), n_snv, replace=False).tolist())
        if contig == hyper[0]:
            positions |= set(rng.choice(np.arange(hyper[1], hyper[2]), 45, replace=False).tolist())
        for pos in sorted(positions):
            ref = int(contigs[contig][pos])
            alt = (ref + int(rng.integers(1, 4))) % 4
            gt = rng.choice(["0/1", "1/1", "0|1", "1|0"], p=[0.45, 0.3, 0.125, 0.125])
            variants.append((contig, pos, "ACGT"[ref], "ACGT"[alt], gt))
        for pos in rng.choice(np.arange(200, size - 200), int(size / 8000), replace=False):
            pos = int(pos)
            ref = "ACGT"[int(contigs[contig][pos])]
            if rng.random() < 0.5:  # deletion: REF = anchor + deleted bases, ALT = anchor
                k = int(rng.integers(1, 6))
                variants.append((contig, pos, ref + "".join("ACGT"[int(x)] for x in contigs[contig][pos + 1: pos + 1 + k]), ref, "0/1"))
            else:  # insertion: REF = anchor, ALT = anchor + inserted bases
                k = int(rng.integers(1, 6))
                variants.append((contig, pos, ref, ref + "".join("ACGT"[int(x)] for x in rng.integers(0, 4, k)), "0/1"))
    variants.sort()
    return {"contigs": contigs, "repeats": sorted(repeats), "segdups": sorted(segdups), "variants": variants,
            "sizes": sizes, "hypervariable": hyper}


def write_genome(world: dict, outdir: str) -> dict:
    """Write the toy world in the file formats of the real sources."""
    os.makedirs(outdir, exist_ok=True)
    paths = {k: os.path.join(outdir, v) for k, v in {
        "fasta": "toy.fa", "rmsk": "toy.fa.out.gz", "segdup": "toy.genomicSuperDups.txt.gz",
        "vcf": "toy_individual.vcf.gz", "confident": "toy_individual.confident.bed", "gff": "toy.gff.gz",
        "report": "toy_assembly_report.txt"}.items()}
    accession = {"chrA": ("TOY_000001.1", "TY000001.1", "A"), "chrB": ("TOY_000002.1", "TY000002.1", "B")}

    with open(paths["fasta"], "w") as fh:  # UCSC-style names
        for name, codes in world["contigs"].items():
            seq = BASES[codes].tobytes().decode()
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i: i + 60] + "\n")

    with gzip.open(paths["rmsk"], "wt") as fh:  # RepeatMasker .out: 1-based inclusive
        fh.write("   SW   perc perc perc  query     position in query    matching  repeat         position in repeat\n")
        fh.write("score   div. del. ins.  sequence  begin end   (left)   repeat    class/family   begin  end    (left)  ID\n\n")
        for i, (contig, s, e, name, cls, div) in enumerate(world["repeats"], start=1):
            left = world["sizes"][contig] - e
            fh.write(f" {2000:5d} {div:5.1f}  0.0  0.0  {contig}  {s + 1} {e}  ({left})  +  {name}  {cls}  1  {e - s}  (0)  {i}\n")

    with gzip.open(paths["segdup"], "wt") as fh:  # UCSC table: bin column first, 0-based half-open
        for i, (contig, s, e, other, os_, oe, identity) in enumerate(world["segdups"], start=1):
            row = [585, contig, s, e, f"{other}:{os_}", 0, "+", other, os_, oe, world["sizes"][other], i, 1000, "N/A",
                   "N/A", "N/A", "N/A", "align.txt", e - s, 0, 0, e - s, int((e - s) * identity),
                   int((e - s) * (1 - identity)), 0, 0, f"{identity:.6f}", f"{identity:.6f}", "0.01", "0.01"]
            fh.write("\t".join(str(x) for x in row) + "\n")

    with gzip.open(paths["vcf"], "wt") as fh:  # Ensembl-style names without "chr", 1-based POS
        fh.write("##fileformat=VCFv4.2\n")
        for name, size in world["sizes"].items():
            fh.write(f"##contig=<ID={accession[name][2]},length={size}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTOY_INDIVIDUAL\n")
        for contig, pos, ref, alt, gt in world["variants"]:
            fh.write(f"{accession[contig][2]}\t{pos + 1}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n")

    # "high-confidence" regions: everything except satellites and near-identical duplications,
    # mirroring how a real benchmark excludes what it cannot certify
    hard = {name: np.zeros(size, bool) for name, size in world["sizes"].items()}
    for contig, s, e, _, cls, _ in world["repeats"]:
        if cls.startswith("Satellite"):
            hard[contig][s:e] = True
    for contig, s, e, *_rest, identity in world["segdups"]:
        if identity >= 0.99:
            hard[contig][s:e] = True
    with open(paths["confident"], "w") as fh:
        for contig, mask in hard.items():
            edges = np.flatnonzero(np.diff(np.concatenate(([True], mask, [True])).astype(int)))
            for s, e in zip(edges[::2], edges[1::2]):
                fh.write(f"{contig}\t{s}\t{e}\n")

    with gzip.open(paths["gff"], "wt") as fh:  # RefSeq-style accessions, 1-based inclusive
        fh.write("##gff-version 3\n")
        gff_type = {"LINE": "mobile_genetic_element", "LTR": "long_terminal_repeat", "Satellite": "repeat_region",
                    "SINE": "mobile_genetic_element"}
        for contig, s, e, name, cls, _ in world["repeats"]:
            kind = gff_type.get(cls.split("/")[0])
            if kind:
                fh.write(f"{accession[contig][0]}\tToy\t{kind}\t{s + 1}\t{e}\t.\t+\t.\tID={name}_{s}\n")

    with open(paths["report"], "w") as fh:  # NCBI assembly report layout
        fh.write("# Assembly name:  ToyGenome1\n# Organism name:  Toy organism\n")
        fh.write("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\t"
                 "Relationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n")
        for name, (refseq, genbank, short) in accession.items():
            fh.write(f"{short}\tassembled-molecule\t{short}\tChromosome\t{genbank}\t=\t{refseq}\tPrimary Assembly\t"
                     f"{world['sizes'][name]}\t{name}\n")
    return paths


# --------------------------------------------------------------------------- #
# hidden instrument models that produce the stand-in "real" reads
# --------------------------------------------------------------------------- #
def illumina_like(seed: int = 1) -> dict:
    """A 2x150 short-read instrument: binned qualities that decay along the
    read, a worse second mate, and substitution-dominated errors."""
    rng = np.random.default_rng(seed)
    length = 150
    qvals = np.array([2, 12, 23, 37])
    q_first, q_trans, err, length_probs = [], [], [], []
    for mate in (0, 1):
        trans = np.zeros((MAX_CYCLES, QMAX, QMAX), np.float32)
        for c in range(1, length):
            decay = 0.004 + 0.035 * (c / length) ** 3 * (1.6 if mate else 1.0)
            for i, q in enumerate(qvals):
                row = np.zeros(4)
                row[i] = 1.0
                if i > 0:  # drop one quality bin
                    row[i] -= decay
                    row[i - 1] += decay
                if i < 3:  # recover one bin
                    row[i] -= 0.3
                    row[i + 1] += 0.3
                trans[c, q, qvals] = row * 1000
        first = np.zeros(QMAX)
        first[qvals] = [0.002, 0.01, 0.05, 0.938]
        q_first.append(first)
        q_trans.append(trans)
        table = np.zeros((CYCLE_BINS, QMAX))
        for b in range(CYCLE_BINS):
            # reported qualities are slightly optimistic, more so late in the read
            table[b] = np.clip(np.power(10.0, -np.arange(QMAX) / 10.0) * (1.3 + 0.4 * b), 0, 0.75)
        err.append(table)
        probs = np.zeros(MAX_CYCLES + 1)
        probs[length] = 1.0
        length_probs.append(probs)
    subst = np.array([[0, 4, 3, 1], [3, 0, 1, 5], [5, 1, 0, 3], [1, 3, 4, 0]], float)
    ins_len = np.zeros(INDEL_LEN_MAX + 1)
    ins_len[1:4] = [80, 15, 5]
    del_len = np.zeros(INDEL_LEN_MAX + 1)
    del_len[1:4] = [70, 20, 10]
    fragments = np.clip(rng.normal(360, 85, 20_000), 60, 900).astype(np.int64)
    meta = {"kind": "short", "paired": True, "source": "toy illumina-like instrument", "ins_rate": 2e-5,
            "del_rate": 5e-5, "mismatch_rate": 0.0, "mean_read_length": float(length)}
    return {"meta": meta, "q_first": np.stack(q_first), "q_trans": np.stack(q_trans),
            "length_probs": np.stack(length_probs), "err_by_q": np.stack(err), "subst": subst,
            "ins_len": ins_len, "del_len": del_len, "fragments": fragments,
            "cycle_err": np.zeros((2, MAX_CYCLES)), "cycle_meanq": np.zeros((2, MAX_CYCLES))}


def _long_profile(lengths, errors, fractions, hp_ins, hp_del, mean_indel, source) -> dict:
    ins_len = np.zeros(INDEL_LEN_MAX + 1)
    del_len = np.zeros(INDEL_LEN_MAX + 1)
    k = np.arange(1, INDEL_LEN_MAX + 1)
    geometric = (1 - 1 / mean_indel) ** (k - 1) / mean_indel
    ins_len[1:] = geometric * 1e6
    del_len[1:] = geometric * 1e6
    mean = float((k * geometric).sum() / geometric.sum())
    meta = {"kind": "long", "paired": False, "source": source, "frac_sub": fractions[0], "frac_ins": fractions[1],
            "frac_del": fractions[2], "mean_ins_len": mean, "mean_del_len": mean}
    return {"meta": meta, "lengths": lengths.astype(np.int64), "errors": errors, "ins_len": ins_len,
            "del_len": del_len, "subst": np.ones((4, 4)) - np.eye(4), "hp_ins_factor": hp_ins,
            "hp_del_factor": hp_del, "hp_den": np.ones(HPMAX + 1), "mean_quals": np.zeros(lengths.size)}


def ont_like(seed: int = 2) -> dict:
    """A nanopore-like instrument: broad read lengths, ~6 % errors that are
    mostly indels, deletions strongly enriched in homopolymers."""
    rng = np.random.default_rng(seed)
    lengths = np.clip(rng.lognormal(np.log(4500), 0.6, 4000), 300, 40_000)
    errors = np.clip(rng.beta(4, 60, 4000), 0.01, 0.3)
    hp = np.arange(HPMAX + 1, dtype=float)
    return _long_profile(lengths, errors, (0.35, 0.22, 0.43), 0.8 + 0.25 * hp, 0.6 + 0.7 * hp, 1.5,
                         "toy nanopore-like instrument")


def hifi_like(seed: int = 3) -> dict:
    """A HiFi-like instrument: tight read lengths, ~0.4 % errors that are
    almost all homopolymer indels."""
    rng = np.random.default_rng(seed)
    lengths = np.clip(rng.normal(8000, 1400, 4000), 2000, 16_000)
    errors = np.clip(rng.gamma(2.0, 0.002, 4000), 0.0003, 0.03)
    hp = np.arange(HPMAX + 1, dtype=float)
    return _long_profile(lengths, errors, (0.10, 0.48, 0.42), 0.3 + 1.6 * hp, 0.3 + 1.4 * hp, 1.2,
                         "toy HiFi-like instrument")
