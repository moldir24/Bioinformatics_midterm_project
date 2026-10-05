"""Genome annotations -> strata (tasks 2 and 7).

Each annotation source (RepeatMasker, segmental duplications, a GFF, any BED)
is parsed into one normalised layout, BED4 with 0-based half-open coordinates
and contig names that match our FASTA:

    contig <tab> start <tab> end <tab> label

A set of such labelled intervals is a "stratum set". Reads are then assigned
to strata by how much of their true interval each label covers.

The parsers are the only place that knows the source formats, and each one
states which coordinate convention the source uses.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterator

import numpy as np

from .coords import ContigResolver, closed1_to_bed

Record = tuple[str, int, int, str]

# RepeatMasker reports dozens of class names; these are the ones large enough
# to stratify on. Everything else is pooled.
RM_MAJOR = ("SINE", "LINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Low_complexity", "Retroposon")
RM_RNA = {"rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "RNA"}
# classes whose divergence from the family consensus measures the age of the copy
RM_INTERSPERSED = {"SINE", "LINE", "LTR", "DNA", "Retroposon"}


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def normalise_rm_class(class_family: str) -> str:
    """'LINE/L1' -> 'LINE', 'DNA?' -> 'DNA', 'Satellite/centr' -> 'Satellite'."""
    cls = class_family.split("/")[0].rstrip("?")
    if cls in RM_MAJOR:
        return cls
    if cls in RM_RNA:
        return "RNA"
    return "Other_repeat"


def age_label(perc_div: float) -> str:
    """Bin a repeat copy by its divergence from the family consensus.

    Low divergence means a recent insertion: its copies are still nearly
    identical to one another, so a read from one copy matches many places.
    """
    if perc_div < 5.0:
        return "young_<5%"
    if perc_div < 15.0:
        return "middle_5-15%"
    return "old_>15%"


# --------------------------------------------------------------------------- #
# parsers: each yields (set_name, contig, start0, end0, label)
# --------------------------------------------------------------------------- #
def parse_rmsk_out(path: str) -> Iterator[tuple[str, str, int, int, str]]:
    """RepeatMasker .out (UCSC hg38.fa.out.gz, hs1.repeatMasker.out.gz).

    Whitespace-separated, three header lines, coordinates 1-based inclusive:
        SW  perc_div  perc_del  perc_ins  query  begin  end  (left)  strand  repeat  class/family ...
    Yields two stratum sets: `repeat_class` and, for interspersed repeats,
    `repeat_age` from the per-copy divergence column.
    """
    with _open(path) as fh:
        for line in fh:
            f = line.split()
            # header and blank lines do not start with an integer score
            if len(f) < 11 or not f[0].isdigit():
                continue
            start0, end0 = closed1_to_bed(int(f[5]), int(f[6]))
            cls = normalise_rm_class(f[10])
            yield "repeat_class", f[4], start0, end0, cls
            if cls in RM_INTERSPERSED:
                yield "repeat_age", f[4], start0, end0, age_label(float(f[1]))


def parse_ucsc_rmsk_table(path: str) -> Iterator[tuple[str, str, int, int, str]]:
    """UCSC database table rmsk.txt.gz. Tab-separated, 0-based half-open:
        bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand repName repClass repFamily ...
    milliDiv is divergence in parts per thousand."""
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            cls = normalise_rm_class(f[11])
            start0, end0 = int(f[6]), int(f[7])
            yield "repeat_class", f[5], start0, end0, cls
            if cls in RM_INTERSPERSED:
                yield "repeat_age", f[5], start0, end0, age_label(int(f[2]) / 10.0)


def sd_label(frac_match: float | None) -> str:
    """Bin a segmental duplication by the identity between its two copies."""
    if frac_match is None:
        return "SD"
    if frac_match >= 0.99:
        return "SD_>=99%"
    if frac_match >= 0.95:
        return "SD_95-99%"
    return "SD_<95%"


def parse_ucsc_superdups(path: str) -> Iterator[tuple[str, str, int, int, str]]:
    """UCSC genomicSuperDups.txt.gz. Tab-separated, 0-based half-open, with a
    leading `bin` column:
        bin chrom chromStart chromEnd name score strand otherChrom ... fracMatch(col 27) ...
    """
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            frac = None
            if len(f) > 26:
                try:
                    frac = float(f[26])
                except ValueError:
                    frac = None
            yield "segdup", f[1], int(f[2]), int(f[3]), sd_label(frac)


def parse_bed(path: str, set_name: str, label: str | None = None, label_col: int | None = None,
              identity_col: int | None = None) -> Iterator[tuple[str, str, int, int, str]]:
    """Any BED file (already 0-based half-open).

    label        fixed label for every interval, or
    label_col    0-based column holding the label, or
    identity_col 0-based column holding a copy identity (fraction or percent),
                 binned like a segmental duplication.
    """
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].isdigit():
                continue
            if identity_col is not None and len(f) > identity_col:
                try:
                    ident = float(f[identity_col])
                    ident = ident / 100.0 if ident > 1.0 else ident
                    lab = sd_label(ident)
                except ValueError:
                    lab = label or set_name
            elif label_col is not None and len(f) > label_col:
                lab = f[label_col]
            else:
                lab = label or set_name
            yield set_name, f[0], int(f[1]), int(f[2]), lab


def parse_gff(path: str, types: list[str], set_name: str = "feature") -> Iterator[tuple[str, str, int, int, str]]:
    """GFF3 features of the requested types (1-based inclusive coordinates).

    For genomes without RepeatMasker tracks (bacteria, yeast) the repeated
    sequence classes are annotated as features: insertion sequences and
    transposons (mobile_genetic_element), rRNA operons, tRNAs, LTRs.
    """
    wanted = set(types)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in wanted:
                continue
            start0, end0 = closed1_to_bed(int(f[3]), int(f[4]))
            yield set_name, f[0], start0, end0, f[2]


PARSERS = {
    "rmsk_out": lambda path, name, opts: parse_rmsk_out(path),
    "ucsc_rmsk": lambda path, name, opts: parse_ucsc_rmsk_table(path),
    "ucsc_superdups": lambda path, name, opts: parse_ucsc_superdups(path),
    "bed": lambda path, name, opts: parse_bed(
        path, opts.get("set", name), opts.get("label"), opts.get("label_col"), opts.get("identity_col")
    ),
    "gff": lambda path, name, opts: parse_gff(path, opts.get("types", []), opts.get("set", "feature")),
}


def build_strata(path: str, fmt: str, name: str, opts: dict, resolver: ContigResolver,
                 contig_lengths: dict[str, int]) -> tuple[dict[str, list[Record]], dict]:
    """Parse one annotation source into normalised records, checking every
    interval against the reference.

    Returns ({set_name: [records]}, stats). Intervals on contigs we do not use
    are dropped and counted; intervals that run past the end of their contig
    raise an error, because that means the annotation belongs to a different
    assembly version than the FASTA.
    """
    if fmt not in PARSERS:
        raise ValueError(f"unknown annotation format {fmt!r}; known: {sorted(PARSERS)}")
    sets: dict[str, list[Record]] = defaultdict(list)
    stats = {"source": name, "format": fmt, "records_in": 0, "records_kept": 0, "off_reference": 0}
    for set_name, contig, start0, end0, label in PARSERS[fmt](path, name, opts):
        stats["records_in"] += 1
        canonical = resolver.resolve(contig)
        if canonical is None:
            stats["off_reference"] += 1
            continue
        if start0 < 0 or end0 <= start0:
            raise ValueError(f"{name}: empty or negative interval {contig}:{start0}-{end0}")
        if end0 > contig_lengths[canonical]:
            raise ValueError(
                f"{name}: interval {contig}:{start0}-{end0} ends after the contig "
                f"({contig_lengths[canonical]} bp). Annotation and FASTA are different assembly versions."
            )
        sets[set_name].append((canonical, start0, end0, label))
        stats["records_kept"] += 1
    stats["contig_names"] = resolver.report()
    if stats["records_in"] and not stats["records_kept"]:
        raise ValueError(f"{name}: no interval matched any reference contig ({resolver.report()})")
    return sets, stats


def write_bed(records: list[Record], path: str) -> None:
    records = sorted(records, key=lambda r: (r[0], r[1], r[2]))
    with open(path, "w") as out:
        for contig, start0, end0, label in records:
            out.write(f"{contig}\t{start0}\t{end0}\t{label}\n")


def read_bed(path: str) -> list[Record]:
    out: list[Record] = []
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                out.append((f[0], int(f[1]), int(f[2]), f[3]))
    return out


# --------------------------------------------------------------------------- #
# interval lookup
# --------------------------------------------------------------------------- #
@dataclass
class _Track:
    """Merged, sorted, non-overlapping intervals of one label on one contig."""

    starts: np.ndarray
    ends: np.ndarray
    cum: np.ndarray = field(init=False)  # covered bases before each interval

    def __post_init__(self):
        lengths = self.ends - self.starts
        self.cum = np.concatenate(([0], np.cumsum(lengths)))

    def covered_before(self, x: np.ndarray) -> np.ndarray:
        """Number of covered bases in [0, x) for each x."""
        idx = np.searchsorted(self.starts, x, side="right") - 1
        safe = np.clip(idx, 0, None)
        inside = np.clip(x - self.starts[safe], 0, self.ends[safe] - self.starts[safe])
        return np.where(idx >= 0, self.cum[safe] + inside, 0)

    def coverage(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """Covered bases inside each query interval [start, end)."""
        return self.covered_before(ends) - self.covered_before(starts)


def _merge(intervals: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    arr = np.array(sorted(intervals), dtype=np.int64).reshape(-1, 2)
    if arr.shape[0] == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    starts, ends = arr[:, 0], arr[:, 1]
    # an interval starts a new merged block when it begins after everything before it ended
    running_end = np.maximum.accumulate(ends)
    new_block = np.concatenate(([True], starts[1:] > running_end[:-1]))
    block_id = np.cumsum(new_block) - 1
    m_starts = starts[new_block]
    m_ends = np.zeros(m_starts.size, dtype=np.int64)
    np.maximum.at(m_ends, block_id, ends)
    return m_starts, m_ends


class StratumSet:
    """One stratum set: labelled intervals with fast coverage queries."""

    def __init__(self, name: str, records: list[Record]):
        self.name = name
        grouped: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
        for contig, start0, end0, label in records:
            grouped[label][contig].append((start0, end0))
        self.labels = sorted(grouped)
        self.tracks: dict[str, dict[str, _Track]] = {
            label: {contig: _Track(*_merge(iv)) for contig, iv in by_contig.items()}
            for label, by_contig in grouped.items()
        }

    def coverage(self, contigs: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> dict[str, np.ndarray]:
        """Covered bases per label for each query interval."""
        out = {label: np.zeros(starts.size, dtype=np.int64) for label in self.labels}
        for contig in np.unique(contigs):
            sel = np.flatnonzero(contigs == contig)
            for label in self.labels:
                track = self.tracks[label].get(contig)
                if track is not None and track.starts.size:
                    out[label][sel] = track.coverage(starts[sel], ends[sel])
        return out

    def assign(self, contigs, starts, ends, min_frac: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
        """Assign each interval to its dominant label.

        Returns (label, fraction). A read gets a label only when that label
        covers at least `min_frac` of it; otherwise it is 'none'. `fraction`
        is the share of the read covered by the best label, which lets the
        analysis separate "read touches a repeat" from "read lies entirely
        inside one".
        """
        contigs = np.asarray(contigs, dtype=object)
        starts = np.asarray(starts, dtype=np.int64)
        ends = np.asarray(ends, dtype=np.int64)
        length = np.maximum(ends - starts, 1)
        best_label = np.full(starts.size, "none", dtype=object)
        best_frac = np.zeros(starts.size, dtype=float)
        for label, cov in self.coverage(contigs, starts, ends).items():
            frac = cov / length
            better = frac > best_frac
            best_frac[better] = frac[better]
            best_label[better] = label
        best_label[best_frac < min_frac] = "none"
        return best_label, best_frac


def load_strata(paths: dict[str, str]) -> dict[str, StratumSet]:
    """{set_name: path to normalised BED} -> {set_name: StratumSet}."""
    return {name: StratumSet(name, read_bed(path)) for name, path in paths.items()}


# --------------------------------------------------------------------------- #
# numeric strata
# --------------------------------------------------------------------------- #
GC_EDGES = (0.30, 0.40, 0.50, 0.60, 0.70)
GC_LABELS = ("<30%", "30-40%", "40-50%", "50-60%", "60-70%", ">=70%")


def gc_bin(gc: np.ndarray) -> np.ndarray:
    """GC fraction -> label. Bin edges sit at round numbers so the extreme
    bins isolate the AT-rich and GC-rich sequence where coverage and accuracy
    are known to drop."""
    idx = np.searchsorted(np.array(GC_EDGES), np.nan_to_num(np.asarray(gc, float), nan=0.45), side="right")
    return np.array(GC_LABELS, dtype=object)[idx]


def variant_density_bin(n_variants: np.ndarray, length: np.ndarray) -> np.ndarray:
    """Variants per kilobase of read -> label.

    Counts are normalised by read length so that a 150 bp read and a 15 kb
    read can share one scale. 1 per kb is roughly the genome-wide human
    average; above 5 per kb is the range of hypervariable loci such as HLA.
    """
    per_kb = np.asarray(n_variants, float) / np.maximum(np.asarray(length, float), 1.0) * 1000.0
    labels = np.array(["0", "0-2/kb", "2-5/kb", "5-10/kb", ">=10/kb"], dtype=object)
    idx = np.where(per_kb <= 0, 0, 1 + np.searchsorted(np.array([2.0, 5.0, 10.0]), per_kb, side="right"))
    return labels[idx]


def length_bin(length: np.ndarray) -> np.ndarray:
    """Read length -> label (meaningful for long reads)."""
    edges = np.array([500, 1000, 2000, 5000, 10000, 20000, 50000])
    labels = np.array(["<500", "0.5-1k", "1-2k", "2-5k", "5-10k", "10-20k", "20-50k", ">=50k"], dtype=object)
    return labels[np.searchsorted(edges, np.asarray(length), side="right")]


def count_in_intervals(positions_by_contig: dict[str, np.ndarray], contigs, starts, ends) -> np.ndarray:
    """How many sorted point positions fall inside each interval [start, end)."""
    contigs = np.asarray(contigs, dtype=object)
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    out = np.zeros(starts.size, dtype=np.int64)
    for contig in np.unique(contigs):
        pos = positions_by_contig.get(contig)
        if pos is None or pos.size == 0:
            continue
        sel = np.flatnonzero(contigs == contig)
        out[sel] = np.searchsorted(pos, ends[sel], side="left") - np.searchsorted(pos, starts[sel], side="left")
    return out
