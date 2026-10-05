"""Coordinate conventions and contig naming (task 2).

Every file format in this project counts positions differently. Mixing them
up shifts every interval by one base, which is invisible in a plot and fatal
for an overlap-based benchmark. All conversions therefore live in this one
module and are unit-tested (tests/test_coords.py).

    format                         first base   end        example: bases 11..20
    -----------------------------  ----------   ---------  ---------------------
    BED, UCSC tables, BAM (pysam)  0            exclusive  start=10  end=20
    RepeatMasker .out              1            inclusive  begin=11  end=20
    GFF3 / GTF                     1            inclusive  start=11  end=20
    VCF                            1            (POS + len(REF) - 1)
    SAM text POS                   1            (computed from CIGAR)

Internally the whole package uses the BED convention: 0-based, half-open
intervals [start, end). The length of an interval is then simply end - start
and two intervals overlap exactly when start_a < end_b and start_b < end_a.
"""
from __future__ import annotations

from typing import Iterable


def closed1_to_bed(start1: int, end1: int) -> tuple[int, int]:
    """1-based inclusive (RepeatMasker .out, GFF) -> 0-based half-open (BED)."""
    if start1 < 1 or end1 < start1:
        raise ValueError(f"not a valid 1-based closed interval: {start1}-{end1}")
    return start1 - 1, end1


def bed_to_closed1(start0: int, end0: int) -> tuple[int, int]:
    """0-based half-open (BED) -> 1-based inclusive (GFF, region strings)."""
    if start0 < 0 or end0 <= start0:
        raise ValueError(f"not a valid 0-based half-open interval: {start0}-{end0}")
    return start0 + 1, end0


def vcf_to_bed(pos1: int, ref: str) -> tuple[int, int]:
    """VCF POS (1-based) plus REF allele -> the 0-based half-open reference
    interval the REF allele occupies."""
    if pos1 < 1 or not ref:
        raise ValueError(f"not a valid VCF record: POS={pos1} REF={ref!r}")
    return pos1 - 1, pos1 - 1 + len(ref)


def region_string(contig: str, start0: int, end0: int) -> str:
    """BED interval -> samtools/tabix region string, which is 1-based inclusive."""
    s1, e1 = bed_to_closed1(start0, end0)
    return f"{contig}:{s1}-{e1}"


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Number of shared bases between two half-open intervals (0 if disjoint)."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


# --------------------------------------------------------------------------- #
# contig names
# --------------------------------------------------------------------------- #
def _variants_of(name: str) -> list[str]:
    """Spellings of one contig name that differ only by archive convention:
    UCSC writes chr1 / chrM, Ensembl and GRC write 1 / MT."""
    out = [name]
    if name.startswith("chr"):
        bare = name[3:]
        out.append("MT" if bare == "M" else bare)
    else:
        out.append("chrM" if name == "MT" else "chr" + name)
    return out


class ContigResolver:
    """Translate contig names from any archive into the names used by our FASTA.

    The resolver knows the reference's own names, the trivial chr-prefix
    variants, and (optionally) an explicit alias table such as an NCBI
    assembly report or a UCSC chromAlias file, which is how GenBank
    (CM000663.2), RefSeq (NC_000001.11) and UCSC (chr1) names for the same
    sequence are tied together.

    Unknown names resolve to None and are counted, so that a silent mismatch
    (an annotation that matches no contig at all) cannot go unnoticed.
    """

    def __init__(self, reference_names: Iterable[str], aliases: dict[str, str] | None = None):
        self.reference = list(reference_names)
        self._map: dict[str, str] = {}
        for name in self.reference:
            for v in _variants_of(name):
                self._map.setdefault(v, name)
        for alias, target in (aliases or {}).items():
            canonical = self._map.get(target)
            if canonical is not None:
                self._map.setdefault(alias, canonical)
        # exact names always win over derived spellings
        for name in self.reference:
            self._map[name] = name
        self.unresolved: dict[str, int] = {}

    def resolve(self, name: str) -> str | None:
        hit = self._map.get(name)
        if hit is None:
            self.unresolved[name] = self.unresolved.get(name, 0) + 1
        return hit

    def report(self) -> str:
        if not self.unresolved:
            return "all contig names resolved"
        worst = sorted(self.unresolved.items(), key=lambda kv: -kv[1])[:8]
        shown = ", ".join(f"{k} ({v})" for k, v in worst)
        return f"{len(self.unresolved)} contig names not in the reference, e.g. {shown}"


class AliasGraph:
    """Groups of names that different archives use for the same sequence.

    Reads either of the two alias tables we meet:

      * NCBI ``*_assembly_report.txt``: one line per sequence with the columns
        Sequence-Name, GenBank-Accn, RefSeq-Accn and UCSC-style-name
        (columns 1, 5, 7 and 10). This is how ``CM000663.2`` (GenBank),
        ``NC_000001.11`` (RefSeq) and ``chr1`` (UCSC) are tied together.
      * UCSC ``*.chromAlias.txt``: every column is a spelling of the same name.
    """

    NCBI_COLUMNS = (0, 4, 6, 9)
    PLACEHOLDERS = {"na", "=", "<>", ""}

    def __init__(self, path: str):
        import gzip

        self.groups: list[set[str]] = []
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            lines = fh.read().splitlines()
        is_ncbi = any(line.startswith("# Sequence-Name") for line in lines)
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if is_ncbi:
                fields = [fields[i] for i in self.NCBI_COLUMNS if i < len(fields)]
            names = {f.strip() for f in fields} - self.PLACEHOLDERS
            if names:
                self.groups.append(names)

    def aliases_for(self, reference_names: Iterable[str]) -> dict[str, str]:
        """{alias: reference name} for every line that names exactly one
        reference sequence (directly or through a chr-prefix variant)."""
        ref = set(reference_names)
        out: dict[str, str] = {}
        for names in self.groups:
            hits = {n for n in names if n in ref}
            if not hits:
                hits = {v for n in names for v in _variants_of(n) if v in ref}
            if len(hits) == 1:
                target = next(iter(hits))
                for n in names:
                    out.setdefault(n, target)
        return out
