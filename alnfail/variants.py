"""Truth-set variants for the simulated sample (tasks 2 and 7).

A real person's reads differ from the reference at their variants, and those
differences are exactly what makes some reads hard to place. Instead of
simulating reads from the bare reference we therefore simulate them from a
*sample haplotype*: the reference with a Genome in a Bottle truth set
applied. Each read is drawn from one of the two haplotypes.

Coordinates stay simple because we never build the haplotype genome. For
each read we take the reference window, apply only the variants inside it,
and carry along, for every base of the result, the reference position it
came from. The true origin of a read is then read straight off that array.

VCF convention: POS is 1-based and refers to the first base of REF. For
indels REF and ALT share their first (anchor) base.
"""
from __future__ import annotations

import gzip
import zlib
from dataclasses import dataclass

import numpy as np

from .coords import ContigResolver, vcf_to_bed
from .seqio import to_array


@dataclass
class Haplotype:
    """Variants carried by one haplotype of one contig, sorted by position and
    non-overlapping."""

    start: np.ndarray  # 0-based start of the REF allele
    ref_len: np.ndarray  # length of the REF allele
    alt: list  # ALT allele as a uint8 array

    def __len__(self) -> int:
        return int(self.start.size)


class VariantSet:
    """All variants of one sample: two haplotypes per contig, plus the plain
    list of variant positions for density counting."""

    def __init__(self):
        self.haps: dict[str, tuple[Haplotype, Haplotype]] = {}
        self.positions: dict[str, np.ndarray] = {}
        self.stats = {"records": 0, "used": 0, "skipped_overlap": 0, "skipped_other": 0, "off_reference": 0}

    @classmethod
    def from_vcf(cls, path: str, resolver: ContigResolver, sample_index: int = 0) -> "VariantSet":
        """Parse a (gzipped) VCF and split genotypes onto two haplotypes.

        Phased genotypes (0|1) are placed as written. Unphased heterozygous
        genotypes (0/1) are placed on a haplotype chosen by a hash of the
        position: arbitrary, but identical on every run, so the simulation
        stays reproducible.
        """
        self = cls()
        raw: dict[str, list[list]] = {}
        points: dict[str, list[int]] = {}
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                self.stats["records"] += 1
                contig = resolver.resolve(f[0])
                if contig is None:
                    self.stats["off_reference"] += 1
                    continue
                ref, alts = f[3], f[4].split(",")
                gt = f[9 + sample_index].split(":")[0]
                phased = "|" in gt
                alleles = gt.replace("|", "/").split("/")
                if any(a == "." for a in alleles) or any(not set(a) <= set("ACGTacgt") for a in alts) \
                        or not set(ref) <= set("ACGTNacgtn"):
                    self.stats["skipped_other"] += 1  # no-calls, symbolic alleles, breakends
                    continue
                alleles = [int(a) for a in alleles]
                if len(alleles) == 1:
                    alleles = alleles * 2
                if not any(alleles):
                    continue
                start0, _ = vcf_to_bed(int(f[1]), ref)
                if not phased and alleles[0] != alleles[1]:
                    if zlib.crc32(f"{contig}:{start0}".encode()) & 1:
                        alleles = alleles[::-1]
                store = raw.setdefault(contig, [[], []])
                for h in (0, 1):
                    if alleles[h] > 0:
                        store[h].append((start0, len(ref), alts[alleles[h] - 1].upper()))
                points.setdefault(contig, []).append(start0)
                self.stats["used"] += 1

        for contig, per_hap in raw.items():
            built = []
            for entries in per_hap:
                entries.sort()
                kept_start, kept_len, kept_alt = [], [], []
                last_end = -1
                for start0, ref_len, alt in entries:
                    if start0 < last_end:  # overlaps the previous variant on this haplotype
                        self.stats["skipped_overlap"] += 1
                        continue
                    kept_start.append(start0)
                    kept_len.append(ref_len)
                    kept_alt.append(to_array(alt))
                    last_end = start0 + ref_len
                built.append(Haplotype(np.array(kept_start, np.int64), np.array(kept_len, np.int64), kept_alt))
            self.haps[contig] = (built[0], built[1])
            self.positions[contig] = np.array(sorted(points[contig]), np.int64)
        return self

    def apply(self, contig: str, hap_index: int, window_start: int, ref: np.ndarray):
        """Apply this haplotype's variants to a reference window.

        ref           uint8 array holding reference bases [window_start, window_start + len)
        returns       (bases, ref_pos, n_variants)

        `ref_pos[i]` is the reference coordinate of output base i. Bases
        created by an insertion inherit the coordinate of the anchor base to
        their left, so min/max of `ref_pos` over any slice gives the slice's
        true reference interval.
        """
        n = ref.size
        positions = np.arange(window_start, window_start + n, dtype=np.int64)
        haps = self.haps.get(contig)
        if haps is None:
            return ref, positions, 0
        hap = haps[hap_index]
        lo = int(np.searchsorted(hap.start, window_start, side="left"))
        hi = int(np.searchsorted(hap.start, window_start + n, side="left"))
        if lo == hi:
            return ref, positions, 0

        seq_parts, pos_parts = [], []
        cursor = 0  # offset into the window
        applied = 0
        for i in range(lo, hi):
            v_off = int(hap.start[i]) - window_start
            v_len = int(hap.ref_len[i])
            if v_off < cursor or v_off + v_len > n:
                continue  # runs off the window end or overlaps: leave reference
            alt = hap.alt[i]
            seq_parts.append(ref[cursor:v_off])
            pos_parts.append(positions[cursor:v_off])
            seq_parts.append(alt)
            # ALT bases map onto the REF allele base by base; any surplus
            # (an insertion) keeps the coordinate of the last REF base
            alt_pos = window_start + v_off + np.minimum(np.arange(alt.size), v_len - 1)
            pos_parts.append(alt_pos.astype(np.int64))
            cursor = v_off + v_len
            applied += 1
        seq_parts.append(ref[cursor:])
        pos_parts.append(positions[cursor:])
        return np.concatenate(seq_parts), np.concatenate(pos_parts), applied


def mask_positions(variants: VariantSet | None, contig: str, pad: int = 2) -> np.ndarray:
    """Sorted reference positions to ignore when measuring sequencing error.

    A real sample's true variants look exactly like sequencing errors in an
    alignment. Counting them would inflate the error profile (human
    heterozygosity is about 1e-3 per base, the same order as Illumina's error
    rate), so known variant sites and a small pad around them are masked.
    """
    if variants is None or contig not in variants.haps:
        return np.zeros(0, np.int64)
    starts = np.concatenate([hap.start for hap in variants.haps[contig]])
    if starts.size == 0:
        return np.zeros(0, np.int64)
    widths = np.concatenate([hap.ref_len for hap in variants.haps[contig]]) + 2 * pad
    # expand every (start, width) into its run of positions without a Python loop
    first = np.repeat(starts - pad, widths)
    within = np.arange(int(widths.sum())) - np.repeat(np.cumsum(widths) - widths, widths)
    return np.unique(first + within)
