"""Coordinate conventions: the place where off-by-one errors are caught."""
import gzip

import numpy as np
import pytest

from alnfail.annot import StratumSet, build_strata, parse_gff, parse_rmsk_out, parse_ucsc_superdups
from alnfail.coords import AliasGraph, ContigResolver, bed_to_closed1, closed1_to_bed, overlap, region_string, vcf_to_bed


def test_closed1_and_bed_describe_the_same_bases():
    # bases 11..20 (1-based, inclusive) are bases 10..19 (0-based), i.e. [10, 20)
    assert closed1_to_bed(11, 20) == (10, 20)
    assert bed_to_closed1(10, 20) == (11, 20)
    for start1, end1 in [(1, 1), (1, 100), (57, 58)]:
        start0, end0 = closed1_to_bed(start1, end1)
        assert end0 - start0 == end1 - start1 + 1  # length is preserved
        assert bed_to_closed1(start0, end0) == (start1, end1)


def test_invalid_intervals_are_rejected():
    with pytest.raises(ValueError):
        closed1_to_bed(0, 5)  # 1-based coordinates cannot be 0
    with pytest.raises(ValueError):
        bed_to_closed1(5, 5)  # empty half-open interval


def test_vcf_position_covers_the_ref_allele():
    assert vcf_to_bed(100, "A") == (99, 100)  # SNV
    assert vcf_to_bed(100, "ACGT") == (99, 103)  # deletion: anchor + 3 deleted bases


def test_region_string_is_one_based_inclusive():
    assert region_string("chr20", 0, 1000) == "chr20:1-1000"


def test_overlap_of_half_open_intervals():
    assert overlap(0, 10, 10, 20) == 0  # touching intervals share no base
    assert overlap(0, 10, 9, 20) == 1
    assert overlap(5, 8, 0, 100) == 3


def test_resolver_handles_chr_prefix_and_mitochondrion():
    r = ContigResolver(["chr1", "chrM"])
    assert r.resolve("1") == "chr1" and r.resolve("chr1") == "chr1"
    assert r.resolve("MT") == "chrM"
    assert r.resolve("chrUn_1") is None and "chrUn_1" in r.unresolved
    r = ContigResolver(["1", "MT"])  # the other way round
    assert r.resolve("chr1") == "1" and r.resolve("chrM") == "MT"


def test_alias_graph_ties_archive_names_together(tmp_path):
    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\t"
        "Relationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1\n"
        "2\tassembled-molecule\t2\tChromosome\tCM000664.2\t=\tNC_000002.12\tPrimary Assembly\t242193529\tchr2\n")
    aliases = AliasGraph(str(report)).aliases_for(["chr1", "chr2"])
    r = ContigResolver(["chr1", "chr2"], aliases)
    assert r.resolve("NC_000001.11") == "chr1"  # RefSeq
    assert r.resolve("CM000664.2") == "chr2"  # GenBank
    assert r.resolve("Chromosome") is None  # a column value, not a sequence name


def test_repeatmasker_out_is_converted_from_one_based(tmp_path):
    path = tmp_path / "rm.out.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("   SW   perc perc perc  query  position in query  matching  repeat  position in repeat\n")
        fh.write("score   div. del. ins.  sequence  begin end (left) repeat class/family begin end (left) ID\n\n")
        fh.write("  463   1.3  0.6  1.7  chr1  10001  10468  (248945954)  +  (TAACCC)n  Simple_repeat  1  471  (0)  1\n")
        fh.write(" 3612  11.4 21.5  1.3  chr1  11505  11675  (248944747)  C  L1MC5a  LINE/L1  (2382)  5648  5452  2\n")
    records = list(parse_rmsk_out(str(path)))
    assert ("repeat_class", "chr1", 10000, 10468, "Simple_repeat") in records
    assert ("repeat_class", "chr1", 11504, 11675, "LINE") in records
    assert ("repeat_age", "chr1", 11504, 11675, "middle_5-15%") in records  # 11.4 % diverged from consensus
    assert not any(r[0] == "repeat_age" and r[4].startswith("young") for r in records)  # simple repeats have no age


def test_superdups_table_is_already_zero_based(tmp_path):
    path = tmp_path / "sd.txt"
    fields = ["585", "chr1", "10000", "20000", "chr2:5", "0", "+", "chr2", "5", "10005", "1000"] + ["x"] * 15 + ["0.995"]
    path.write_text("\t".join(fields) + "\n")
    assert list(parse_ucsc_superdups(str(path))) == [("segdup", "chr1", 10000, 20000, "SD_>=99%")]


def test_gff_is_converted_from_one_based(tmp_path):
    path = tmp_path / "a.gff"
    path.write_text("##gff-version 3\nNC_1\tRefSeq\tmobile_genetic_element\t101\t200\t.\t+\t.\tID=x\n"
                    "NC_1\tRefSeq\tgene\t1\t50\t.\t+\t.\tID=y\n")
    assert list(parse_gff(str(path), ["mobile_genetic_element"])) == [("feature", "NC_1", 100, 200, "mobile_genetic_element")]


def test_annotation_from_another_assembly_is_refused(tmp_path):
    path = tmp_path / "b.bed"
    path.write_text("chr1\t900\t1100\n")  # runs past the end of a 1,000 bp contig
    with pytest.raises(ValueError, match="different assembly"):
        build_strata(str(path), "bed", "x", {}, ContigResolver(["chr1"]), {"chr1": 1000})


def test_stratum_assignment_counts_covered_bases():
    sset = StratumSet("rep", [("chr1", 100, 200, "SINE"), ("chr1", 150, 260, "SINE"), ("chr1", 500, 600, "LINE")])
    contigs = np.array(["chr1"] * 4, dtype=object)
    starts = np.array([0, 100, 250, 550])
    ends = np.array([100, 200, 350, 650])
    label, frac = sset.assign(contigs, starts, ends)
    # [0,100) touches nothing; [100,200) lies inside the merged SINE block [100,260);
    # [250,350) is covered for 10 bases only; [550,650) is half inside the LINE
    assert list(label) == ["none", "SINE", "none", "LINE"]
    assert np.allclose(frac, [0.0, 1.0, 0.1, 0.5])
