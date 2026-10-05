"""The simulator must know exactly where every read came from."""
import numpy as np
import pysam
import pytest

from alnfail import toy
from alnfail.coords import ContigResolver
from alnfail.seqio import BASES, CODE, homopolymer_run_lengths, revcomp, to_array
from alnfail.simulate import _apply_indels, apply_divergence, simulate_long, simulate_short
from alnfail.variants import VariantSet


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """A 30 kb random genome and a small VCF."""
    root = tmp_path_factory.mktemp("genome")
    rng = np.random.default_rng(1)
    codes = rng.integers(0, 4, 30_000).astype(np.uint8)
    seq = BASES[codes].tobytes().decode()
    fasta = root / "g.fa"
    fasta.write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    pysam.faidx(str(fasta))
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    for pos in range(500, 29_000, 700):  # homozygous SNVs, insertions and deletions in turn
        ref = seq[pos]
        kind = (pos // 700) % 3
        if kind == 0:
            lines.append(f"1\t{pos + 1}\t.\t{ref}\t{'ACGT'[('ACGT'.index(ref) + 1) % 4]}\t.\t.\t.\tGT\t1/1")
        elif kind == 1:
            lines.append(f"1\t{pos + 1}\t.\t{ref}\t{ref}TTT\t.\t.\t.\tGT\t1/1")
        else:
            lines.append(f"1\t{pos + 1}\t.\t{seq[pos:pos + 4]}\t{ref}\t.\t.\t.\tGT\t1/1")
    vcf = root / "v.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    return {"fasta": str(fasta), "seq": seq, "vcf": str(vcf), "root": root}


def test_reverse_complement_and_homopolymers():
    assert revcomp(to_array("AACGT")).tobytes() == b"ACGTT"
    assert list(homopolymer_run_lengths(to_array("AAACGG"))) == [3, 3, 3, 1, 2, 2]


def test_indels_keep_reference_positions_in_step():
    codes = np.array([0, 1, 2, 3, 0, 1], np.uint8)
    rpos = np.arange(100, 106)
    out, pos = _apply_indels(codes, rpos, np.array([1]), np.array([2]), np.array([4]), np.array([2]),
                             np.array([3, 3], np.uint8))
    assert list(out) == [0, 3, 0, 3, 3, 1]  # bases 101-102 deleted, two bases inserted after 104
    assert list(pos) == [100, 103, 104, 104, 104, 105]  # inserted bases carry their anchor's position


def test_divergence_rate_is_what_was_asked_for():
    rng = np.random.default_rng(2)
    codes = rng.integers(0, 4, 200_000).astype(np.uint8)
    _, _, events = apply_divergence(codes, np.arange(codes.size), 0.05, rng)
    assert 0.045 < events / codes.size < 0.055


def test_variants_change_bases_but_not_coordinates(genome):
    variants = VariantSet.from_vcf(genome["vcf"], ContigResolver(["chr1"]))  # VCF says "1", FASTA says "chr1"
    assert variants.stats["used"] == len(variants.positions["chr1"]) > 30
    ref = to_array(genome["seq"][400:1400])
    seq, pos, n = variants.apply("chr1", 0, 400, ref)
    assert n >= 1 and pos[0] == 400 and pos[-1] == 1399
    assert np.all(np.diff(pos) >= 0)  # reference positions never go backwards
    assert seq.size != ref.size or not np.array_equal(seq, ref)


def test_short_reads_come_from_where_the_truth_table_says(genome):
    root = genome["root"]
    prof = toy.illumina_like()
    info = simulate_short(genome["fasta"], prof, 400, str(root / "s_1.fq.gz"), str(root / "s_2.fq.gz"),
                          str(root / "s.parquet"), 7, "t")
    assert info["fragments"] == 400 and info["paired"]
    import pandas as pd

    from alnfail.seqio import read_fastq

    truth = pd.read_parquet(root / "s.parquet").set_index(["qname", "mate"])
    assert len(truth) == 800
    for mate, path in ((1, root / "s_1.fq.gz"), (2, root / "s_2.fq.gz")):
        for name, read, qual in read_fastq(str(path)):
            row = truth.loc[(name, mate)]
            assert len(read) == len(qual) == row["length"]
            origin = to_array(genome["seq"][row["start"]:row["end"]])
            if row["strand"] == "-":
                origin = revcomp(origin)
            if origin.size == len(read):  # no indel error in this read: compare base by base
                mismatches = int((CODE[origin] != CODE[to_array(read)]).sum())
                assert mismatches <= max(8, 0.1 * len(read)), (name, mate, mismatches)
    # mates of a pair face each other on opposite strands
    pairs = truth.reset_index().pivot(index="qname", columns="mate", values="strand")
    assert (pairs[1] != pairs[2]).all()


def test_long_reads_come_from_where_the_truth_table_says(genome):
    root = genome["root"]
    prof = toy.hifi_like()
    prof["lengths"] = np.clip(prof["lengths"], 1000, 4000)
    simulate_long(genome["fasta"], prof, 60, str(root / "l.fq.gz"), str(root / "l.parquet"), 7, "t")
    import pandas as pd

    from alnfail.seqio import read_fastq

    truth = pd.read_parquet(root / "l.parquet").set_index("qname")
    for name, read, _ in read_fastq(str(root / "l.fq.gz")):
        row = truth.loc[name]
        assert 0 <= row["start"] < row["end"] <= 30_000
        # a HiFi-like read spans about as many reference bases as it has bases
        assert abs((row["end"] - row["start"]) - len(read)) < 0.05 * len(read) + 20
        origin = genome["seq"][row["start"]:row["end"]]
        if row["strand"] == "-":
            origin = revcomp(to_array(origin)).tobytes().decode()
        # a large share of the 15-mers of the read occur in its recorded origin (almost none would in a random window)
        kmers = {origin[i:i + 15] for i in range(len(origin) - 14)}
        shared = sum(read[i:i + 15] in kmers for i in range(0, len(read) - 14, 5))
        assert shared > 0.4 * len(range(0, len(read) - 14, 5)), name


def test_same_seed_gives_identical_reads(genome):
    root = genome["root"]
    prof = toy.ont_like()
    prof["lengths"] = np.clip(prof["lengths"], 500, 3000)
    for tag in ("a", "b"):
        simulate_long(genome["fasta"], prof, 20, str(root / f"{tag}.fq.gz"), str(root / f"{tag}.parquet"), 11, "same")
    import gzip

    assert gzip.open(root / "a.fq.gz").read() == gzip.open(root / "b.fq.gz").read()
