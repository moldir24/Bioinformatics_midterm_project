"""Build the small test dataset shipped in test_data/ (see toy.py).

Run once with

    python -m alnfail make-test-data --outdir test_data

and commit the result. Needs minimap2 and samtools on the PATH, because the
stand-in reads are delivered the way the real ones are: as sorted, indexed
BAM files that the pipeline slices by region.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile

from . import toy
from .coords import ContigResolver
from .simulate import simulate_long, simulate_short
from .variants import VariantSet

COVERAGE = {"illumina": 18, "ont": 14, "hifi": 12}


def _align_sorted(fasta: str, preset: str, fastqs: list[str], out_bam: str) -> None:
    cmd = (f"minimap2 -a -x {preset} -t 2 '{fasta}' {' '.join(repr(f) for f in fastqs)} 2>/dev/null"
           f" | samtools sort -@ 2 -o '{out_bam}' - && samtools index '{out_bam}'")
    subprocess.run(["bash", "-c", "set -euo pipefail; " + cmd], check=True)


def make_test_data(outdir: str, seed: int = 26) -> None:
    world = toy.build_genome(seed)
    paths = toy.write_genome(world, outdir)
    import pysam

    pysam.faidx(paths["fasta"])
    genome_size = sum(world["sizes"].values())
    resolver = ContigResolver(world["sizes"])
    individual = VariantSet.from_vcf(paths["vcf"], resolver)
    reads_dir = os.path.join(outdir, "reads")
    os.makedirs(reads_dir, exist_ok=True)

    summary = {"genome_bp": genome_size, "repeats": len(world["repeats"]), "segdups": len(world["segdups"]),
               "variants": len(world["variants"]), "note": "synthetic test world, see alnfail/toy.py"}
    with tempfile.TemporaryDirectory() as scratch:
        prefix = os.path.join(scratch, "illumina")
        n = COVERAGE["illumina"] * genome_size // 300
        simulate_short(paths["fasta"], toy.illumina_like(), n, prefix + "_1.fastq.gz", prefix + "_2.fastq.gz",
                       prefix + ".parquet", seed, "TOYILMN", variants=individual,
                       adapters=(toy.ADAPTER_R1, toy.ADAPTER_R2))
        _align_sorted(paths["fasta"], "sr", [prefix + "_1.fastq.gz", prefix + "_2.fastq.gz"],
                      os.path.join(reads_dir, "toy_illumina.bam"))
        summary["illumina_pairs"] = int(n)

        for name, model, preset, mean_len in (("ont", toy.ont_like(), "map-ont", 5400), ("hifi", toy.hifi_like(), "map-hifi", 8000)):
            prefix = os.path.join(scratch, name)
            n = COVERAGE[name] * genome_size // mean_len
            simulate_long(paths["fasta"], model, n, prefix + ".fastq.gz", prefix + ".parquet", seed,
                          f"TOY{name.upper()}", variants=individual)
            _align_sorted(paths["fasta"], preset, [prefix + ".fastq.gz"], os.path.join(reads_dir, f"toy_{name}.bam"))
            summary[f"{name}_reads"] = int(n)

    os.remove(paths["fasta"] + ".fai")  # the pipeline builds its own index
    with open(os.path.join(outdir, "test_data.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
