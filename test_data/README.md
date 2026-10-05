# Test dataset

A miniature world for checking the pipeline end to end in a few minutes,
without a network connection:

```
snakemake --cores 4 --configfile config/test.yaml
```

**Everything in this directory is synthetic.** It is a stand-in that has the
same file formats and the same kinds of features as the real inputs, so that
every rule of the workflow is exercised. No number produced from it is a
result, and none belongs in the report.

| File | Stands in for | Format (as in the real source) |
|---|---|---|
| `toy.fa` | a reference genome: two contigs, 350 kb, with GC-poor and GC-rich blocks | FASTA, UCSC-style names (`chrA`, `chrB`) |
| `toy.fa.out.gz` | RepeatMasker annotation: SINE, LINE, LTR, DNA, satellite, simple repeats, low complexity, with per-copy divergence | RepeatMasker `.out`, 1-based inclusive |
| `toy.genomicSuperDups.txt.gz` | segmental duplications at 99.5 %, 96 % and 92 % identity | UCSC table, 0-based, leading `bin` column |
| `toy_individual.vcf.gz` | a truth set for one diploid individual: SNVs, indels, one hypervariable region | VCF, 1-based, contig names without `chr` (`A`, `B`) |
| `toy_individual.confident.bed` | high-confidence regions (everything except satellites and the 99.5 % duplication) | BED |
| `toy.gff.gz` | feature annotation, as for bacteria and yeast | GFF3, 1-based, RefSeq-style accessions (`TOY_000001.1`) |
| `toy_assembly_report.txt` | the table tying the three naming schemes together | NCBI assembly report |
| `reads/toy_{illumina,ont,hifi}.bam` (+ `.bai`) | a data provider's aligned reads of that individual on three platforms | sorted, indexed BAM |

The contig names deliberately differ between files (`chrA`, `A`,
`TOY_000001.1`), so the test also checks identifier mapping.

The reads were produced by three hidden "instrument" models defined in
`alnfail/toy.py`: a 2×150 short-read instrument with binned qualities,
adapter read-through and a worse second mate; a nanopore-like instrument with
about 6 % indel-rich errors; a HiFi-like instrument with about 0.4 % errors
concentrated in homopolymers. The pipeline does not know those models. It has
to measure them from the BAM files, which is what makes this a real test of
the calibration step.

To rebuild the directory (needs minimap2 and samtools):

```
python -m alnfail make-test-data --outdir test_data
```
