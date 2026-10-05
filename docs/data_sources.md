# Data sources, formats and coordinate conventions (task 2)

Everything is retrieved by the workflow (`workflow/rules/data.smk`); nothing is
downloaded by hand. After a run, `results/tables/provenance.tsv` lists every
source with URL, size, MD5 and retrieval time, and
`results/tables/strata_sources.tsv` shows how many intervals of each annotation
matched our contigs. Copy the provenance table into the report's data section.

**Before anything else, run the source check** (handbook, week 1: "confirm your
accessions resolve"):

```
snakemake check_sources --cores 1 --configfile config/config.yaml
```

It contacts every URL and every BAM header in the config and prints size and
status. The URLs in `config/config.yaml` were taken from the providers'
documentation but were *not* test-downloaded when this repository was written;
this command is the test. If one has moved, fix it in the config.

## Four archives

| Archive | What we take | Why this archive |
|---|---|---|
| UCSC Genome Browser downloads | GRCh38 (`hg38`) and T2T-CHM13v2.0 (`hs1`) sequence, RepeatMasker output, segmental duplications | sequence and annotation with matching `chr` names, one RepeatMasker format for both assemblies |
| NCBI RefSeq | *E. coli* K-12 MG1655 (GCF_000005845.2) and *S. cerevisiae* S288C (GCF_000146045.2): sequence, GFF3, assembly report | curated reference assemblies with feature annotation |
| Genome in a Bottle (NIST), hosted on the NCBI FTP site | HG002 v4.2.1 truth variants and high-confidence regions; HG002 Illumina, PacBio HiFi and ONT reads as aligned BAMs | the best-characterised human sample; paths and checksums are listed in the `giab_data_indexes` repository |
| T2T consortium (AWS) | CHM13 segmental duplications | produced by the assembly's authors |

## Identifier mapping

The same chromosome has a different name in each archive:

| Archive style | Chromosome 20 of GRCh38 | *E. coli* chromosome |
|---|---|---|
| UCSC | `chr20` | — |
| Ensembl / GRC / many VCFs | `20` | `Chromosome` |
| GenBank | `CM000682.2` | `U00096.3` |
| RefSeq | `NC_000020.11` | `NC_000913.3` |

`alnfail/coords.py` resolves names in two steps: trivial spelling variants
(`chr20` ↔ `20`, `chrM` ↔ `MT`), then an explicit alias table when the config
gives one (an NCBI `*_assembly_report.txt` or a UCSC `*.chromAlias.txt`).
Names that cannot be resolved are counted and reported, never silently
dropped; an annotation in which *no* interval matches a contig stops the run.
The contig names of a remote BAM are resolved the same way before slicing.

## Coordinate conventions

| Source file | Convention | Handled in |
|---|---|---|
| RepeatMasker `.out` (`hg38.fa.out.gz`, `hs1.repeatMasker.out.gz`) | 1-based, end inclusive | `annot.parse_rmsk_out` |
| UCSC table `genomicSuperDups.txt.gz` | 0-based, end exclusive, leading `bin` column | `annot.parse_ucsc_superdups` |
| BED (GIAB regions, CHM13 SDs) | 0-based, end exclusive | `annot.parse_bed` |
| GFF3 (RefSeq features) | 1-based, end inclusive | `annot.parse_gff` |
| VCF (GIAB truth set) | 1-based POS of the first REF base; indels carry an anchor base | `coords.vcf_to_bed`, `variants.py` |
| BAM via pysam | 0-based, end exclusive | used as is |
| samtools region strings (`chr20:1-1000`) | 1-based, end inclusive | `coords.region_string` |

Internally everything is 0-based, half-open. The conversions are unit-tested
(`tests/test_coords.py`), including a real RepeatMasker line. Two further
guards catch version mix-ups: an annotation interval or a variant that lies
beyond the end of its contig aborts the run with "different assembly".

## Reference genomes

| Key | Assembly | Used sequences | Annotations |
|---|---|---|---|
| `grch38` | GRCh38 / hg38 | chr20, chr21, chr22 | RepeatMasker, segmental duplications, GIAB v4.2.1 HG002 regions and variants |
| `chm13` | T2T-CHM13v2.0 / hs1 | chr20, chr21, chr22 | RepeatMasker, segmental duplications |
| `ecoli` | ASM584v2 | all | GFF features: insertion sequences and prophages (`mobile_genetic_element`), rRNA, tRNA, REP elements |
| `yeast` | R64 | all | GFF features: Ty elements, LTRs, rRNA, tRNA, telomeres |

Why three human chromosomes and not the genome: see the comment in
`config/config.yaml` and the limitations section of `docs/methods.md`.

## Real reads

All three read sets are HG002 (NA24385), from GIAB's aligned BAM files:

| Key | Platform | GIAB data set | How it is used |
|---|---|---|---|
| `hg002_illumina` | Illumina HiSeq, 2×148 bp, 300× | `NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam` | windows of chr20–22, back to FASTQ |
| `hg002_hifi` | PacBio Sequel II HiFi, 15 kb + 20 kb libraries | `PacBio_CCS_15kb_20kb_chemistry2/GRCh38/…haplotag.10x.bam` | same |
| `hg002_ont` | ONT ultra-long, Guppy 3.2.4 | `Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/…phased.bam` | same |

Retrieval is by **remote slicing**: `samtools view` reads the BAM index over
HTTPS and downloads only the blocks overlapping our windows (40 × 20 kb for
Illumina, 30 × 100–200 kb for long reads), a few hundred megabytes instead of
hundreds of gigabytes. The slice is converted back to raw FASTQ, so every
aligner starts from the same unaligned reads. GIAB's own placement of each
read is kept in a side table and used as an independent opinion in task 8.

Known consequences, all to be stated in the report:

- The reads were selected because GIAB's pipeline (novoalign, pbmm2 or
  minimap2 against the whole genome) placed them in our windows. Reads that
  pipeline misplaced elsewhere are missing; reads it misplaced into the
  windows are included.
- The checksum GIAB publishes covers the whole BAM, so a slice cannot be
  verified against it; the expected MD5 is recorded for reference.
- The three data sets date from 2015–2020. The nanopore set in particular
  (R9.4 pores, Guppy 3) has a far higher error rate than current chemistry, so
  the "ONT" conclusions describe that generation. Newer data can be swapped in
  by editing `real_reads` (a `bam_slice` URL, an ENA run via `source: ena`, or
  local FASTQ files via `source: local`).

### Adding a read set from ENA

```yaml
real_reads:
  my_run:
    genome: ecoli
    platform: ont
    source: ena
    accession: <run accession>   # SRR…, ERR… or DRR…
    n_reads: 20000
    skip: 20000                  # skip the first reads of the file (flow-cell edge)
```

The workflow asks the ENA portal API for the run's FASTQ URLs, MD5 and
metadata, checks that the archive's `instrument_platform` and
`library_layout` agree with the platform you assigned (a mismatch is written
to the provenance table and printed), and streams only the requested reads.
ENA, NCBI SRA and DDBJ hold the same runs; ENA is used because it serves
FASTQ directly over HTTPS.
