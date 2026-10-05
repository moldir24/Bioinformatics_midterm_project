# Methods, design decisions and limitations

This file explains *why* the pipeline does what it does. It is the raw
material for the Methods and Limitations sections of the report, and the list
of things either partner may be asked about at the defence.

## The experiment in one paragraph

Real reads from three platforms are used to measure each platform's error
profile. A simulator driven by those profiles produces reads whose reference
origin is known. Seven aligners place those reads; each placement is scored as
correct, misplaced or unmapped, and each reported mapping quality is compared
with the error rate actually observed. Failures are then broken down by the
genomic context of the read's true origin. Finally, the held-out half of the
real reads is aligned by the same tools, and what can be measured without
truth is compared between simulation and reality.

## Why a simulator of our own

Standard simulators (ART, Mason, PBSIM, Badread) exist and are good. We wrote
a small one for three reasons that matter for this project:

1. **Truth stays exact when the sample differs from the reference.** We
   simulate reads from an individual (the reference with GIAB truth variants
   applied) and from diverged genomes. Our simulator carries the reference
   coordinate of every base through every edit, so the true interval is read
   off directly. With an external simulator one has to build a modified genome
   and lift coordinates back, which is where truth sets usually go wrong.
2. **One procedure measures real and simulated reads.** The profile learned
   from real reads and the profile of the simulated reads come from the same
   code, so "does the simulation match?" is a like-for-like comparison.
3. **Every line can be explained.** Under 500 lines of numpy in
   `alnfail/simulate.py`.

## Task 3: calibration is a loop

`learn-profile` measures, on real reads aligned to their reference:

- short reads: read lengths, a per-cycle Markov chain of base qualities, the
  mismatch probability given quality and position in the read, the
  substitution spectrum, indel rates and lengths, fragment lengths;
- long reads: the joint sample of (read length, per-read error rate), the
  split of errors into substitutions, insertions and deletions, indel length
  distributions, and the enrichment of indels in homopolymers by run length.

Known variant sites of the individual (GIAB truth set, ±2 bp) are masked so
that true variants are not counted as sequencing errors. Only primary
alignments with MAPQ ≥ 20 are used.

An error rate *measured through an alignment* is not the rate at which errors
were *made*: the aligner turns an adjacent insertion and deletion into a
mismatch, moves indels within homopolymers, clips bad read ends. Injecting the
measured rates therefore does not reproduce them. So calibration closes the
loop: a pilot simulation with the measured rates is aligned and measured by
the same procedure, and each error type is rescaled by real ÷ pilot
(`calibrate`). The benchmark simulations use the corrected profile, and the
final check (`compare-profiles`, `results/tables/calibration_check.tsv`,
`fig_calibration_*.png`) compares real and simulated profiles metric by
metric, flagging everything outside the tolerance.

## Tasks 4 and 5: what is simulated

- Windows are drawn uniformly over the non-N sequence of the chosen contigs.
- On GRCh38 each read is drawn from one of the two HG002 haplotypes (GIAB
  variants applied; unphased heterozygous sites are assigned to a haplotype by
  a hash of their position, so runs are reproducible).
- Optional divergence adds substitutions and short indels at a set rate (90 %
  substitutions, 10 % indels), modelling reads from a related genome.
- Short reads: paired-end, fragment and read lengths from the real
  distributions, mates on opposite strands, errors from the quality chain.
- Long reads: length and error rate drawn together from one real read; errors
  are indel-dominated and homopolymer-dependent; half of inserted bases repeat
  the preceding base.
- Every dataset is seeded by `seed` and its own id: same config, same reads.

## Task 6: scoring

- Only the **primary** alignment of each read is scored. Secondary records are
  alternative guesses; supplementary records are other pieces of a split read.
- **Correct**: same contig, and the aligned reference interval overlaps the
  true interval by more than 10 % of their union. This is the criterion of
  minimap2's `paftools.js mapeval`, chosen so that our numbers are comparable
  with published benchmarks. It tolerates soft clipping at noisy read ends.
  Stricter thresholds (50 %, 90 %) are kept in `accuracy.tsv` to show how much
  the result depends on the choice.
- **Misplaced**: mapped, but not correct. Subdivided into other contig, same
  contig elsewhere, same locus with low overlap (`misplaced_kinds.tsv`).
- **Unmapped**: flagged unmapped, or absent from the output.
- For pairs, each mate is scored separately.
- **MAPQ calibration**: per aligner and reported MAPQ, alignments and wrong
  alignments are counted; rates carry 95 % Wilson intervals. With zero errors
  observed, the upper interval bound is reported as a lower bound on the MAPQ
  earned.

## Task 7: strata

Each read is labelled by the context of its **true** origin:

| Stratum set | Source | Meaning |
|---|---|---|
| `repeat_class` | RepeatMasker | SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity, … (the class covering at least half the read) |
| `repeat_age` | RepeatMasker divergence column | young (< 5 % from consensus), middle, old; young copies are nearly identical to each other |
| `segdup` | segmental duplications | binned by identity between the copies (≥ 99 %, 95–99 %, < 95 %) |
| `feature` | GFF (bacteria, yeast) | insertion sequences, rRNA operons, tRNA, LTRs, telomeres |
| `giab_confident` | GIAB regions | inside or outside the regions GIAB considers characterisable |
| `gc_content` | the read's true reference interval | six bins from < 30 % to ≥ 70 % |
| `variant_density` | GIAB variants in the interval | per kilobase, so short and long reads share a scale |
| `read_length` | long reads | length bins |
| `*_cover` | all interval sets | whether the read only touches the element or lies entirely inside it |

A read spanning a repeat with unique flanks is easy; a read lying entirely
inside it is hard. The `_cover` strata separate those two cases, which matters
most for long reads.

"Divergence" appears twice, on purpose: divergence of a *repeat copy from its
family consensus* (`repeat_age`), and divergence of the *sample from the
reference* (the divergence series).

## Task 8: validation on real reads

Real reads have no known origin. Two kinds of quantity can still be compared:

1. **Observables**: mapped fraction, MAPQ distribution, edit distance per
   base, clipping, split alignments, per aligner. If the simulation is
   realistic, an aligner treats simulated reads as it treats real ones.
2. **Disagreement with the consensus**: the placement supported by a strict
   majority of aligners stands in for truth. On simulated reads, where truth
   is known, we first measure how well that stand-in tracks real error
   (`proxy_validity_spearman`). Then the per-stratum disagreement rates are
   compared between simulation and reality (`transfer_spearman`,
   `fig_transfer.png`).

GIAB's own whole-genome alignment of the same reads gives a third reference
(`provider_agreement.tsv`). The GIAB truth set contributes the variant-density
and confident-region strata on both sides.

## Engineering

- **Stream, do not load**: FASTQ and VCF are read record by record; reference
  windows are fetched through the FASTA index; remote BAMs are sliced through
  their index.
- **Columnar tables**: every per-read table is Parquet with zstd compression.
- **Statistical power instead of coverage**: accuracy needs a number of reads,
  not a depth. With `n` reads, rates down to about `3/n` can be told from zero.
- **Measured cost**: every simulate, index and map job is timed
  (`results/tables/benchmarks.tsv`, `fig_resources.png`).
- **Workflow**: Snakemake; one rule per step; rerunning continues where it
  stopped; any step can be run by hand via `python -m alnfail …`.

## Limitations (to be stated without being asked)

1. **Three chromosomes, not the genome.** Paralogs elsewhere are absent from
   the reference, so misplacement between chromosomes is under-estimated and
   MAPQ looks better than on a full genome. chr20–22 were chosen to keep seven
   aligners runnable on a laptop; `contigs: all` removes the limit.
2. **The simulator has no chimeric reads, no contamination, no GC or
   sequence-specific coverage bias, no systematic context-specific errors**
   beyond homopolymers. Real unmapped and split-read fractions will be higher
   than simulated ones; the transfer table quantifies this.
3. **Error profiles are measured through one aligner** (minimap2) on reads it
   places confidently, which favours easy regions. The pilot loop corrects the
   three overall rates (mismatch, insertion, deletion) but not this selection,
   and not the shape of the homopolymer enrichment.
4. **Long-read base qualities are simulated as one constant per read.** The
   long-read aligners tested do not use base qualities for placement, so this
   does not affect the benchmark, but the FASTQ is not suitable for testing
   quality-aware tools.
5. **Real reads were pre-selected by GIAB's aligner** (see
   `docs/data_sources.md`), and the platforms are represented by 2015–2020
   chemistry.
6. **Variants only for GRCh38.** CHM13, *E. coli* and yeast reads are
   simulated from the bare reference (plus optional uniform divergence).
   Uniform divergence is not how real genomes diverge: real differences
   cluster and include structural variants.
7. **Consensus is not truth.** Aligners sharing an algorithm (BWA-MEM2 and
   minimap2 both rely on exact seeds) can be wrong together; disagreement
   measures difference, not error.
8. **Default parameters.** Each tool can be tuned to do better in a given
   context; we benchmark how the tools are normally run.
9. **Rare events need many reads.** A MAPQ-60 error rate of one in a million
   cannot be confirmed with 400,000 reads; the intervals in the tables show
   what can and cannot be concluded.
