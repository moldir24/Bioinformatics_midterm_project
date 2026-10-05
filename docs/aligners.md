# The aligners: what they do and why they are run this way

Seven tools, two groups. The split is not a matter of taste: short-read and
long-read aligners solve differently shaped problems.

| | Short reads | Long reads |
|---|---|---|
| Read | 100–250 bp, error ≈ 0.1–1 %, mostly substitutions | 1–100+ kb, error 0.1–15 %, mostly indels |
| Seeds | a short read has few seeds; all of them count | thousands of seeds; most can be thrown away |
| Hard part | too little information: a 150 bp read fits many repeats | too much noise: indels break exact matches |
| Core idea | exact seeds (FM-index or hashes) then full dynamic programming | sparse seeds (minimizers), *chaining* of collinear seeds, banded alignment between them |
| Helped by | the mate: pairing resolves many ambiguities | length: a read that spans a repeat anchors itself in unique flanks |

## Short-read aligners

### BWA-MEM2 (`bwamem2`), defaults
Seeds are super-maximal exact matches found in an FM-index (a compressed
suffix array of the reference). Seeds are chained, extended by banded
Smith-Waterman, and the mate is "rescued" by a local search near its partner.
BWA-MEM2 is a vectorised re-implementation of BWA-MEM producing the same
alignments, and is the de-facto standard for human short-read data, so it is
the baseline every other tool is compared against. Defaults are used because
that is how almost everyone runs it; it estimates the fragment-length
distribution from the data itself.

Memory: building the index needs about 28 bytes per reference base (about
5 GB for chr20–22, about 85 GB for a whole human genome). On a whole genome
with less memory, put `bwa` in `aligners.short` instead: same algorithm,
about 5 GB.

### Bowtie 2 (`bowtie2`), defaults plus `-X` from the data
Also FM-index based, but seeds are fixed-length substrings sampled at
intervals, extended by SIMD dynamic programming in end-to-end mode (the whole
read must align; no soft clipping). The one option we set is `-X`, the largest
fragment length considered a valid pair. Its default is 500 bp. Our measured
fragment distribution is used instead (`fragment_max` in the profile), since a
real library with longer fragments would otherwise have its long pairs
reported as discordant, which costs the aligner the pairing information.

### minimap2 `-ax sr`
The long-read aligner in its short-read preset: minimizer seeds (k = 21,
w = 11), chaining, then base-level alignment. Included because it is widely
used as a fast alternative and because the same tool appears in the long-read
group, which lets us see how one algorithm behaves on both kinds of data.

### strobealign (`strobealign`), defaults
Seeds are *strobemers*: two k-mers linked at a variable distance, which stay
matchable when an indel or mismatch falls between them. Very fast, index built
on the fly. It picks its seed parameters from the read length it observes.

## Long-read aligners

### minimap2 `-ax map-ont` / `-ax map-hifi`
Minimizer index, chaining by dynamic programming over collinear seeds, banded
alignment in the gaps. The presets differ in seed size and scoring:
`map-ont` (k = 15) tolerates the higher, indel-rich error of nanopore reads;
`map-hifi` (k = 19, stricter gap penalties) exploits the accuracy of HiFi
reads to separate near-identical repeat copies. Using the wrong preset is a
classic mistake; the preset is therefore taken from the platform in the config.
(With minimap2 ≥ 2.27 and modern Q20+ nanopore data, `-ax lr:hq` is the
recommended preset; our ONT set is older R9.4 data, for which `map-ont` is
right.)

### Winnowmap2 (`winnowmap`), `-ax map-ont` / `-ax map-pb`
A minimap2 derivative for repetitive sequence. Standard minimizer sampling
discards the most frequent k-mers, which in a satellite array means discarding
almost all seeds. Winnowmap instead *down-weights* the most frequent 0.02 % of
15-mers (counted with meryl; this is the `index_winnowmap` rule) so that they
are sampled less but not lost. Expectation to test: better placement in
satellites and segmental duplications, at some cost in speed.

### NGMLR (`ngmlr`), `-x ont` / `-x pacbio`
Built for structural-variant calling: it splits a read into sub-segments,
places them independently and joins them with a convex gap cost, so that a
read spanning a large deletion or inversion is reported as a split alignment
rather than forced into one. It is much slower, and it does not print reads it
cannot place (our evaluation counts reads missing from the output as unmapped).

### BWA-MEM `-x ont2d` / `-x pacbio` (`bwa`)
The short-read algorithm with relaxed seed length and gap penalties. This is
how long reads were aligned before minimizer-based tools existed. It is in the
panel as the historical baseline: it shows what a seed-and-extend design costs
on indel-rich reads.

## What is deliberately the same for every tool

- identical reads: every aligner of a group receives the same FASTQ files;
- identical thread count (`threads` in the config), timed with Snakemake's
  `benchmark` (wall time, CPU time, peak resident memory);
- identical post-processing: output converted to BAM in input order, primary
  alignment taken as the aligner's answer;
- no tool-specific tuning beyond the documented preset. A benchmark that tunes
  one tool and not the others measures the tuner.

## Reading the MAPQ of each tool

See `docs/01_mapq.md`. In short: BWA, minimap2, Winnowmap, strobealign and
NGMLR report 0–60; Bowtie 2 reports 0–42 from a lookup table. A filter such as
"MAPQ ≥ 30" therefore selects different things for different tools, which is
exactly what `fig_mapq_threshold.png` shows.

## Versions

The versions actually used are pinned in `environment.yml` and recorded in the
logs under `results/logs/`. Quote them in the report: aligner behaviour,
including MAPQ, changes between releases.
