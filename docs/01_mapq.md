# Task 1. What a mapping quality is supposed to mean

## The definition

An aligner answers the question "where in the reference did this read come from?".
The mapping quality (MAPQ, column 5 of a SAM record) is the aligner's stated
confidence in that answer, written on the Phred scale:

```
MAPQ = -10 * log10( P(the reported position is wrong) )
```

| MAPQ | claimed probability of a wrong position | in words |
|---|---|---|
| 0 | 1 (in practice: 0.5 or more) | "I picked one of several equally good places" |
| 3 | 0.5 | two equally good places |
| 10 | 0.1 | 1 in 10 |
| 20 | 0.01 | 1 in 100 |
| 30 | 0.001 | 1 in 1,000 |
| 60 | 0.000001 | 1 in a million |

The scale is the same one used for base qualities, and it was introduced for
alignments by Li, Ruan and Durbin (2008) in the MAQ paper. The SAM
specification fixes the meaning ("MAPping Quality. It equals −10 log10
Pr{mapping position is wrong}, rounded to the nearest integer. A value 255
indicates that the mapping quality is not available"), but not how to compute it.

## Where the probability comes from

The probability is a Bayesian posterior. Suppose read `r` could have come from
candidate positions `u1, u2, …`. For each candidate we can compute a likelihood
`P(r | u)`: the probability of seeing exactly this read if it really came from
there. With base error rate `e`, a matching base contributes `(1 − e)` and a
mismatching base `e / 3` (the base was misread, and into this particular one of
three wrong bases). If every position is equally likely beforehand,

```
P(u_best is right | r) = P(r | u_best) / Σ_u P(r | u)
```

and the MAPQ is the Phred score of one minus that.

Two consequences follow directly, and both are worth knowing by heart
(`results/tables/mapq_reference.tsv` reproduces them):

1. **MAPQ measures ambiguity, not read quality.** A read with five mismatches
   in a unique region deserves a high MAPQ; a perfect read in a two-copy repeat
   deserves MAPQ 3. What matters is how much better the best position explains
   the read than the second best.
2. **Each extra mismatch at the competing position is worth a fixed number of
   MAPQ points.** One additional mismatch multiplies the likelihood by
   `(e/3) / (1 − e)`. At `e = 1 %` that is a factor of about 300, i.e. about 25
   MAPQ points per distinguishing mismatch; at `e = 0.1 %`, about 35 points.

`alnfail/mapq_theory.py` checks this with the smallest possible genome: a locus
and one paralogous copy that differs at `d` sites. Reads are drawn from the
true copy, assigned to the copy that explains them better, and labelled with
the exact posterior. Counting how often reads with each MAPQ were actually
wrong puts every point on the diagonal (`fig_mapq_theory.png`). **An aligner
that reports the true posterior is perfectly calibrated by construction.** That
diagonal is the yardstick for task 6.

## Why real aligners cannot report the true posterior

The sum in the denominator runs over every position in the genome. A real
aligner never evaluates it:

- it only looks at positions that share a seed (a k-mer, a minimizer, a
  maximal exact match) with the read, so candidates without a seed hit are
  invisible;
- it stops extending and scoring candidates as soon as heuristics say enough;
- its alignment score is not a log-likelihood under the true error model (the
  scoring matrix is fixed, base qualities are mostly ignored, indels are
  penalised by convention rather than by their real frequency).

So every aligner replaces the posterior by a **heuristic that grows with the
gap between the best and the second-best alignment** and is then squeezed
into a fixed range. They differ in what they measure the gap on, how they scale
it, and where they cap it. That is why the same read can get MAPQ 60 from one
tool and 42 from another, and why MAPQ values cannot be compared across tools
without a calibration like ours.

| Aligner | What the MAPQ is computed from | Range | Notes |
|---|---|---|---|
| MAQ (historical) | sum of base qualities at mismatches of the best hit versus the second-best hit, plus the number of second-best hits | 0–99 | the original approximation of the posterior |
| BWA-MEM / BWA-MEM2 | Smith-Waterman score of the best hit minus the best sub-optimal score, scaled (≈ 6.02 × score gap ÷ match score), reduced for many sub-optimal hits and for repetitive seeds | 0–60 | BWA-MEM2 is a faster re-implementation and reports the same values |
| Bowtie 2 | how far the best alignment score is above the minimum acceptable score, and how far above the second-best, looked up in a fixed table | 0–42, only certain values occur | not a formula but a lookup table; 42 is its "unique" |
| minimap2 | chaining scores: `40 · (1 − f2/f1) · min(1, m/10) · ln f1`, with `f1`, `f2` the best and second-best chain scores and `m` the number of seeds in the best chain | 0–60 | computed before base-level alignment; few seeds cap the confidence |
| Winnowmap2 | same code base and formula as minimap2, but seeds are chosen with repetitive k-mers down-weighted, which changes which chains exist | 0–60 | designed for satellite and duplicated sequence |
| strobealign | a minimap2-style formula on the scores of its best and second-best seed matches; for pairs, the gap between the two best pair scores | 0–60 | |
| NGMLR | relative gap between the best and second-best alignment score of each read segment, scaled to 0–60 | 0–60 | split reads get one MAPQ per segment |

These summaries come from the papers and source code listed below. Before the
defence, open the cited source for the aligners you talk about: the exact
constants are in `bwamem.c` (`mem_approx_mapq_se`), Bowtie 2's `unique.h`
(class `BowtieMapq2`), and section 2.1.3 of the minimap2 paper.

Further reasons the numbers differ:

- **Paired-end rescue.** For read pairs, the mate's position is evidence. A
  read that is ambiguous alone can be placed confidently if its mate is
  unique. Each tool folds this in differently.
- **Caps.** 60 (BWA, minimap2) and 42 (Bowtie 2) are ceilings, not
  probabilities. "MAPQ 60" does not mean the tool computed one in a million;
  it means "as confident as I ever get".
- **MAPQ 0 conventions.** BWA and minimap2 give 0 to a read with two equal
  best hits and report one at random. Bowtie 2 gives 0 or 1.

## What MAPQ can never account for

The posterior is taken over positions *in the reference the aligner was
given*. If the read's true origin is missing from that reference (an
unassembled centromere, a paralog on a chromosome that was left out, a
sequence present in the sample but not in the reference), the best available
position can look unique and receive a high MAPQ while being wrong. This is
the mechanism behind "blind spots": the aligner is confidently wrong because
the alternative was never on the table. It is why this project runs the same
reads against GRCh38 and against T2T-CHM13, and why restricting the human
reference to three chromosomes is listed as a limitation.

## How task 6 tests the claim

For every aligner and every MAPQ value `q` we count, on simulated reads whose
origin is known, the alignments reported with that MAPQ (`n`) and how many of
them are wrong (`k`). The claim is `k / n ≤ 10^(−q/10)`. We report `k / n` with
a 95 % Wilson interval and convert it back to the Phred scale ("MAPQ earned").

The headline question, *does MAPQ 30 mean one in a thousand?*, has a
sample-size trap. To show that an error rate is at most 1 in 1,000 you need
well over a thousand alignments at that MAPQ; with zero errors in `n`
alignments the 95 % upper bound is about `3 / n`. That is why the tables carry
intervals and why "no errors seen" is drawn as a lower bound, not as "perfect".

## Sources

- Li H, Ruan J, Durbin R. Mapping short DNA sequencing reads and calling variants using mapping quality scores. *Genome Research* 2008.
- The SAM/BAM Format Specification Working Group. Sequence Alignment/Map Format Specification.
- Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997, 2013. Source: `bwamem.c`.
- Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *IPDPS* 2019.
- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nature Methods* 2012. Source: `unique.h`.
- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 2018.
- Jain C, Rhie A, Hansen NF, Koren S, Phillippy AM. Long-read mapping to repetitive reference sequences using Winnowmap2. *Nature Methods* 2022.
- Sahlin K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. *Genome Biology* 2022.
- Sedlazeck FJ et al. Accurate detection of complex structural variations using single-molecule sequencing. *Nature Methods* 2018 (NGMLR).
