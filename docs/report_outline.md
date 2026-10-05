# Report outline (8 to 12 pages)

The handbook asks for a report "written for a scientist who is not a
programmer", structured as: problem and biological framing, data sources with
every accession, methods, results, limitations, conclusion. Below, each section
lists what belongs in it and which file supplies the numbers. Write the text
yourselves, from your own results: the digest (`results/report/summary.md`)
gives the numbers, not the interpretation.

Page budget is a suggestion. Figures need axis labels, units and captions
(they have labels; you write the captions).

## 1. Problem and biological framing (about 1 page) — rubric 1

- What read alignment is for: every downstream result (variants, expression,
  assembly polishing) inherits the aligner's mistakes.
- The two ways to fail and why they are not equal: an unmapped read is lost
  information, a misplaced read is false evidence (a false variant call in a
  paralog, wrong coverage in a duplicated gene).
- Why genomes make this hard: repeats, segmental duplications, satellites; the
  share of the human genome they occupy; medically relevant genes inside them.
- Why an independent benchmark: authors benchmark their own tools; simulated
  reads give ground truth; the simulation must be shown to be realistic.
- The question, stated precisely: for each aligner, in which genomic contexts
  does it misplace reads, and is its reported confidence (MAPQ) honest?

## 2. Data sources (about 1 page) — rubric 2

- Table from `results/tables/provenance.tsv`: every URL or accession, size,
  checksum, retrieval date.
- The four archives and what each contributes (`docs/data_sources.md`).
- Identifier mapping and coordinate conventions: one paragraph plus the small
  table of conventions; mention the two guards (contig-name resolution report,
  "different assembly" check) and what `strata_sources.tsv` showed.
- Any problem you met and how you dealt with it (a moved URL, a contig naming
  mismatch, the unverifiable checksum of a BAM slice). The rubric rewards
  problems "found, documented and dealt with".

## 3. Platforms (about 1 page) — rubrics 3, 4, 5

- One table, measured on your real reads (`results/profiles/*.raw.json`):
  read length (mean, N50), error rate, split into substitutions / insertions /
  deletions, homopolymer enrichment.
- What each instrument's chemistry explains about those numbers, and what each
  makes possible or impossible (a 150 bp read cannot span an Alu; a 20 kb read
  spans most repeats but not a satellite array).
- QC: the decision table (`results/tables/qc_decisions.tsv`): each filter, its
  threshold, how much it removed and why it exists. Short-read artefacts
  (adapter read-through, poly-G) and why read ends were *not* quality-trimmed.

## 4. Methods (about 2 pages) — rubrics 6, 7, 9

- Figure: the pipeline as a flow diagram (draw it yourselves).
- What MAPQ is supposed to mean (`docs/01_mapq.md`, condensed to half a page,
  with `fig_mapq_theory.png`).
- Calibration loop and realism check (`fig_calibration_*.png`,
  `calibration_check.tsv`).
- The aligners and why each is in the panel (`docs/aligners.md`), versions.
- Scoring rule, with the reason for the 10 % overlap criterion.
- Strata definitions.
- Engineering in a few sentences: streaming, Parquet, remote slicing,
  Snakemake, measured runtime and memory (`fig_resources.png`).

## 5. Results (3 to 4 pages) — rubric 8

Open with the headline, then justify it.

1. **Accuracy** (`fig_accuracy.png`, `headline.tsv`): correct / misplaced /
   unmapped per aligner and platform. Which tools trade misplacement for
   unmapped reads?
2. **Is MAPQ honest?** (`fig_mapq_calibration.png`, `fig_mapq_threshold.png`,
   digest section 2). Answer the handbook's question directly: for each tool,
   the observed error at MAPQ ≥ 30 with its interval, against 1 in 1,000.
3. **Where each tool breaks** (`fig_strata_*.png`, `strata.tsv`): the hardest
   contexts, the contexts where tools differ most, and the
   confident-but-wrong cells (misplaced at MAPQ ≥ 30). These are the "blind
   spots nobody advertises".
4. **GRCh38 versus CHM13** and **bacterium / yeast versus human**: what changes
   when the reference is complete, and when the genome is small.
5. **Divergence** (`fig_divergence.png`): how fast each tool degrades.
6. **Short versus long reads**: the same contexts, a different outcome; say
   why (read length versus repeat length).
7. **Does it hold on real data?** (`fig_transfer.png`, digest section 5):
   what transfers, what does not, and the likely reason for each failure.
8. **Negative results**: anything you expected and did not find. Report it.

## 6. Limitations (about 1 page) — read closely by the markers

Start from the list in `docs/methods.md` and rewrite it for *your* results:
which limitations actually bite, and how they would change the conclusions.
Add anything you observed that the list does not cover.

## 7. Conclusion (about half a page)

Advice a non-specialist can act on: which aligner for which data, which MAPQ
threshold means what for which tool, which regions to distrust regardless.

## Contribution statement (one paragraph, required)

Name who led which component. It will be checked against the Git history, so
commit your own work under your own account as you go.

## Appendix: AI assistance (required)

State which tools were used and for what. For this repository, for example:
the pipeline code, documentation and unit tests were drafted with an AI
assistant (Claude) and then reviewed, run, debugged and modified by the
authors; list what you changed, what you verified and how. The handbook's rule
applies: you are responsible for every line, and either of you must be able to
explain any part of it.
