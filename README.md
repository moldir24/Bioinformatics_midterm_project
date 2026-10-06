# Where Do Aligners Actually Fail?

An independent benchmark of read aligners, built on reads whose true origin is
known. Project 26 of the AITU *Introduction to Bioinformatics* course project.

Aligner benchmarks are usually published by the authors of the aligner. This
one calibrates a read simulator on real Genome in a Bottle data, simulates
short and long reads from four genomes, runs seven aligners on identical
reads, and reports for each of them:

- how often a read is **correctly placed, misplaced, or left unmapped**;
- whether the reported **mapping quality is honest** (does MAPQ 30 mean one
  error in a thousand?);
- **where** it fails: by repeat class, repeat age, segmental duplication,
  GC content, variant density, read length and divergence;
- whether those conclusions **hold on real reads**.

## Quick start

Linux x86-64 (on Windows use WSL2; on Apple silicon use the Dockerfile).
Conda or mamba is the only prerequisite.

```bash
# 1. get the code
git clone <your repository URL> && cd <repository>

# 2. create the environment
mamba env create -f environment.yml        # or: conda env create -f environment.yml
conda activate alnfail

# 3. check the installation on the bundled toy dataset (a few minutes, no network)
python -m pytest -q
snakemake --cores 4 --configfile config/test.yaml

# 4. check that every real data source resolves from your network
snakemake check_sources --cores 1 --configfile config/config.yaml

# 5. run the real analysis: this one command reproduces the headline result
snakemake --cores 8 --configfile config/config.yaml
```

Step 5 downloads about 3 GB and needs about 8 GB of RAM. It has not been
timed on the real data yet. For orientation: the toy run takes six minutes on
two cores, and a 5 Mb bacterial genome with 200,000 read pairs and 4,000 to
6,000 long reads per dataset (nine simulated and three real datasets) took 25
minutes on two cores with a peak of 1 GB. The run can be interrupted and
restarted; Snakemake continues where it stopped. Open
`results/report/report.html` in a browser to read the results; the headline
table is `results/tables/headline.tsv`, the figures and the plain-text digest
are in `results/report/`.

If conda cannot be used, `docker build -t alnfail .` builds the same
environment in a container (see the `Dockerfile` for the run command).

## What you get

```
results/
├── report/
│   ├── report.html                 everything on one page, figures included (start here)
│   ├── summary.md                  the same numbers as plain Markdown tables
│   ├── fig_accuracy.png            misplaced and unmapped reads per aligner, genome, platform
│   ├── fig_mapq_calibration.png    reported MAPQ against observed error
│   ├── fig_mapq_threshold.png      reads kept versus errors kept when filtering by MAPQ
│   ├── fig_strata_<genome>_<platform>.png   the failure map
│   ├── fig_divergence.png          accuracy as the sample diverges from the reference
│   ├── fig_calibration_<reads>.png simulator against real reads
│   ├── fig_transfer.png            simulation against real reads, per genomic context
│   ├── fig_resources.png           runtime and memory
│   └── fig_mapq_theory.png         what a calibrated aligner looks like
└── tables/
    ├── headline.tsv                one row per platform and aligner
    ├── accuracy.tsv                per dataset and aligner, with intervals
    ├── mapq_calibration*.tsv       per reported MAPQ value
    ├── strata.tsv                  per genomic context
    ├── misplaced_kinds.tsv         where misplaced reads went
    ├── calibration_check.tsv       is the simulator realistic?
    ├── transfer_*.tsv, consensus.tsv, observables.tsv, provider_agreement.tsv   real-data validation
    ├── qc_decisions.tsv            every QC filter, its threshold and what it removed
    ├── provenance.tsv              every source: URL, size, checksum, date
    ├── strata_sources.tsv          how each annotation matched the reference
    └── benchmarks.tsv              wall time, CPU time and peak memory of every job
```

## The eight tasks and where they live

| Task | What | Code | Explained in |
|---|---|---|---|
| 1 | What MAPQ means and why aligners compute it differently | `alnfail/mapq_theory.py` | `docs/01_mapq.md` |
| 2 | References, repeat and variation annotations, coordinate conventions | `alnfail/fetch.py`, `coords.py`, `annot.py`, `variants.py`; `workflow/rules/data.smk` | `docs/data_sources.md` |
| 3 | Calibrate the simulator on real reads | `alnfail/qc.py`, `profile.py`; `workflow/rules/simulate.smk` | `docs/methods.md` |
| 4 | Short reads: simulate, three outcomes | `alnfail/simulate.py`, `evaluate.py` | `docs/methods.md` |
| 5 | Long reads, indel-dominated error, different aligners | `alnfail/simulate.py` | `docs/aligners.md` |
| 6 | At least four aligners on identical data; MAPQ calibration | `workflow/rules/align.smk`; `alnfail/samtable.py`, `summarise.py` | `docs/aligners.md`, `docs/01_mapq.md` |
| 7 | Stratify failures | `alnfail/annot.py`, `summarise.py` | `docs/methods.md` |
| 8 | Do the conclusions transfer to real GIAB data? | `alnfail/consensus.py` | `docs/methods.md` |

`docs/report_outline.md` maps tables and figures to the sections of the
written report.

## Repository layout

```
alnfail/            the analysis package; every step is `python -m alnfail <command>`
                    (`html.py` turns the digest and figures into the one-page report)
workflow/           Snakemake workflow (Snakefile and rules/, one file per group of tasks)
config/config.yaml  the real analysis: genomes, platforms, read sets, aligners
config/test.yaml    the same workflow on the toy dataset
test_data/          small synthetic dataset for verification (see its README)
tests/              unit tests: coordinates, simulator truth, scoring, calibration arithmetic, the HTML page
docs/               task 1 essay, data sources, aligners, methods and limitations, report outline
environment.yml     conda environment; Dockerfile builds the same in a container
```

## Changing the experiment

Everything is driven by `config/config.yaml`:

- **More statistical power**: raise `n_reads` under `platforms`. The smallest
  error rate distinguishable from zero is about 3 ÷ (number of reads).
- **Whole human genome**: set `contigs: all` for `grch38` and `chm13` (needs
  ≥ 64 GB RAM; replace `bwamem2` by `bwa` to stay near 16 GB for the others).
- **Another read set**: add it under `real_reads` (`bam_slice`, `ena` or
  `local`) and point a platform's `calibrate_on` at it.
- **Another aligner**: add an index and a map rule to
  `workflow/rules/align.smk`, its name to `aligners`, and a colour to
  `alnfail/plots.py`.

Useful commands:

```bash
snakemake -n --configfile config/config.yaml                 # dry run: list what would be done
snakemake --cores 8 --configfile config/config.yaml --rerun-incomplete
python -m alnfail --help                                     # every step can be run by hand
python -m alnfail html --report results/report               # rebuild report.html from summary.md and the figures
```

## Status of this code

- The whole workflow has been run end to end on the toy dataset and on a
  5 Mb bacterial genome with synthetic stand-in reads; the unit tests pass.
- **It has not yet been run on the real sources in `config/config.yaml`**
  (they were not reachable from the machine it was written on). Expect to fix
  small things on the first real run: a moved URL, an unexpected column in an
  annotation file. `check_sources` and the built-in format checks are there to
  make those failures loud and early.
- Developed against: Python 3.13, Snakemake 9.27, pandas 3.0, numpy 2.5,
  samtools 1.19.2, minimap2 2.26, BWA 0.7.17, BWA-MEM2 2.2.1, Bowtie 2 2.5.2,
  strobealign 0.16.1, Winnowmap 2.03, NGMLR 0.2.7, fastp 0.23.4.

## Acknowledgements and disclosure

Data: Genome in a Bottle Consortium (NIST), Telomere-to-Telomere Consortium,
UCSC Genome Browser, NCBI RefSeq. The pipeline code was drafted with an AI
assistant; see the AI-assistance appendix of the report for what was generated
and what the authors verified and changed.
