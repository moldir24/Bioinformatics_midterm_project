# Tasks 4-8: score alignments, stratify failures, validate on real reads, report.


def strata_dir(wc):
    return f"{OUT}/resources/{DS[wc.ds]['genome']}/strata"


def long_flag(wc):
    return "" if is_short(wc.ds) else "--long"


rule bam_table:
    """Primary alignments -> columnar table."""
    input:
        OUT + "/aln/{ds}/{aligner}.bam",
    output:
        OUT + "/tables/alignments/{ds}/{aligner}.parquet",
    shell:
        "{ALNFAIL} bam-table --bam {input} --out {output}"


rule annotate_truth:
    """Label every simulated read with the genomic context of its true origin."""
    input:
        truth=OUT + "/data/sim/{ds}/truth.parquet",
        strata=strata_dir,
        fasta=lambda wc: genome_fasta(DS[wc.ds]["genome"]),
        vcf=lambda wc: genome_vcf(DS[wc.ds]["genome"]),
    output:
        OUT + "/data/sim/{ds}/truth.annotated.parquet",
    params:
        vcf=lambda wc: vcf_flag(DS[wc.ds]["genome"]),
        long=long_flag,
        divergence=lambda wc: DS[wc.ds]["divergence"],
    shell:
        "{ALNFAIL} annotate-truth --truth {input.truth} --strata {input.strata} --fasta {input.fasta} "
        "{params.vcf} {params.long} --divergence {params.divergence} --out {output}"


rule evaluate:
    """Correct, misplaced or unmapped, for every read."""
    input:
        truth=OUT + "/data/sim/{ds}/truth.parquet",
        table=OUT + "/tables/alignments/{ds}/{aligner}.parquet",
    output:
        OUT + "/tables/outcomes/{ds}/{aligner}.parquet",
    params:
        min_overlap=EVAL.get("min_overlap", 0.1),
    shell:
        "{ALNFAIL} evaluate --truth {input.truth} --table {input.table} --min-overlap {params.min_overlap} --out {output}"


def meta(wc):
    info = DS[wc.ds]
    return (f"dataset={wc.ds} kind={info['kind']} genome={info['genome']} platform={info['platform']} "
            f"divergence={fmt_div(info['divergence'])}")


rule summarise_sim:
    """Accuracy, MAPQ calibration and stratified failures for one simulated dataset."""
    input:
        truth=OUT + "/data/sim/{ds}/truth.annotated.parquet",
        evals=lambda wc: expand(OUT + "/tables/outcomes/{ds}/{aligner}.parquet", ds=wc.ds, aligner=aligners_of(wc.ds)),
        tables=lambda wc: expand(OUT + "/tables/alignments/{ds}/{aligner}.parquet", ds=wc.ds, aligner=aligners_of(wc.ds)),
        strata=strata_dir,
        fasta=lambda wc: genome_fasta(DS[wc.ds]["genome"]),
        vcf=lambda wc: genome_vcf(DS[wc.ds]["genome"]),
    output:
        directory(OUT + "/datasets/{ds}"),
    wildcard_constraints:
        ds=constraint(SIM_DS),
    params:
        evals=lambda wc: [f"{a}={OUT}/tables/outcomes/{wc.ds}/{a}.parquet" for a in aligners_of(wc.ds)],
        tables=lambda wc: [f"{a}={OUT}/tables/alignments/{wc.ds}/{a}.parquet" for a in aligners_of(wc.ds)],
        vcf=lambda wc: vcf_flag(DS[wc.ds]["genome"]),
        long=long_flag,
        meta=meta,
        min_overlap=EVAL.get("min_overlap", 0.1),
    shell:
        "{ALNFAIL} summarise --truth {input.truth} --evals {params.evals} --tables {params.tables} "
        "--strata {input.strata} --fasta {input.fasta} {params.vcf} {params.long} --meta {params.meta} "
        "--min-overlap {params.min_overlap} --outdir {output}"


rule summarise_real:
    """Task 8 on one real dataset: aligner consensus, strata, truth-free observables."""
    input:
        tables=lambda wc: expand(OUT + "/tables/alignments/{ds}/{aligner}.parquet", ds=wc.ds, aligner=aligners_of(wc.ds)),
        strata=strata_dir,
        fasta=lambda wc: genome_fasta(DS[wc.ds]["genome"]),
        vcf=lambda wc: genome_vcf(DS[wc.ds]["genome"]),
        provider=lambda wc: f"{OUT}/data/real/{DS[wc.ds]['rd']}/raw.provider.parquet",
    output:
        directory(OUT + "/datasets/{ds}"),
    wildcard_constraints:
        ds=constraint(REAL_DS),
    params:
        tables=lambda wc: [f"{a}={OUT}/tables/alignments/{wc.ds}/{a}.parquet" for a in aligners_of(wc.ds)],
        vcf=lambda wc: vcf_flag(DS[wc.ds]["genome"]),
        long=long_flag,
        meta=meta,
        # the provider aligned to the read set's own genome; on another assembly the coordinates differ
        provider=lambda wc, input: f"--provider {input.provider}" if DS[wc.ds]["genome"] == REAL[DS[wc.ds]["rd"]]["genome"] else "",
        min_overlap=EVAL.get("min_overlap", 0.1),
    shell:
        "{ALNFAIL} real --tables {params.tables} --strata {input.strata} --fasta {input.fasta} {params.vcf} "
        "{params.long} {params.provider} --meta {params.meta} --min-overlap {params.min_overlap} --outdir {output}"


rule collect:
    """Merge all datasets into the result tables, including the simulation-versus-real comparison."""
    input:
        datasets=expand(OUT + "/datasets/{ds}", ds=SIM_DS + REAL_DS),
        calibration=expand(OUT + "/calibration/{rd}.tsv", rd=CALIBRATORS),
        qc=expand(OUT + "/qc/{rd}.tsv", rd=list(REAL)),
        strata=expand(OUT + "/resources/{genome}/strata", genome=list(GENOMES)),
    output:
        headline=f"{OUT}/tables/headline.tsv",
        accuracy=f"{OUT}/tables/accuracy.tsv",
        strata=f"{OUT}/tables/strata.tsv",
        calibration=f"{OUT}/tables/mapq_calibration.tsv",
        meta=f"{OUT}/tables/report_meta.json",
    params:
        pairs=transfer_pairs(),
        strata_stats=lambda wc, input: [f"{d}/strata_stats.tsv" for d in input.strata],
        min_reads=EVAL.get("min_reads_per_stratum", 100),
        outdir=f"{OUT}/tables",
        meta=json.dumps({
            "labels": labels(),
            "profiles": [{"dataset": rd, "platform": REAL[rd]["platform"], "real": f"{OUT}/profiles/{rd}.raw.npz",
                          "sim": f"{OUT}/profiles/{rd}.sim.npz"} for rd in CALIBRATORS],
        }),
    shell:
        "{ALNFAIL} collect --datasets {input.datasets} --pairs {params.pairs} --calibration {input.calibration} "
        "--qc {input.qc} --strata-stats {params.strata_stats} --min-reads {params.min_reads} --outdir {params.outdir} && "
        "printf '%s' {params.meta:q} > {output.meta}"


rule benchmarks:
    """Runtime and peak memory of every timed step, in one table."""
    input:
        datasets=expand(OUT + "/datasets/{ds}", ds=SIM_DS + REAL_DS),
    output:
        f"{OUT}/tables/benchmarks.tsv",
    params:
        pattern=f"{OUT}/benchmarks",
    shell:
        "{ALNFAIL} benchmarks --files $(find {params.pattern} -name '*.tsv' | sort) --out {output}"


rule mapq_theory:
    """Task 1: the exact posterior as a yardstick."""
    output:
        f"{OUT}/tables/mapq_ideal_calibration.tsv",
        f"{OUT}/tables/mapq_reference.tsv",
    params:
        outdir=f"{OUT}/tables",
    shell:
        "{ALNFAIL} mapq-theory --outdir {params.outdir}"


rule report:
    """Figures, the results digest (summary.md) and the same digest as one readable page (report.html)."""
    input:
        f"{OUT}/tables/headline.tsv",
        f"{OUT}/tables/mapq_ideal_calibration.tsv",
        meta=f"{OUT}/tables/report_meta.json",
        benchmarks=f"{OUT}/tables/benchmarks.tsv",
    output:
        f"{OUT}/report/summary.md",
        f"{OUT}/report/report.html",
    params:
        tables=f"{OUT}/tables",
        outdir=f"{OUT}/report",
    shell:
        "{ALNFAIL} report --tables {params.tables} --outdir {params.outdir} --meta {input.meta} --benchmarks {input.benchmarks}"
