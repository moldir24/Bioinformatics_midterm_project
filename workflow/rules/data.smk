# Task 2: references, annotations, truth sets and real reads.
# Every file is fetched by a rule and leaves a provenance record.

rule check_sources:
    """Week 1: confirm that every configured source resolves, and how big it is."""
    output:
        report=f"{OUT}/tables/sources_check.tsv",
        spec=temp(f"{OUT}/tables/sources_check.json"),
    params:
        spec=json.dumps(sources_to_check()),
    shell:
        "printf '%s' {params.spec:q} > {output.spec} && "
        "{ALNFAIL} check-sources --sources {output.spec} --out {output.report}"


rule fetch_file:
    output:
        file=OUT + "/resources/{genome}/downloads/{fname}",
        provenance=OUT + "/provenance/{genome}__{fname}.json",
    params:
        url=lambda wc, output: DOWNLOADS[output.file]["url"],
        md5=lambda wc, output: DOWNLOADS[output.file]["md5"],
        label=lambda wc, output: DOWNLOADS[output.file]["label"],
    retries: 2
    shell:
        "{ALNFAIL} fetch-file --url {params.url:q} --out {output.file} --md5 '{params.md5}' "
        "--label {params.label:q} --provenance {output.provenance}"


rule prepare_reference:
    """Keep the configured contigs and index them."""
    input:
        fasta=lambda wc: FASTA_DL[wc.genome],
        alias=lambda wc: ALIAS_DL[wc.genome] or [],
    output:
        fasta=OUT + "/resources/{genome}/genome.fa",
        fai=OUT + "/resources/{genome}/genome.fa.fai",
        contigs=OUT + "/resources/{genome}/contigs.tsv",
    params:
        contigs=lambda wc: " ".join(contigs_of(wc.genome)),
        alias=lambda wc: f"--alias {ALIAS_DL[wc.genome]}" if ALIAS_DL[wc.genome] else "",
    shell:
        "{ALNFAIL} prepare-reference --fasta {input.fasta} --contigs {params.contigs} {params.alias} "
        "--out {output.fasta} --contigs-out {output.contigs}"


rule build_strata:
    """Repeat, duplication and feature annotations -> normalised strata (0-based BED)."""
    input:
        fasta=OUT + "/resources/{genome}/genome.fa",
        sources=lambda wc: [s["path"] for s in STRATA_SPEC[wc.genome]],
        alias=lambda wc: ALIAS_DL[wc.genome] or [],
    output:
        directory(OUT + "/resources/{genome}/strata"),
    params:
        spec=lambda wc: json.dumps(STRATA_SPEC[wc.genome]),
        alias=lambda wc: f"--alias {ALIAS_DL[wc.genome]}" if ALIAS_DL[wc.genome] else "",
    shell:
        "mkdir -p {output} && printf '%s' {params.spec:q} > {output}/sources.json && "
        "{ALNFAIL} build-strata --fasta {input.fasta} --spec {output}/sources.json {params.alias} --outdir {output}"


rule prepare_variants:
    """Truth-set VCF restricted to our contigs and renamed to our contig names."""
    input:
        vcf=lambda wc: VCF_DL[wc.genome],
        fasta=OUT + "/resources/{genome}/genome.fa",
        alias=lambda wc: ALIAS_DL[wc.genome] or [],
    output:
        vcf=OUT + "/resources/{genome}/variants.vcf.gz",
        stats=OUT + "/resources/{genome}/variants.json",
    params:
        alias=lambda wc: f"--alias {ALIAS_DL[wc.genome]}" if ALIAS_DL[wc.genome] else "",
    shell:
        "{ALNFAIL} prepare-variants --vcf {input.vcf} --fasta {input.fasta} {params.alias} "
        "--out {output.vcf} --stats {output.stats}"


def read_spec(wc):
    cfg = dict(REAL[wc.rd])
    cfg["paired"] = PLATFORMS[cfg["platform"]]["type"] == "short"
    cfg["label"] = wc.rd
    cfg["contigs"] = [c for c in contigs_of(cfg["genome"]) if c != "all"]
    cfg.setdefault("seed", SEED)
    if "url" in cfg:
        cfg["url"] = src(cfg["url"])
    if "fastq" in cfg:
        cfg["fastq"] = [src(f) for f in cfg["fastq"]]
    return json.dumps(cfg)


rule fetch_reads_short:
    """A subsample of real paired-end reads (ENA stream, remote BAM slice, or local)."""
    input:
        fasta=lambda wc: genome_fasta(REAL[wc.rd]["genome"]),
        contigs=lambda wc: f"{OUT}/resources/{REAL[wc.rd]['genome']}/contigs.tsv",
    output:
        r1=OUT + "/data/real/{rd}/raw_1.fastq.gz",
        r2=OUT + "/data/real/{rd}/raw_2.fastq.gz",
        provider=OUT + "/data/real/{rd}/raw.provider.parquet",
        provenance=OUT + "/provenance/reads__{rd}.json",
        spec=OUT + "/data/real/{rd}/source.json",
    wildcard_constraints:
        rd=constraint(SHORT_RD),
    params:
        spec=read_spec,
        prefix=OUT + "/data/real/{rd}/raw",
    retries: 2
    shell:
        "printf '%s' {params.spec:q} > {output.spec} && "
        "{ALNFAIL} fetch-reads --spec {output.spec} --fasta {input.fasta} --out-prefix {params.prefix} "
        "--provenance {output.provenance}"


use rule fetch_reads_short as fetch_reads_long with:
    output:
        r1=OUT + "/data/real/{rd}/raw_1.fastq.gz",
        provider=OUT + "/data/real/{rd}/raw.provider.parquet",
        provenance=OUT + "/provenance/reads__{rd}.json",
        spec=OUT + "/data/real/{rd}/source.json",
    wildcard_constraints:
        rd=constraint(LONG_RD),


rule fastp:
    """Short-read cleaning: adapters, poly-G, low-quality and N-rich reads.

    No quality trimming of read ends: the decay of quality along the read is
    part of the error profile we want to learn, not something to cut away."""
    input:
        r1=OUT + "/data/real/{rd}/raw_1.fastq.gz",
        r2=OUT + "/data/real/{rd}/raw_2.fastq.gz",
    output:
        r1=temp(OUT + "/data/real/{rd}/clean_1.fastq.gz"),
        r2=temp(OUT + "/data/real/{rd}/clean_2.fastq.gz"),
        json=OUT + "/qc/{rd}.fastp.json",
        html=OUT + "/qc/{rd}.fastp.html",
    params:
        min_length=lambda wc: PLATFORMS[REAL[wc.rd]["platform"]].get("min_length", 50),
        extra=lambda wc: PLATFORMS[REAL[wc.rd]["platform"]].get("fastp_extra", ""),
    threads: min(int(config.get("threads", 4)), 8)
    log:
        OUT + "/logs/fastp/{rd}.log",
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} "
        "--detect_adapter_for_pe --length_required {params.min_length} {params.extra} "
        "--thread {threads} --json {output.json} --html {output.html} 2> {log}"


rule qc_short:
    """Decision table from fastp, then calibration / held-out split."""
    input:
        r1=OUT + "/data/real/{rd}/clean_1.fastq.gz",
        r2=OUT + "/data/real/{rd}/clean_2.fastq.gz",
        json=OUT + "/qc/{rd}.fastp.json",
    output:
        OUT + "/data/real/{rd}/calib_1.fastq.gz",
        OUT + "/data/real/{rd}/calib_2.fastq.gz",
        OUT + "/data/real/{rd}/heldout_1.fastq.gz",
        OUT + "/data/real/{rd}/heldout_2.fastq.gz",
        table=OUT + "/qc/{rd}.tsv",
    wildcard_constraints:
        rd=constraint(SHORT_RD),
    params:
        prefix=OUT + "/data/real/{rd}",
        min_length=lambda wc: PLATFORMS[REAL[wc.rd]["platform"]].get("min_length", 50),
    shell:
        "{ALNFAIL} qc-short --fastp-json {input.json} --r1 {input.r1} --r2 {input.r2} "
        "--calib-prefix {params.prefix}/calib --heldout-prefix {params.prefix}/heldout "
        "--min-length {params.min_length} --dataset {wildcards.rd} --out {output.table}"


rule qc_long:
    """Long-read filter by length and mean quality, then calibration / held-out split."""
    input:
        OUT + "/data/real/{rd}/raw_1.fastq.gz",
    output:
        calib=OUT + "/data/real/{rd}/calib_1.fastq.gz",
        heldout=OUT + "/data/real/{rd}/heldout_1.fastq.gz",
        table=OUT + "/qc/{rd}.tsv",
    wildcard_constraints:
        rd=constraint(LONG_RD),
    params:
        min_length=lambda wc: PLATFORMS[REAL[wc.rd]["platform"]].get("min_length", 1000),
        min_q=lambda wc: PLATFORMS[REAL[wc.rd]["platform"]].get("min_mean_q", 7),
    shell:
        "{ALNFAIL} qc-long --fastq {input} --calib {output.calib} --heldout {output.heldout} "
        "--min-length {params.min_length} --min-mean-q {params.min_q} --dataset {wildcards.rd} --out {output.table}"


rule provenance:
    """One table listing every source: URL or accession, size, checksum, retrieval date."""
    input:
        [OUT + "/provenance/" + os.path.relpath(p, f"{OUT}/resources").replace("/downloads/", "__") + ".json" for p in DOWNLOADS],
        expand(OUT + "/provenance/reads__{rd}.json", rd=list(REAL)),
    output:
        f"{OUT}/tables/provenance.tsv",
    shell:
        "{ALNFAIL} provenance --records {input} --out {output}"
