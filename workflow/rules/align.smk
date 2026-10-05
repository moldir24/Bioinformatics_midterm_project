# Task 6: every aligner on identical reads.
#
# Each aligner has an index rule (per genome) and a map rule (per dataset).
# All map rules write an unsorted BAM in input order and are timed by
# Snakemake's `benchmark`, so runtime and peak memory are measured, not guessed.
# Options come from config.aligner_options and are justified in docs/aligners.md.

THREADS = int(config.get("threads", 4))
IDX = OUT + "/resources/{genome}/index"


def idx(aligner):
    return lambda wc: f"{OUT}/resources/{DS[wc.ds]['genome']}/index/{aligner}"


def fasta_of(wc):
    return genome_fasta(DS[wc.ds]["genome"])


def reads(wc):
    return reads_of(wc.ds)


# ----------------------------------------------------------------- BWA-MEM2
rule index_bwamem2:
    input:
        OUT + "/resources/{genome}/genome.fa",
    output:
        multiext(IDX + "/bwamem2/genome", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    params:
        prefix=IDX + "/bwamem2/genome",
    log:
        OUT + "/logs/index/bwamem2/{genome}.log",
    benchmark:
        OUT + "/benchmarks/index/bwamem2/{genome}.tsv"
    shell:
        "bwa-mem2 index -p {params.prefix} {input} > {log} 2>&1"


rule map_bwamem2:
    input:
        reads=reads,
        index=lambda wc: multiext(idx("bwamem2")(wc) + "/genome", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        OUT + "/aln/{ds}/bwamem2.bam",
    params:
        prefix=lambda wc: idx("bwamem2")(wc) + "/genome",
        opts=lambda wc: option("bwamem2", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/bwamem2/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/bwamem2/{ds}.tsv"
    shell:
        "bwa-mem2 mem -t {threads} {params.opts} {params.prefix} {input.reads} 2> {log} | samtools view -b -o {output} -"


# ---------------------------------------------------------------------- BWA
rule index_bwa:
    input:
        OUT + "/resources/{genome}/genome.fa",
    output:
        multiext(IDX + "/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        prefix=IDX + "/bwa/genome",
    log:
        OUT + "/logs/index/bwa/{genome}.log",
    benchmark:
        OUT + "/benchmarks/index/bwa/{genome}.tsv"
    shell:
        "bwa index -p {params.prefix} {input} > {log} 2>&1"


rule map_bwa:
    input:
        reads=reads,
        index=lambda wc: multiext(idx("bwa")(wc) + "/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        OUT + "/aln/{ds}/bwa.bam",
    params:
        prefix=lambda wc: idx("bwa")(wc) + "/genome",
        opts=lambda wc: option("bwa", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/bwa/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/bwa/{ds}.tsv"
    shell:
        "bwa mem -t {threads} {params.opts} {params.prefix} {input.reads} 2> {log} | samtools view -b -o {output} -"


# ------------------------------------------------------------------ Bowtie 2
rule index_bowtie2:
    input:
        OUT + "/resources/{genome}/genome.fa",
    output:
        multiext(IDX + "/bowtie2/genome", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    params:
        prefix=IDX + "/bowtie2/genome",
    threads: THREADS
    log:
        OUT + "/logs/index/bowtie2/{genome}.log",
    benchmark:
        OUT + "/benchmarks/index/bowtie2/{genome}.tsv"
    shell:
        "bowtie2-build --threads {threads} {input} {params.prefix} > {log} 2>&1"


rule map_bowtie2:
    """-X (maximum fragment length) is taken from the measured fragment-length
    distribution: Bowtie 2's default of 500 bp would call every longer real
    fragment discordant."""
    input:
        reads=reads,
        index=lambda wc: multiext(idx("bowtie2")(wc) + "/genome", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        profile=lambda wc: f"{OUT}/profiles/{PLATFORMS[DS[wc.ds]['platform']]['calibrate_on']}.raw.json",
    output:
        OUT + "/aln/{ds}/bowtie2.bam",
    params:
        prefix=lambda wc: idx("bowtie2")(wc) + "/genome",
        opts=lambda wc: option("bowtie2", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/bowtie2/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/bowtie2/{ds}.tsv"
    shell:
        "MAXINS=$(python -c \"import json; print(max(500, int(json.load(open('{input.profile}'))['fragment_max'])))\") && "
        "bowtie2 -p {threads} -X $MAXINS {params.opts} -x {params.prefix} -1 {input.reads[0]} -2 {input.reads[1]} 2> {log} "
        "| samtools view -b -o {output} -"


# ---------------------------------------------------------------- strobealign
rule map_strobealign:
    """strobealign builds its index on the fly (seconds to a minute), so its
    map time includes indexing; this is how the tool is meant to be run."""
    input:
        reads=reads,
        fasta=fasta_of,
    output:
        OUT + "/aln/{ds}/strobealign.bam",
    params:
        opts=lambda wc: option("strobealign", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/strobealign/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/strobealign/{ds}.tsv"
    shell:
        "strobealign -t {threads} {params.opts} {input.fasta} {input.reads} 2> {log} | samtools view -b -o {output} -"


# ------------------------------------------------------------------ minimap2
rule map_minimap2:
    """One tool, three presets: -ax sr for short reads, map-ont and map-hifi
    for long reads. The minimizer index is built on the fly because it depends
    on the preset."""
    input:
        reads=reads,
        fasta=fasta_of,
    output:
        OUT + "/aln/{ds}/minimap2.bam",
    params:
        opts=lambda wc: option("minimap2", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/minimap2/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/minimap2/{ds}.tsv"
    shell:
        "minimap2 -t {threads} {params.opts} {input.fasta} {input.reads} 2> {log} | samtools view -b -o {output} -"


# ----------------------------------------------------------------- Winnowmap
rule index_winnowmap:
    """Winnowmap down-weights the most frequent k-mers, which first have to be counted."""
    input:
        OUT + "/resources/{genome}/genome.fa",
    output:
        kmers=IDX + "/winnowmap/repetitive_k15.txt",
        db=temp(directory(IDX + "/winnowmap/merylDB")),
    log:
        OUT + "/logs/index/winnowmap/{genome}.log",
    benchmark:
        OUT + "/benchmarks/index/winnowmap/{genome}.tsv"
    shell:
        "meryl count k=15 output {output.db} {input} > {log} 2>&1 && "
        "meryl print greater-than distinct=0.9998 {output.db} > {output.kmers} 2>> {log}"


rule map_winnowmap:
    input:
        reads=reads,
        fasta=fasta_of,
        kmers=lambda wc: idx("winnowmap")(wc) + "/repetitive_k15.txt",
    output:
        OUT + "/aln/{ds}/winnowmap.bam",
    params:
        opts=lambda wc: option("winnowmap", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/winnowmap/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/winnowmap/{ds}.tsv"
    shell:
        "winnowmap -W {input.kmers} -t {threads} {params.opts} {input.fasta} {input.reads} 2> {log} "
        "| samtools view -b -o {output} -"


# --------------------------------------------------------------------- NGMLR
rule index_ngmlr:
    """NGMLR writes its index next to the reference on first use. Building it
    once here keeps parallel map jobs from racing to create the same files."""
    input:
        OUT + "/resources/{genome}/genome.fa",
    output:
        fasta=IDX + "/ngmlr/genome.fa",
        enc=IDX + "/ngmlr/genome.fa-enc.2.ngm",
        table=IDX + "/ngmlr/genome.fa-ht-13-2.2.ngm",
    log:
        OUT + "/logs/index/ngmlr/{genome}.log",
    benchmark:
        OUT + "/benchmarks/index/ngmlr/{genome}.tsv"
    shell:
        "ln -sf $(realpath {input}) {output.fasta} && "
        "printf '@seed\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' > {output.fasta}.seed.fq && "
        "ngmlr -t 1 -r {output.fasta} -q {output.fasta}.seed.fq -o /dev/null > {log} 2>&1 && rm -f {output.fasta}.seed.fq"


rule map_ngmlr:
    input:
        reads=reads,
        fasta=lambda wc: idx("ngmlr")(wc) + "/genome.fa",
        index=lambda wc: [idx("ngmlr")(wc) + "/genome.fa-enc.2.ngm", idx("ngmlr")(wc) + "/genome.fa-ht-13-2.2.ngm"],
    output:
        bam=OUT + "/aln/{ds}/ngmlr.bam",
        sam=temp(OUT + "/aln/{ds}/ngmlr.sam"),
    params:
        opts=lambda wc: option("ngmlr", wc.ds),
    threads: THREADS
    log:
        OUT + "/logs/map/ngmlr/{ds}.log",
    benchmark:
        OUT + "/benchmarks/map/ngmlr/{ds}.tsv"
    shell:
        "ngmlr -t {threads} {params.opts} -r {input.fasta} -q {input.reads} -o {output.sam} > {log} 2>&1 && "
        "samtools view -b -o {output.bam} {output.sam}"
