# Task 3: calibrate the simulator on real reads.  Tasks 4-5: simulate.
#
#   real calibration reads --minimap2--> raw profile
#   raw profile --simulate pilot--> minimap2 --> pilot profile
#   raw + pilot --> calibrated profile        (injected rates corrected)
#   calibrated profile --> benchmark simulations
#   benchmark simulation --minimap2--> simulated profile --> compared with raw

def profile_kind(rd):
    return PLATFORMS[REAL[rd]["platform"]]["type"]


rule learn_profile_raw:
    """Measure the platform's error profile on real reads (variant sites masked)."""
    input:
        bam=OUT + "/aln/calib-{rd}/minimap2.bam",
        fasta=lambda wc: genome_fasta(REAL[wc.rd]["genome"]),
        vcf=lambda wc: genome_vcf(REAL[wc.rd]["genome"]),
    output:
        npz=OUT + "/profiles/{rd}.raw.npz",
        json=OUT + "/profiles/{rd}.raw.json",
    params:
        kind=lambda wc: profile_kind(wc.rd),
        vcf=lambda wc: vcf_flag(REAL[wc.rd]["genome"]),
        max_reads=lambda wc: 100000 if profile_kind(wc.rd) == "short" else 20000,
    shell:
        "{ALNFAIL} learn-profile --bam {input.bam} --fasta {input.fasta} --kind {params.kind} {params.vcf} "
        "--max-reads {params.max_reads} --source {wildcards.rd} --out {output.npz}"


rule learn_profile_pilot:
    """Measure a pilot simulation exactly as the real reads were measured."""
    input:
        bam=OUT + "/aln/pilot-{rd}/minimap2.bam",
        fasta=lambda wc: genome_fasta(REAL[wc.rd]["genome"]),
        vcf=lambda wc: genome_vcf(REAL[wc.rd]["genome"]),
    output:
        npz=OUT + "/profiles/{rd}.pilot.npz",
        json=OUT + "/profiles/{rd}.pilot.json",
    params:
        kind=lambda wc: profile_kind(wc.rd),
        vcf=lambda wc: vcf_flag(REAL[wc.rd]["genome"]),
    shell:
        "{ALNFAIL} learn-profile --bam {input.bam} --fasta {input.fasta} --kind {params.kind} {params.vcf} "
        "--source pilot-{wildcards.rd} --out {output.npz}"


rule calibrate_profile:
    """Correct the injected error rates so that measured-on-simulation matches measured-on-real."""
    input:
        raw=OUT + "/profiles/{rd}.raw.npz",
        pilot=OUT + "/profiles/{rd}.pilot.npz",
    output:
        npz=OUT + "/profiles/{rd}.npz",
        json=OUT + "/profiles/{rd}.json",
    shell:
        "{ALNFAIL} calibrate --profile {input.raw} --pilot {input.pilot} --out {output.npz}"


def sim_profile(wc):
    info = DS[wc.ds]
    if info["kind"] == "pilot":
        return f"{OUT}/profiles/{info['rd']}.raw.npz"
    return profile_of(info["platform"])


def sim_n(wc):
    info = DS[wc.ds]
    n = int(PLATFORMS[info["platform"]]["n_reads"])
    # the pilot only has to pin down three error rates, so a fraction of the reads is enough
    return max(n // 5, 500) if info["kind"] == "pilot" else n


rule simulate_short:
    """Simulate paired-end reads with known origin."""
    input:
        fasta=lambda wc: genome_fasta(DS[wc.ds]["genome"]),
        profile=sim_profile,
        vcf=lambda wc: genome_vcf(DS[wc.ds]["genome"]),
    output:
        r1=OUT + "/data/sim/{ds}/reads_1.fastq.gz",
        r2=OUT + "/data/sim/{ds}/reads_2.fastq.gz",
        truth=OUT + "/data/sim/{ds}/truth.parquet",
    wildcard_constraints:
        ds=constraint([d for d in SHORT_DS if DS[d]["kind"] in ("sim", "pilot")]),
    params:
        n=sim_n,
        vcf=lambda wc: vcf_flag(DS[wc.ds]["genome"]),
        divergence=lambda wc: DS[wc.ds]["divergence"],
        prefix=OUT + "/data/sim/{ds}/reads",
    benchmark:
        OUT + "/benchmarks/simulate/short/{ds}.tsv"
    shell:
        "{ALNFAIL} simulate --fasta {input.fasta} --profile {input.profile} --n {params.n} --id {wildcards.ds} "
        "--seed {SEED} --divergence {params.divergence} {params.vcf} --out-prefix {params.prefix} --truth {output.truth}"


use rule simulate_short as simulate_long with:
    output:
        r1=OUT + "/data/sim/{ds}/reads_1.fastq.gz",
        truth=OUT + "/data/sim/{ds}/truth.parquet",
    wildcard_constraints:
        ds=constraint([d for d in LONG_DS if DS[d]["kind"] in ("sim", "pilot")]),
    benchmark:
        OUT + "/benchmarks/simulate/long/{ds}.tsv"


def calibration_sim(rd):
    """The benchmark simulation that mirrors a real calibration set."""
    cfg = REAL[rd]
    return f"sim-{cfg['genome']}-{cfg['platform']}-d0"


rule learn_profile_sim:
    """Measure the final simulation for the realism check."""
    input:
        bam=lambda wc: f"{OUT}/aln/{calibration_sim(wc.rd)}/minimap2.bam",
        fasta=lambda wc: genome_fasta(REAL[wc.rd]["genome"]),
        vcf=lambda wc: genome_vcf(REAL[wc.rd]["genome"]),
    output:
        npz=OUT + "/profiles/{rd}.sim.npz",
        json=OUT + "/profiles/{rd}.sim.json",
    params:
        kind=lambda wc: profile_kind(wc.rd),
        vcf=lambda wc: vcf_flag(REAL[wc.rd]["genome"]),
    shell:
        "{ALNFAIL} learn-profile --bam {input.bam} --fasta {input.fasta} --kind {params.kind} {params.vcf} "
        "--source sim-{wildcards.rd} --out {output.npz}"


rule compare_profiles:
    """Task 3 verdict: does the simulated error profile match the real one?"""
    input:
        real=OUT + "/profiles/{rd}.raw.npz",
        sim=OUT + "/profiles/{rd}.sim.npz",
    output:
        OUT + "/calibration/{rd}.tsv",
    params:
        tolerance=EVAL.get("calibration_tolerance", 0.25),
    shell:
        "{ALNFAIL} compare-profiles --real {input.real} --sim {input.sim} --tolerance {params.tolerance} "
        "--dataset {wildcards.rd} --out {output}"
