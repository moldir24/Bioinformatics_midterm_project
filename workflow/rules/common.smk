# Shared definitions: paths, the dataset table, and lookups used by every rule.
#
# The whole experiment is described by the config file. This file turns that
# description into a table of datasets (DS) so that every rule can ask
# "which genome, platform and reads belong to this dataset id?".
import json
import os
import re
from pathlib import Path

ROOT = Path(workflow.basedir).parent  # repository root
OUT = config.get("outdir", "results")  # everything this workflow writes goes below OUT
# run the package from the repository without installing it
ALNFAIL = f"PYTHONPATH={ROOT}${{PYTHONPATH:+:$PYTHONPATH}} python -m alnfail"

GENOMES = config["genomes"]
PLATFORMS = config["platforms"]
REAL = config.get("real_reads") or {}
ALIGNERS = config["aligners"]
OPTIONS = config.get("aligner_options") or {}
SIM = config.get("simulation") or {}
EVAL = config.get("evaluation") or {}
SEED = int(config.get("seed", 26))

def validate_config():
    """Fail early, with a readable message, on a config that cannot work."""
    pattern = re.compile(r"^[A-Za-z0-9_]+$")
    for kind, names in (("genome", GENOMES), ("platform", PLATFORMS), ("real_reads", REAL)):
        for name in names:
            if not pattern.match(name):
                raise ValueError(f"{kind} name {name!r}: use only letters, digits and underscore")
    for name, cfg in PLATFORMS.items():
        if cfg["type"] not in ("short", "long"):
            raise ValueError(f"platform {name}: type must be 'short' (paired-end) or 'long'")
        if cfg["calibrate_on"] not in REAL:
            raise ValueError(f"platform {name}: calibrate_on must name an entry of real_reads")
        if REAL[cfg["calibrate_on"]]["platform"] != name:
            raise ValueError(f"platform {name}: calibrate_on points at reads of another platform")
    for name, cfg in REAL.items():
        for genome in [cfg["genome"]] + list(cfg.get("also_align_to") or []):
            if genome not in GENOMES:
                raise ValueError(f"real_reads {name}: unknown genome {genome!r}")
        if cfg["platform"] not in PLATFORMS:
            raise ValueError(f"real_reads {name}: unknown platform {cfg['platform']!r}")


validate_config()


def src(url):
    """Remote URLs pass through; relative paths are relative to the repository."""
    if "://" in url or os.path.isabs(url):
        return url
    return str(ROOT / url)


# ---------------------------------------------------------------- downloads
DOWNLOADS = {}  # output path -> {url, md5, label}


def register(genome, name, entry):
    base = os.path.basename(entry["url"].split("?")[0])
    path = f"{OUT}/resources/{genome}/downloads/{name}__{base}"
    DOWNLOADS[path] = {"url": src(entry["url"]), "md5": entry.get("md5", ""), "label": f"{genome}:{name}"}
    return path


FASTA_DL, ALIAS_DL, VCF_DL, STRATA_SPEC = {}, {}, {}, {}


def register_genomes():
    for g, cfg in GENOMES.items():
        FASTA_DL[g] = register(g, "fasta", cfg["fasta"])
        ALIAS_DL[g] = register(g, "alias", cfg["alias"]) if cfg.get("alias") else None
        VCF_DL[g] = register(g, "variants", cfg["variants"]["vcf"]) if cfg.get("variants") else None
        STRATA_SPEC[g] = [
            {"name": name, "format": s["format"], "path": register(g, name, s),
             "opts": {k: v for k, v in s.items() if k not in ("format", "url", "md5")}}
            for name, s in (cfg.get("strata") or {}).items()
        ]


register_genomes()


def genome_fasta(g):
    return f"{OUT}/resources/{g}/genome.fa"


def genome_vcf(g):
    """Truth-set variants of this genome's sample, or nothing."""
    return f"{OUT}/resources/{g}/variants.vcf.gz" if VCF_DL[g] else []


def vcf_flag(g):
    return f"--vcf {OUT}/resources/{g}/variants.vcf.gz" if VCF_DL[g] else ""


def contigs_of(g):
    contigs = GENOMES[g].get("contigs", "all")
    return ["all"] if contigs in ("all", None) else list(contigs)


# ----------------------------------------------------------------- datasets
def fmt_div(d):
    return f"{float(d):g}"


def build_datasets():
    """The experiment as a table: dataset id -> what it is.

    sim    simulated reads, one per genome x platform x divergence level
    calib  the calibration half of a real read set (measured to learn the profile)
    pilot  a small simulation used to correct the injected error rates
    real   the held-out half of a real read set, aligned for task 8
    """
    table = {}
    series = SIM.get("divergence_series") or {}
    for g in GENOMES:
        levels = [0.0]
        if g in (series.get("genomes") or []):
            levels += [float(x) for x in series.get("levels", [])]
        for p in PLATFORMS:
            for d in levels:
                table[f"sim-{g}-{p}-d{fmt_div(d)}"] = {"kind": "sim", "genome": g, "platform": p, "divergence": d}
    for rd, cfg in REAL.items():
        p = cfg["platform"]
        table[f"calib-{rd}"] = {"kind": "calib", "genome": cfg["genome"], "platform": p, "rd": rd, "divergence": 0.0}
        if PLATFORMS[p]["calibrate_on"] == rd:
            table[f"pilot-{rd}"] = {"kind": "pilot", "genome": cfg["genome"], "platform": p, "rd": rd, "divergence": 0.0}
        for g in [cfg["genome"]] + list(cfg.get("also_align_to") or []):
            table[f"real-{rd}-on-{g}"] = {"kind": "real", "genome": g, "platform": p, "rd": rd, "divergence": 0.0}
    return table


DS = build_datasets()

SIM_DS = [d for d, v in DS.items() if v["kind"] == "sim"]
REAL_DS = [d for d, v in DS.items() if v["kind"] == "real"]


def ptype(ds):
    return PLATFORMS[DS[ds]["platform"]]["type"]


def is_short(ds):
    return ptype(ds) == "short"


def profile_of(platform):
    """The calibrated error profile a platform's simulations are drawn from."""
    return f"{OUT}/profiles/{PLATFORMS[platform]['calibrate_on']}.npz"


def reads_of(ds):
    """FASTQ file(s) of a dataset: two for short (paired-end), one for long."""
    info = DS[ds]
    mates = (1, 2) if is_short(ds) else (1,)
    if info["kind"] in ("sim", "pilot"):
        return [f"{OUT}/data/sim/{ds}/reads_{m}.fastq.gz" for m in mates]
    half = "calib" if info["kind"] == "calib" else "heldout"
    return [f"{OUT}/data/real/{info['rd']}/{half}_{m}.fastq.gz" for m in mates]


def aligners_of(ds):
    """Benchmark datasets go through every aligner of their read type; the
    calibration helpers only need one alignment to measure error rates."""
    if DS[ds]["kind"] in ("calib", "pilot"):
        return ["minimap2"]
    return list(ALIGNERS[ptype(ds)])


def option(aligner, ds):
    """Aligner options for this dataset's platform, from config.aligner_options."""
    return (OPTIONS.get(aligner) or {}).get(DS[ds]["platform"], "")


def constraint(names):
    return "|".join(re.escape(n) for n in names) or "NONE"


SHORT_DS = [d for d in DS if is_short(d)]
LONG_DS = [d for d in DS if not is_short(d)]
SHORT_RD = [rd for rd, c in REAL.items() if PLATFORMS[c["platform"]]["type"] == "short"]
LONG_RD = [rd for rd, c in REAL.items() if PLATFORMS[c["platform"]]["type"] == "long"]
CALIBRATORS = sorted({p["calibrate_on"] for p in PLATFORMS.values()})

wildcard_constraints:
    genome=constraint(GENOMES),
    rd=constraint(REAL),
    ds=constraint(DS),
    aligner="[a-z0-9]+",


def labels():
    """Display names and lookups handed to the report."""
    return {
        "genomes": {g: c.get("label", g) for g, c in GENOMES.items()},
        "platforms": {p: c.get("label", p) for p, c in PLATFORMS.items()},
        "genome_order": list(GENOMES), "platform_order": list(PLATFORMS),
        "dataset_platform": {d: v["platform"] for d, v in DS.items() if v["kind"] in ("sim", "real")},
    }


def transfer_pairs():
    """(simulated, real) dataset pairs of the same genome and platform."""
    pairs = []
    for real_id in REAL_DS:
        info = DS[real_id]
        sim_id = f"sim-{info['genome']}-{info['platform']}-d0"
        if sim_id in DS:
            pairs.append(f"{sim_id}={real_id}")
    return pairs


def sources_to_check():
    """Every external source in the config, for the week-1 resolvability check."""
    items = [{"label": v["label"], "kind": "download", "url": v["url"]} for v in DOWNLOADS.values()]
    for rd, cfg in REAL.items():
        kind = cfg["source"]
        item = {"label": rd, "kind": kind, "platform": cfg["platform"],
                "paired": PLATFORMS[cfg["platform"]]["type"] == "short"}
        if kind == "ena":
            item.update(accession=cfg["accession"], api=cfg.get("api"))
        elif kind == "bam_slice":
            item.update(url=src(cfg["url"]))
        else:
            item.update(kind="download", url=src(cfg["fastq"][0]))
        items.append(item)
    return items
