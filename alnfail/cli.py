"""Command-line entry points, one per workflow step.

    python -m alnfail <command> --help

The Snakemake workflow calls these commands; each can also be run by hand,
which is the easiest way to debug a single step.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import sys

import pandas as pd


def _kv(items: list[str] | None) -> dict[str, str]:
    """['a=1', 'b=x'] -> {'a': '1', 'b': 'x'}"""
    out = {}
    for item in items or []:
        key, _, value = item.partition("=")
        out[key] = value
    return out


def _makedirs(path: str) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)


def _resolver(fasta: str, alias: str | None):
    import pysam

    from .coords import AliasGraph, ContigResolver

    fa = pysam.FastaFile(fasta)
    names = list(fa.references)
    lengths = {n: fa.get_reference_length(n) for n in names}
    aliases = AliasGraph(alias).aliases_for(names) if alias else None
    return ContigResolver(names, aliases), lengths


def _load_variants(vcf: str | None, fasta: str):
    if not vcf:
        return None
    from .variants import VariantSet

    resolver, _ = _resolver(fasta, None)
    return VariantSet.from_vcf(vcf, resolver)


def _load_strata(directory: str | None):
    from .annot import load_strata

    if not directory:
        return {}
    paths = {os.path.basename(p)[:-4]: p for p in sorted(glob.glob(os.path.join(directory, "*.bed")))}
    return load_strata(paths)


# --------------------------------------------------------------------------- #
# data retrieval and preparation
# --------------------------------------------------------------------------- #
def cmd_fetch_file(a):
    """Download one file and record its provenance."""
    from .fetch import download, write_provenance

    record = download(a.url, a.out, a.md5 or None, a.label)
    write_provenance(a.provenance, record)
    print(f"{a.label}: {record['bytes']:,} bytes, md5 {record['md5']}")


def cmd_check_sources(a):
    """Week-1 check: does every configured source resolve, and how big is it?"""
    from .fetch import bam_contigs, check_url, check_ena_metadata, ena_filereport

    with open(a.sources) as fh:
        sources = json.load(fh)
    rows = []
    for src in sources:
        kind = src["kind"]
        try:
            if kind == "ena":
                row = ena_filereport(src["accession"], src.get("api") or "https://www.ebi.ac.uk/ena/portal/api/filereport")[0]
                problems = check_ena_metadata(row, src["platform"], src["paired"])
                size = sum(int(x) for x in row.get("fastq_bytes", "").split(";") if x)
                detail = f"{row.get('scientific_name')}; {row.get('instrument_model')}; {row.get('library_layout')}"
                rows.append({"label": src["label"], "kind": kind, "source": src["accession"], "ok": not problems,
                             "bytes": size, "detail": detail + ("; PROBLEM: " + "; ".join(problems) if problems else "")})
            elif kind == "bam_slice":
                contigs = bam_contigs(src["url"])
                rows.append({"label": src["label"], "kind": kind, "source": src["url"], "ok": bool(contigs),
                             "bytes": 0, "detail": f"header readable, {len(contigs)} contigs, first: {list(contigs)[:3]}"})
            else:
                result = check_url(src["url"])
                rows.append({"label": src["label"], "kind": kind, "source": src["url"], "ok": result["ok"],
                             "bytes": result["bytes"], "detail": result["detail"]})
        except Exception as exc:  # report every source, do not stop at the first failure
            rows.append({"label": src["label"], "kind": kind, "source": src.get("url") or src.get("accession"),
                         "ok": False, "bytes": 0, "detail": f"{type(exc).__name__}: {exc}"})
    table = pd.DataFrame(rows)
    table["gigabytes"] = (table["bytes"] / 1e9).round(3)
    _makedirs(a.out)
    table.to_csv(a.out, sep="\t", index=False)
    print(table[["label", "kind", "ok", "gigabytes", "detail"]].to_string(index=False))
    bad = table[~table["ok"]]
    if len(bad):
        sys.exit(f"\n{len(bad)} source(s) did not resolve. Fix config before running the pipeline.")
    print(f"\nall {len(table)} sources resolve; total size of whole-file downloads {table['gigabytes'].sum():.2f} GB")


def cmd_prepare_reference(a):
    """Extract the chosen contigs from a (gzipped) FASTA and index the result."""
    import pysam

    from .coords import AliasGraph

    wanted = None if a.contigs in (None, [], ["all"]) else list(a.contigs)
    # the config may name contigs in another archive's spelling than the FASTA uses
    aliases = AliasGraph(a.alias).aliases_for(wanted) if a.alias and wanted else {}
    opener = gzip.open if a.fasta.endswith(".gz") else open
    found, first_names = [], []
    keep = False
    _makedirs(a.out)
    with opener(a.fasta, "rt") as src, open(a.out, "w") as out:
        for line in src:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if len(first_names) < 5:
                    first_names.append(name)
                if wanted is not None:
                    name = name if name in wanted else aliases.get(name, name)
                keep = wanted is None or name in wanted
                if keep:
                    found.append(name)
                    out.write(f">{name}\n")  # drop the description: some aligners keep it in the contig name
            elif keep:
                # UCSC FASTA marks repeats in lower case; we take repeats from the annotation instead
                out.write(line.upper())
    if wanted and set(wanted) - set(found):
        os.remove(a.out)
        sys.exit(f"contigs {sorted(set(wanted) - set(found))} are not in {a.fasta}; "
                 f"its first sequence names are {first_names} (check the naming convention)")
    pysam.faidx(a.out)
    fa = pysam.FastaFile(a.out)
    with open(a.contigs_out, "w") as fh:
        fh.write("contig\tlength\n")
        for name in fa.references:
            fh.write(f"{name}\t{fa.get_reference_length(name)}\n")
    print(f"{len(found)} contigs, {sum(fa.lengths):,} bp -> {a.out}")


def cmd_build_strata(a):
    """Parse annotation sources into normalised stratum BED files."""
    from .annot import build_strata, write_bed

    resolver, lengths = _resolver(a.fasta, a.alias)
    with open(a.spec) as fh:
        spec = json.load(fh)
    os.makedirs(a.outdir, exist_ok=True)
    stats_rows = []
    merged: dict[str, list] = {}
    for source in spec:
        resolver.unresolved.clear()
        sets, stats = build_strata(source["path"], source["format"], source["name"], source.get("opts", {}),
                                   resolver, lengths)
        for set_name, records in sets.items():
            merged.setdefault(set_name, []).extend(records)
            stats_rows.append({**stats, "stratum_set": set_name, "intervals": len(records),
                               "bases": sum(e - s for _, s, e, _ in records)})
    for set_name, records in merged.items():
        write_bed(records, os.path.join(a.outdir, f"{set_name}.bed"))
    pd.DataFrame(stats_rows).to_csv(os.path.join(a.outdir, "strata_stats.tsv"), sep="\t", index=False)
    print(pd.DataFrame(stats_rows).to_string(index=False) if stats_rows else "no annotation sources configured")


def cmd_prepare_variants(a):
    """Rewrite a truth-set VCF with our contig names, keeping only our contigs."""
    resolver, lengths = _resolver(a.fasta, a.alias)
    kept = dropped = 0
    opener = gzip.open if a.vcf.endswith(".gz") else open
    _makedirs(a.out)
    with opener(a.vcf, "rt") as src, gzip.open(a.out, "wt", compresslevel=4) as out:
        for line in src:
            if line.startswith("##contig"):
                continue
            if line.startswith("#"):
                out.write(line)
                continue
            contig, rest = line.split("\t", 1)
            name = resolver.resolve(contig)
            if name is None:
                dropped += 1
                continue
            pos = int(rest.split("\t", 1)[0])
            if pos > lengths[name]:
                sys.exit(f"variant at {contig}:{pos} lies beyond the end of {name} ({lengths[name]} bp): "
                         "the VCF is for a different assembly than the FASTA")
            out.write(f"{name}\t{rest}")
            kept += 1
    with open(a.stats, "w") as fh:
        json.dump({"variants_kept": kept, "variants_on_other_contigs": dropped,
                   "contig_names": resolver.report()}, fh, indent=2)
    print(f"{kept:,} variants kept, {dropped:,} on contigs we do not use")


def cmd_fetch_reads(a):
    """Retrieve a subsample of real reads (ENA stream, remote BAM slice, or local files)."""
    from . import fetch

    with open(a.spec) as fh:
        spec = json.load(fh)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    paired = bool(spec["paired"])
    source = spec["source"]
    if source == "ena":
        record = fetch.fetch_ena(spec["accession"], spec["platform"], paired, int(spec["n_reads"]),
                                 int(spec.get("skip", 0)), a.out_prefix, spec.get("api") or fetch.ENA_API)
    elif source == "bam_slice":
        record = fetch.fetch_bam_slice(spec["url"], a.fasta, spec["contigs"], spec["platform"], paired,
                                       int(spec["n_reads"]), int(spec["n_windows"]), int(spec["window_bp"]),
                                       int(spec.get("seed", 26)), a.out_prefix, spec.get("md5"))
    elif source == "local":
        record = fetch.fetch_local(spec["fastq"], int(spec["n_reads"]), a.out_prefix)
    else:
        sys.exit(f"unknown read source {source!r} (use ena, bam_slice or local)")
    record["label"] = spec.get("label", record["label"])
    fetch.write_provenance(a.provenance, record)
    if record.get("metadata_problems"):
        print("WARNING, archive metadata disagrees with the config:", "; ".join(record["metadata_problems"]))
    if not os.path.exists(a.out_prefix + ".provider.parquet"):  # only BAM slices carry the provider's alignment
        from .fetch import pa, pq
        pq.write_table(pa.table({"qname": pa.array([], pa.string()), "mate": pa.array([], pa.int8()),
                                 "contig": pa.array([], pa.string()), "start": pa.array([], pa.int64()),
                                 "end": pa.array([], pa.int64()), "mapq": pa.array([], pa.int16())}),
                       a.out_prefix + ".provider.parquet")
    print(json.dumps({k: record[k] for k in ("label", "kind", "bytes")}, indent=None))


def cmd_provenance(a):
    """Merge provenance records into one table."""
    from .fetch import merge_provenance

    _makedirs(a.out)
    merge_provenance(a.records, a.out)


# --------------------------------------------------------------------------- #
# QC
# --------------------------------------------------------------------------- #
def cmd_qc_long(a):
    """Filter long reads by length and mean quality; split into calibration and held-out halves."""
    from .qc import filter_long

    rows = filter_long(a.fastq, a.calib, a.heldout, a.min_length, a.min_mean_q, a.dataset)
    _makedirs(a.out)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)


def cmd_qc_short(a):
    """Summarise fastp's report as a decision table; split pairs into calibration and held-out halves."""
    from .qc import fastp_decisions, split_pairs

    rows = fastp_decisions(a.fastp_json, a.dataset, a.min_length)
    n = split_pairs(a.r1, a.r2, a.calib_prefix, a.heldout_prefix)
    rows.append({"dataset": a.dataset, "step": "split", "metric": "fragments_after_qc", "value": n, "threshold": "",
                 "why": "alternate fragments go to the calibration half and the held-out half"})
    _makedirs(a.out)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)


# --------------------------------------------------------------------------- #
# calibration and simulation
# --------------------------------------------------------------------------- #
def cmd_learn_profile(a):
    """Learn a platform error profile from aligned reads."""
    from .profile import learn_long, learn_short, save_profile

    variants = _load_variants(a.vcf, a.fasta)
    learn = learn_short if a.kind == "short" else learn_long
    prof = learn(a.bam, a.fasta, variants, max_reads=a.max_reads, min_mapq=a.min_mapq, source=a.source)
    _makedirs(a.out)
    save_profile(a.out, prof)
    print(json.dumps(prof["meta"], indent=2))


def cmd_simulate(a):
    """Simulate reads with known origin from a profile."""
    from .profile import load_profile
    from .simulate import simulate_long, simulate_short

    prof = load_profile(a.profile)
    variants = _load_variants(a.vcf, a.fasta)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    if prof["meta"]["kind"] == "short":
        r2 = f"{a.out_prefix}_2.fastq.gz" if prof["meta"]["paired"] else None
        info = simulate_short(a.fasta, prof, a.n, f"{a.out_prefix}_1.fastq.gz", r2, a.truth, a.seed, a.id,
                              variants=variants, divergence=a.divergence)
    else:
        info = simulate_long(a.fasta, prof, a.n, f"{a.out_prefix}_1.fastq.gz", a.truth, a.seed, a.id,
                             variants=variants, divergence=a.divergence)
    print(json.dumps(info))


def cmd_calibrate(a):
    """Close the loop: correct the injected error rates using a pilot simulation."""
    from .profile import load_profile, refine_profile, save_profile

    prof = refine_profile(load_profile(a.profile), load_profile(a.pilot))
    _makedirs(a.out)
    save_profile(a.out, prof)
    print(json.dumps(prof["meta"]["inject"], indent=2))


def cmd_compare_profiles(a):
    """Compare the error profile of real reads with that of simulated reads."""
    from .profile import compare_profiles, load_profile

    rows = compare_profiles(load_profile(a.real), load_profile(a.sim), a.tolerance)
    table = pd.DataFrame(rows)
    table.insert(0, "dataset", a.dataset)
    _makedirs(a.out)
    table.to_csv(a.out, sep="\t", index=False)
    print(table.to_string(index=False))


# --------------------------------------------------------------------------- #
# evaluation
# --------------------------------------------------------------------------- #
def cmd_bam_table(a):
    """Reduce a BAM to a Parquet table of primary alignments."""
    from .samtable import bam_to_table

    _makedirs(a.out)
    print(json.dumps(bam_to_table(a.bam, a.out)))


def cmd_annotate_truth(a):
    """Attach stratum labels to the true intervals of simulated reads."""
    from .summarise import annotate_intervals

    truth = pd.read_parquet(a.truth)
    variants = _load_variants(a.vcf, a.fasta)
    out = annotate_intervals(truth, _load_strata(a.strata), a.fasta,
                             variants.positions if variants else None, long_reads=a.long, divergence=a.divergence)
    _makedirs(a.out)
    out.to_parquet(a.out, compression="zstd", index=False)


def cmd_evaluate(a):
    """Score each simulated read: correct, misplaced or unmapped."""
    from .evaluate import evaluate

    result = evaluate(pd.read_parquet(a.truth, columns=["qname", "mate", "contig", "start", "end", "strand"]),
                      pd.read_parquet(a.table), a.min_overlap)
    _makedirs(a.out)
    result.to_parquet(a.out, compression="zstd", index=False)
    print(result["outcome"].value_counts().to_string())


def _strata_names(directory: str | None) -> list[str]:
    if not directory:
        return []
    return [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(directory, "*.bed")))]


def _truth_free(a, tables: dict, meta: dict, truth=None):
    """Consensus placement, its strata, and the truth-free observables.

    Shared by simulated and real datasets so that task 8 compares numbers
    produced by one identical procedure.
    """
    from .consensus import consensus, consensus_summary, observables
    from .summarise import annotate_intervals, stratum_columns

    cons = consensus(tables, a.min_overlap)
    variants = _load_variants(a.vcf, a.fasta)
    placed = cons[cons["has_consensus"]]
    annotated = annotate_intervals(placed, _load_strata(a.strata), a.fasta,
                                   variants.positions if variants else None, long_reads=a.long)
    cols = stratum_columns(annotated, _strata_names(a.strata))
    full = pd.concat([annotated, cons[~cons["has_consensus"]]], ignore_index=True)
    consensus_summary(full, list(tables), meta, cols, truth=truth, min_overlap=a.min_overlap).to_csv(
        os.path.join(a.outdir, "consensus.tsv"), sep="\t", index=False)
    observables(tables, meta).to_csv(os.path.join(a.outdir, "observables.tsv"), sep="\t", index=False)


def cmd_summarise(a):
    """Tasks 4-7 on one simulated dataset, plus its truth-free measurements."""
    from .summarise import stratum_columns, summarise_dataset

    meta = _kv(a.meta)
    truth = pd.read_parquet(a.truth)
    cols = stratum_columns(truth, _strata_names(a.strata))
    evals = {k: pd.read_parquet(v) for k, v in _kv(a.evals).items()}
    tables = {k: pd.read_parquet(v) for k, v in _kv(a.tables).items()}
    os.makedirs(a.outdir, exist_ok=True)
    for name, frame in summarise_dataset(truth, evals, meta, cols).items():
        frame.to_csv(os.path.join(a.outdir, f"{name}.tsv"), sep="\t", index=False)
    _truth_free(a, tables, meta, truth=truth)


def cmd_real(a):
    """Task 8 on one real dataset: consensus, strata, observables."""
    from .consensus import provider_agreement

    meta = _kv(a.meta)
    tables = {k: pd.read_parquet(v) for k, v in _kv(a.tables).items()}
    os.makedirs(a.outdir, exist_ok=True)
    _truth_free(a, tables, meta)
    provider = pd.read_parquet(a.provider) if a.provider and os.path.exists(a.provider) else pd.DataFrame()
    agreement = provider_agreement(tables, provider, meta, a.min_overlap) if len(provider) else pd.DataFrame(
        columns=["dataset", "aligner", "n_compared", "n_mapped", "n_agree_provider", "n_confident",
                 "n_confident_disagree_provider"])
    agreement.to_csv(os.path.join(a.outdir, "provider.tsv"), sep="\t", index=False)


def cmd_collect(a):
    """Merge per-dataset tables and derive the report tables."""
    from .consensus import transfer_observables, transfer_strata
    from .summarise import add_rates, calibration_table, headline

    os.makedirs(a.outdir, exist_ok=True)

    def gather(name):
        frames = []
        for directory in a.datasets:
            path = os.path.join(directory, f"{name}.tsv")
            if os.path.exists(path) and os.path.getsize(path) > 1:
                frame = pd.read_csv(path, sep="\t", dtype={"stratum": str, "key": str, "divergence": str})
                if len(frame):
                    frames.append(frame)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    accuracy = add_rates(gather("accuracy"))
    mapq = gather("mapq")
    strata = add_rates(gather("strata"))
    accuracy.to_csv(os.path.join(a.outdir, "accuracy.tsv"), sep="\t", index=False)
    strata.to_csv(os.path.join(a.outdir, "strata.tsv"), sep="\t", index=False)
    gather("misplaced_kinds").to_csv(os.path.join(a.outdir, "misplaced_kinds.tsv"), sep="\t", index=False)

    calibration_table(mapq, ["dataset", "genome", "platform", "divergence", "aligner"]).to_csv(
        os.path.join(a.outdir, "mapq_calibration.tsv"), sep="\t", index=False)
    # the headline pools the genomes at zero divergence: reads from the reference's own species
    base = mapq[mapq["divergence"].astype(float) == 0]
    pooled = calibration_table(base, ["platform", "aligner"])
    pooled.to_csv(os.path.join(a.outdir, "mapq_calibration_pooled.tsv"), sep="\t", index=False)
    head = headline(accuracy[accuracy["divergence"].astype(float) == 0], pooled, ["platform", "aligner"])
    head.to_csv(os.path.join(a.outdir, "headline.tsv"), sep="\t", index=False)

    cons, obs = gather("consensus"), gather("observables")
    cons.to_csv(os.path.join(a.outdir, "consensus.tsv"), sep="\t", index=False)
    obs.to_csv(os.path.join(a.outdir, "observables.tsv"), sep="\t", index=False)
    gather("provider").to_csv(os.path.join(a.outdir, "provider_agreement.tsv"), sep="\t", index=False)

    pairs = [tuple(p.split("=", 1)) for p in a.pairs or []]
    t_obs = transfer_observables(obs, pairs) if pairs and len(obs) else pd.DataFrame()
    joined, t_strata = transfer_strata(strata, cons, pairs, a.min_reads) if pairs and len(cons) else (pd.DataFrame(), pd.DataFrame())
    t_obs.to_csv(os.path.join(a.outdir, "transfer_observables.tsv"), sep="\t", index=False)
    joined.to_csv(os.path.join(a.outdir, "transfer_strata_points.tsv"), sep="\t", index=False)
    t_strata.to_csv(os.path.join(a.outdir, "transfer_strata.tsv"), sep="\t", index=False)

    for name, files in (("calibration_check", a.calibration), ("qc_decisions", a.qc), ("strata_sources", a.strata_stats)):
        frames = [pd.read_csv(f, sep="\t") for f in files or [] if os.path.getsize(f) > 1]
        (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()).to_csv(
            os.path.join(a.outdir, f"{name}.tsv"), sep="\t", index=False)
    print(head.to_string(index=False))


def cmd_benchmarks(a):
    """Collect Snakemake benchmark files into one resource table."""
    rows = []
    for path in a.files:
        frame = pd.read_csv(path, sep="\t")
        parts = os.path.normpath(path).split(os.sep)
        frame.insert(0, "target", os.path.splitext(parts[-1])[0])
        frame.insert(0, "group", parts[-2])
        frame.insert(0, "step", parts[-3])
        rows.append(frame)
    table = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    _makedirs(a.out)
    table.to_csv(a.out, sep="\t", index=False)


def cmd_report(a):
    """Figures and the results digest."""
    from .report import build_report

    build_report(a.tables, a.outdir, a.meta, a.benchmarks)


def cmd_html(a):
    """One self-contained HTML page from summary.md and the figures."""
    from .html import build_html

    print(f"written: {build_html(a.report, a.out)}")


def cmd_mapq_theory(a):
    """Tables behind the task-1 explanation of MAPQ."""
    from .mapq_theory import reference_table, two_copy_experiment

    os.makedirs(a.outdir, exist_ok=True)
    reference_table().to_csv(os.path.join(a.outdir, "mapq_reference.tsv"), sep="\t", index=False)
    # the exact posterior, checked by simulation: the yardstick for task 6
    two_copy_experiment(error_rate=a.error_rate, reads_per_setting=a.reads).to_csv(
        os.path.join(a.outdir, "mapq_ideal_calibration.tsv"), sep="\t", index=False)


def cmd_make_test_data(a):
    """Rebuild the toy dataset in test_data/."""
    from .testdata import make_test_data

    make_test_data(a.outdir, a.seed)


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="alnfail", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def add(name, func, *arguments):
        p = sub.add_parser(name, help=(func.__doc__ or "").strip().split("\n")[0])
        for flags, kwargs in arguments:
            p.add_argument(*flags, **kwargs)
        p.set_defaults(func=func)

    def arg(*flags, **kwargs):
        return flags, kwargs

    add("fetch-file", cmd_fetch_file, arg("--url", required=True), arg("--out", required=True), arg("--md5", default=""),
        arg("--label", default=""), arg("--provenance", required=True))
    add("check-sources", cmd_check_sources, arg("--sources", required=True), arg("--out", required=True))
    add("prepare-reference", cmd_prepare_reference, arg("--fasta", required=True), arg("--contigs", nargs="*"),
        arg("--alias"),
        arg("--out", required=True), arg("--contigs-out", required=True))
    add("build-strata", cmd_build_strata, arg("--fasta", required=True), arg("--spec", required=True),
        arg("--alias"), arg("--outdir", required=True))
    add("prepare-variants", cmd_prepare_variants, arg("--vcf", required=True), arg("--fasta", required=True),
        arg("--alias"), arg("--out", required=True), arg("--stats", required=True))
    add("fetch-reads", cmd_fetch_reads, arg("--spec", required=True), arg("--fasta"), arg("--out-prefix", required=True),
        arg("--provenance", required=True))
    add("provenance", cmd_provenance, arg("--records", nargs="+", required=True), arg("--out", required=True))
    add("qc-long", cmd_qc_long, arg("--fastq", required=True), arg("--calib", required=True),
        arg("--heldout", required=True), arg("--min-length", type=int, default=1000),
        arg("--min-mean-q", type=float, default=7.0), arg("--dataset", default=""), arg("--out", required=True))
    add("qc-short", cmd_qc_short, arg("--fastp-json", required=True), arg("--r1", required=True), arg("--r2"),
        arg("--calib-prefix", required=True), arg("--heldout-prefix", required=True),
        arg("--min-length", type=int, default=50), arg("--dataset", default=""), arg("--out", required=True))
    add("learn-profile", cmd_learn_profile, arg("--bam", required=True), arg("--fasta", required=True),
        arg("--kind", choices=["short", "long"], required=True), arg("--vcf"), arg("--max-reads", type=int, default=100_000),
        arg("--min-mapq", type=int, default=20), arg("--source", default=""), arg("--out", required=True))
    add("simulate", cmd_simulate, arg("--fasta", required=True), arg("--profile", required=True),
        arg("--n", type=int, required=True), arg("--id", required=True), arg("--seed", type=int, default=26),
        arg("--divergence", type=float, default=0.0), arg("--vcf"), arg("--out-prefix", required=True),
        arg("--truth", required=True))
    add("calibrate", cmd_calibrate, arg("--profile", required=True), arg("--pilot", required=True),
        arg("--out", required=True))
    add("compare-profiles", cmd_compare_profiles, arg("--real", required=True), arg("--sim", required=True),
        arg("--tolerance", type=float, default=0.25), arg("--dataset", default=""), arg("--out", required=True))
    add("bam-table", cmd_bam_table, arg("--bam", required=True), arg("--out", required=True))
    add("annotate-truth", cmd_annotate_truth, arg("--truth", required=True), arg("--strata"), arg("--fasta", required=True),
        arg("--vcf"), arg("--long", action="store_true"), arg("--divergence", type=float), arg("--out", required=True))
    add("evaluate", cmd_evaluate, arg("--truth", required=True), arg("--table", required=True),
        arg("--min-overlap", type=float, default=0.1), arg("--out", required=True))
    add("summarise", cmd_summarise, arg("--truth", required=True), arg("--evals", nargs="+", required=True),
        arg("--tables", nargs="+", required=True), arg("--strata"), arg("--fasta", required=True), arg("--vcf"),
        arg("--long", action="store_true"), arg("--meta", nargs="*"),
        arg("--min-overlap", type=float, default=0.1), arg("--outdir", required=True))
    add("real", cmd_real, arg("--tables", nargs="+", required=True), arg("--strata"), arg("--fasta", required=True),
        arg("--vcf"), arg("--provider"), arg("--long", action="store_true"), arg("--meta", nargs="*"),
        arg("--min-overlap", type=float, default=0.1), arg("--outdir", required=True))
    add("collect", cmd_collect, arg("--datasets", nargs="+", required=True), arg("--pairs", nargs="*"),
        arg("--calibration", nargs="*"), arg("--qc", nargs="*"), arg("--strata-stats", nargs="*"),
        arg("--min-reads", type=int, default=100), arg("--outdir", required=True))
    add("benchmarks", cmd_benchmarks, arg("--files", nargs="+", required=True), arg("--out", required=True))
    add("report", cmd_report, arg("--tables", required=True), arg("--outdir", required=True), arg("--meta"),
        arg("--benchmarks"))
    add("html", cmd_html, arg("--report", required=True), arg("--out"))
    add("mapq-theory", cmd_mapq_theory, arg("--outdir", required=True), arg("--error-rate", type=float, default=0.05),
        arg("--reads", type=int, default=1_000_000))
    add("make-test-data", cmd_make_test_data, arg("--outdir", required=True), arg("--seed", type=int, default=26))
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
