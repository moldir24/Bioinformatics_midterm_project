"""Scripted, logged data retrieval (task 2).

Nothing in this project is downloaded by hand. Every file that enters the
pipeline comes through one of the functions below, and every one of them
writes a provenance record: where the bytes came from, when, how many, and
their checksum. The records are merged into results/provenance.tsv, which is
the accession list the report must contain.

Three kinds of source:

  download   a whole file by URL (references, annotations, truth sets)
  ena        the first N reads of a run, found through the ENA portal API by
             run accession and streamed, so that a 100 GB run costs a few
             hundred megabytes of traffic
  bam_slice  reads overlapping chosen windows of a remote, indexed BAM. The
             index lets samtools fetch only those byte ranges ("index rather
             than scan"), which is how a 300x human genome becomes usable
             from a laptop.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import tempfile
import urllib.parse
import urllib.request

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .coords import ContigResolver, region_string
from .seqio import CODE, read_fastq, to_array

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = ("run_accession,study_accession,sample_accession,scientific_name,instrument_platform,"
              "instrument_model,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes")

# what the ENA metadata should say for each of our platform labels
EXPECTED_PLATFORM = {"illumina": "ILLUMINA", "ont": "OXFORD_NANOPORE", "hifi": "PACBIO_SMRT"}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5sum(path: str, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def is_local(url: str) -> bool:
    return "://" not in url or url.startswith("file://")


def local_path(url: str) -> str:
    return url[7:] if url.startswith("file://") else url


def write_provenance(path: str, record: dict) -> None:
    with open(path, "w") as fh:
        json.dump(record, fh, indent=2)


def download(url: str, dest: str, md5: str | None = None, label: str = "") -> dict:
    """Fetch one file and return its provenance record.

    curl resumes partial downloads (-C -) and retries transient failures, so a
    dropped connection in the middle of a 1 GB reference does not restart it.
    When an expected MD5 is given, a mismatch is an error: a truncated or
    silently replaced file must stop the pipeline, not flow into the results.
    """
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    if is_local(url):
        shutil.copyfile(local_path(url), dest)
    else:
        part = dest + ".part"
        cmd = ["curl", "--fail", "--location", "--silent", "--show-error", "--retry", "5", "--retry-delay", "10",
               "--continue-at", "-", "--output", part, url]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"download failed ({result.returncode}) for {url}\n{result.stderr.strip()}")
        os.replace(part, dest)
    observed = md5sum(dest)
    if md5 and observed != md5:
        raise RuntimeError(f"checksum mismatch for {url}: expected {md5}, got {observed}")
    return {"label": label, "kind": "download", "url": url, "path": dest, "bytes": os.path.getsize(dest),
            "md5": observed, "md5_expected": md5 or "", "md5_verified": bool(md5), "retrieved_utc": now()}


def check_url(url: str, timeout: int = 60) -> dict:
    """Does this source resolve, and how big is it? (week-1 check)"""
    if is_local(url):
        path = local_path(url)
        ok = os.path.exists(path)
        return {"url": url, "ok": ok, "bytes": os.path.getsize(path) if ok else 0, "detail": "local file"}
    result = subprocess.run(["curl", "--head", "--location", "--silent", "--show-error", "--fail",
                             "--max-time", str(timeout), url], capture_output=True, text=True)
    size = 0
    for line in result.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            size = int(line.split(":")[1].strip() or 0)
    return {"url": url, "ok": result.returncode == 0, "bytes": size,
            "detail": result.stderr.strip() or "resolves"}


# --------------------------------------------------------------------------- #
# ENA
# --------------------------------------------------------------------------- #
def ena_filereport(accession: str, api: str = ENA_API) -> list[dict]:
    """Ask the ENA portal API what files a run accession has.

    ENA, NCBI SRA and DDBJ mirror the same runs, but ENA serves ready-made
    FASTQ over HTTPS together with per-file MD5 and the submitter's metadata,
    which is what lets us both stream the reads and sanity-check them.
    """
    query = urllib.parse.urlencode({"accession": accession, "result": "read_run", "fields": ENA_FIELDS,
                                    "format": "tsv"})
    url = f"{api}?{query}"
    if is_local(api):  # tests point this at a saved report
        text = open(local_path(api)).read()
    else:
        with urllib.request.urlopen(url, timeout=120) as response:
            text = response.read().decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA knows no run called {accession} ({url})")
    return rows


def check_ena_metadata(row: dict, platform: str, paired: bool) -> list[str]:
    """Compare what the archive says about a run with what we assume.

    Returns human-readable problems (empty when consistent). A run labelled
    with another instrument, or single-end where we expect pairs, would
    silently calibrate the wrong error model.
    """
    problems = []
    expected = EXPECTED_PLATFORM.get(platform)
    if expected and row.get("instrument_platform") != expected:
        problems.append(f"platform is {row.get('instrument_platform')}, expected {expected}")
    layout = row.get("library_layout", "")
    if paired and layout != "PAIRED":
        problems.append(f"library layout is {layout}, expected PAIRED")
    n_files = len([u for u in row.get("fastq_ftp", "").split(";") if u])
    if paired and n_files < 2:
        problems.append(f"only {n_files} FASTQ file(s) for a paired-end run")
    if n_files == 0:
        problems.append("ENA lists no FASTQ files for this run")
    return problems


def _stream_head(url: str, dest: str, n_reads: int, skip: int) -> None:
    """Write reads [skip, skip + n_reads) of a remote FASTQ.gz to `dest`.

    The download is cut off as soon as enough reads have arrived. `skip`
    exists because the first reads of a file come from the edge of the flow
    cell and are of below-average quality.
    """
    first = 4 * skip + 1
    if is_local(url):
        source = f"gzip -dc '{local_path(url)}'"
    else:
        source = f"curl --fail --location --silent --retry 5 '{url}' | gzip -dc"
    # head closing the pipe makes the upstream commands exit with SIGPIPE; that is the intended early stop
    cmd = f"set +o pipefail; {source} 2>/dev/null | tail -n +{first} | head -n {4 * n_reads} | gzip -c > '{dest}'"
    subprocess.run(["bash", "-c", cmd], check=True)
    if os.path.getsize(dest) < 100:
        raise RuntimeError(f"no reads could be streamed from {url}")


def fetch_ena(accession: str, platform: str, paired: bool, n_reads: int, skip: int, out_prefix: str,
              api: str = ENA_API) -> dict:
    """Stream a subsample of one ENA run. Writes <prefix>_1.fastq.gz (and _2)."""
    row = ena_filereport(accession, api)[0]
    problems = check_ena_metadata(row, platform, paired)
    # ENA lists FASTQ locations without a scheme (ftp.sra.ebi.ac.uk/vol1/...); the same paths are served over HTTPS
    urls = [u if "://" in u or u.startswith("/") else "https://" + u for u in row["fastq_ftp"].split(";") if u]
    if paired:
        # a paired run may also carry a third file of orphan reads; keep _1 and _2
        urls = sorted(u for u in urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz")))[:2] or urls[:2]
    else:
        urls = urls[:1]
    if not urls:
        raise RuntimeError(f"{accession}: " + "; ".join(problems))
    outputs = []
    for i, url in enumerate(urls, start=1):
        dest = f"{out_prefix}_{i}.fastq.gz"
        _stream_head(url, dest, n_reads, skip)
        outputs.append(dest)
    return {"label": accession, "kind": "ena", "accession": accession, "url": ";".join(urls),
            "path": ";".join(outputs), "bytes": sum(os.path.getsize(p) for p in outputs),
            "md5": "", "md5_expected": row.get("fastq_md5", ""), "md5_verified": False,
            "note": f"streamed reads {skip}..{skip + n_reads} only, so the archive MD5 of the whole file cannot be checked",
            "scientific_name": row.get("scientific_name", ""), "instrument_model": row.get("instrument_model", ""),
            "library_layout": row.get("library_layout", ""), "run_read_count": row.get("read_count", ""),
            "study_accession": row.get("study_accession", ""), "sample_accession": row.get("sample_accession", ""),
            "metadata_problems": problems, "retrieved_utc": now()}


# --------------------------------------------------------------------------- #
# remote BAM slicing
# --------------------------------------------------------------------------- #
def _bam_location(url: str) -> str:
    """Remote URLs pass through; local paths become absolute, because samtools
    is run from a scratch directory."""
    return os.path.abspath(local_path(url)) if is_local(url) else url


def bam_contigs(url: str) -> dict[str, int]:
    """Contig names and lengths from a (remote) BAM header."""
    url = _bam_location(url)
    result = subprocess.run(["samtools", "view", "-H", url], capture_output=True, text=True, cwd=tempfile.gettempdir())
    if result.returncode != 0:
        raise RuntimeError(f"cannot read BAM header of {url}\n{result.stderr.strip()}")
    contigs = {}
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            tags = dict(f.split(":", 1) for f in line.split("\t")[1:])
            contigs[tags["SN"]] = int(tags["LN"])
    return contigs


def choose_windows(fasta_path: str, contigs: list[str], n_windows: int, window_bp: int, seed: int,
                   max_n_fraction: float = 0.1) -> list[tuple[str, int, int]]:
    """Pick windows spread evenly over the chosen contigs, skipping gaps.

    The genome is cut into n equal segments and one window is drawn at a
    random position inside each. That covers every part of every contig
    (unlike purely random windows, which can cluster) while staying
    reproducible through the seed. Windows that are mostly N (centromere and
    short-arm gaps) are redrawn, since no read can come from them.
    """
    rng = np.random.default_rng(seed)
    fasta = pysam.FastaFile(fasta_path)
    lengths = {c: fasta.get_reference_length(c) for c in contigs}
    total = sum(lengths.values())
    windows = []
    for contig in contigs:
        length = lengths[contig]
        share = max(1, round(n_windows * length / total))
        size = min(window_bp, length)
        segment = length / share
        for i in range(share):
            low, high = int(i * segment), max(int(i * segment), int((i + 1) * segment) - size)
            for _ in range(20):
                start = int(rng.integers(low, high + 1))
                seq = CODE[to_array(fasta.fetch(contig, start, start + size))]
                if (seq == 4).mean() <= max_n_fraction:
                    windows.append((contig, start, start + size))
                    break
    return windows


def slice_bam(url: str, windows: list[tuple[str, int, int]], out_bam: str) -> dict:
    """Download only the alignments overlapping the given windows.

    Window contig names are translated into the BAM's own names first, since
    the provider may spell them differently from our FASTA (chr20 vs 20).
    Secondary and supplementary records are dropped: we want each read once.
    """
    url = _bam_location(url)
    header = bam_contigs(url)
    resolver = ContigResolver(header)
    regions = []
    for contig, start, end in windows:
        name = resolver.resolve(contig)
        if name is None:
            raise RuntimeError(f"contig {contig} not found in the BAM header of {url} "
                               f"(header starts with {list(header)[:5]})")
        regions.append(region_string(name, start, min(end, header[name])))
    os.makedirs(os.path.dirname(os.path.abspath(out_bam)), exist_ok=True)
    # run in a scratch directory: samtools saves the remote index into the working directory
    with tempfile.TemporaryDirectory() as scratch:
        cmd = ["samtools", "view", "-b", "-F", "0x900", "-o", os.path.abspath(out_bam), url, *regions]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=scratch)
        if result.returncode != 0:
            raise RuntimeError(f"samtools could not slice {url}\n{result.stderr.strip()}")
    return {"regions": regions, "bam_contig_names": resolver.report()}


def bam_to_fastq(bam_path: str, paired: bool, out_prefix: str, n_reads: int, seed: int) -> tuple[list[str], int]:
    """Turn a BAM slice back into raw FASTQ, as if it came off the sequencer.

    samtools fastq restores the original orientation of reverse-strand reads.
    For paired data only pairs with both mates inside the slice are kept.
    `collate` shuffles reads by a hash of their name, so taking the first
    n_reads afterwards is a random subsample rather than the left-most
    windows.
    """
    with tempfile.TemporaryDirectory() as scratch:
        collated = os.path.join(scratch, "collated.bam")
        subprocess.run(["samtools", "collate", "-u", "-o", collated, bam_path, os.path.join(scratch, "tmp")], check=True)
        if paired:
            r1, r2 = os.path.join(scratch, "r1.fq"), os.path.join(scratch, "r2.fq")
            subprocess.run(["samtools", "fastq", "-n", "-F", "0x900", "-1", r1, "-2", r2, "-0", "/dev/null",
                            "-s", "/dev/null", collated], check=True, capture_output=True)
            sources, outputs = [r1, r2], [f"{out_prefix}_1.fastq.gz", f"{out_prefix}_2.fastq.gz"]
        else:
            r0 = os.path.join(scratch, "r.fq")
            with open(r0, "w") as fh:
                subprocess.run(["samtools", "fastq", "-n", "-F", "0x900", collated], check=True, stdout=fh,
                               stderr=subprocess.DEVNULL)
            sources, outputs = [r0], [f"{out_prefix}_1.fastq.gz"]
        kept = 0
        for source, output in zip(sources, outputs):
            kept = 0
            with gzip.open(output, "wt", compresslevel=4) as out:
                for name, seq, qual in read_fastq(source):
                    if kept >= n_reads:
                        break
                    out.write(f"@{name}\n{seq}\n+\n{qual}\n")
                    kept += 1
    return outputs, kept


def provider_alignments(bam_path: str, out_parquet: str) -> int:
    """Keep where the data provider's own pipeline placed each read.

    GIAB aligned these reads against the whole genome with their own choice
    of aligner. That placement is not ground truth, but it is an independent
    opinion that our benchmark can be compared against in task 8.
    """
    names, mates, contigs, starts, ends, mapqs = [], [], [], [], [], []
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            names.append(read.query_name)
            mates.append(1 if read.is_read1 else 2 if read.is_read2 else 0)
            contigs.append(read.reference_name)
            starts.append(read.reference_start)
            ends.append(read.reference_end)
            mapqs.append(read.mapping_quality)
    table = pa.table({"qname": pa.array(names, pa.string()), "mate": pa.array(mates, pa.int8()),
                      "contig": pa.array(contigs, pa.string()).dictionary_encode(),
                      "start": pa.array(starts, pa.int64()), "end": pa.array(ends, pa.int64()),
                      "mapq": pa.array(mapqs, pa.int16())})
    pq.write_table(table, out_parquet, compression="zstd")
    return len(names)


def fetch_bam_slice(url: str, fasta_path: str, contigs: list[str], platform: str, paired: bool, n_reads: int,
                    n_windows: int, window_bp: int, seed: int, out_prefix: str, md5: str | None = None) -> dict:
    """Windows -> remote slice -> FASTQ subsample, with provenance."""
    contigs = list(contigs) or list(pysam.FastaFile(fasta_path).references)
    windows = choose_windows(fasta_path, contigs, n_windows, window_bp, seed)
    slice_path = out_prefix + ".slice.bam"
    info = slice_bam(url, windows, slice_path)
    n_provider = provider_alignments(slice_path, out_prefix + ".provider.parquet")
    outputs, kept = bam_to_fastq(slice_path, paired, out_prefix, n_reads, seed)
    record = {"label": os.path.basename(url), "kind": "bam_slice", "url": url, "path": ";".join(outputs),
              "bytes": os.path.getsize(slice_path), "md5": "", "md5_expected": md5 or "", "md5_verified": False,
              "note": f"{len(windows)} windows of {window_bp} bp sliced through the BAM index; "
                      "the archive MD5 covers the whole BAM and cannot be checked on a slice",
              "windows": [region_string(*w) for w in windows], "alignments_in_slice": n_provider,
              "reads_kept": kept, "bam_contig_names": info["bam_contig_names"], "platform": platform,
              "retrieved_utc": now()}
    os.remove(slice_path)
    return record


def fetch_local(paths: list[str], n_reads: int, out_prefix: str) -> dict:
    """Use FASTQ files that are already on disk (first n_reads of each)."""
    outputs = []
    for i, path in enumerate(paths, start=1):
        dest = f"{out_prefix}_{i}.fastq.gz"
        _stream_head(path, dest, n_reads, 0)
        outputs.append(dest)
    return {"label": os.path.basename(paths[0]), "kind": "local", "url": ";".join(paths), "path": ";".join(outputs),
            "bytes": sum(os.path.getsize(p) for p in outputs), "md5": ";".join(md5sum(p) for p in paths),
            "md5_expected": "", "md5_verified": False, "retrieved_utc": now()}


def merge_provenance(json_paths: list[str], out_tsv: str) -> None:
    """All provenance records -> one table for the report's data section."""
    columns = ["label", "kind", "url", "bytes", "md5", "md5_expected", "md5_verified", "retrieved_utc", "note",
               "scientific_name", "instrument_model", "library_layout", "study_accession", "sample_accession",
               "metadata_problems", "path"]
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        for path in sorted(json_paths):
            with open(path) as src:
                data = json.load(src)
            for record in data if isinstance(data, list) else [data]:
                row = []
                for column in columns:
                    value = record.get(column, "")
                    row.append("; ".join(value) if isinstance(value, list) else value)
                writer.writerow(row)
