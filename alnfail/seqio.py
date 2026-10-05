"""Sequence helpers: base encoding, reverse complement, FASTQ in and out.

Sequences are handled as numpy uint8 arrays of ASCII codes. That keeps the
simulator and the profiler vectorised (one numpy call per read instead of one
Python operation per base), which is what makes a million reads feasible in
plain Python.
"""
from __future__ import annotations

import gzip
import io
import subprocess
from contextlib import contextmanager

import numpy as np

A, C, G, T, N = (ord(c) for c in "ACGTN")

# ASCII -> 0..3 for ACGT (upper or lower case), 4 for everything else
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    CODE[ord(_c)] = _i
    CODE[ord(_c.lower())] = _i
BASES = np.frombuffer(b"ACGTN", dtype=np.uint8)

# ASCII complement table (IUPAC-aware enough for ACGTN, case preserved as upper)
_COMP = np.arange(256, dtype=np.uint8)
for _a, _b in zip("ACGTNacgtn", "TGCANTGCAN"):
    _COMP[ord(_a)] = ord(_b)

_UPPER = np.arange(256, dtype=np.uint8)
for _c in "acgtn":
    _UPPER[ord(_c)] = ord(_c.upper())


def to_array(seq: str | bytes) -> np.ndarray:
    """Sequence string -> upper-case uint8 array (a writable copy)."""
    if isinstance(seq, str):
        seq = seq.encode("ascii")
    return _UPPER[np.frombuffer(seq, dtype=np.uint8)]


def to_str(arr: np.ndarray) -> str:
    return arr.tobytes().decode("ascii")


def revcomp(arr: np.ndarray) -> np.ndarray:
    """Reverse complement of a uint8 base array."""
    return _COMP[arr[::-1]]


def gc_fraction(arr: np.ndarray) -> float:
    """GC fraction among A/C/G/T bases (N ignored). NaN for an empty sequence."""
    codes = CODE[arr]
    acgt = int((codes < 4).sum())
    if acgt == 0:
        return float("nan")
    return float(((codes == 1) | (codes == 2)).sum()) / acgt


def homopolymer_run_lengths(arr: np.ndarray) -> np.ndarray:
    """For every base, the length of the homopolymer run it belongs to.

    AAACG -> [3, 3, 3, 1, 1]. Used both to learn and to reproduce the
    homopolymer-dependent indel errors of long-read platforms.
    """
    n = arr.size
    if n == 0:
        return np.zeros(0, dtype=np.int32)
    change = np.flatnonzero(arr[1:] != arr[:-1]) + 1
    starts = np.concatenate(([0], change))
    lengths = np.diff(np.concatenate((starts, [n])))
    return np.repeat(lengths, lengths).astype(np.int32)


# --------------------------------------------------------------------------- #
# FASTQ
# --------------------------------------------------------------------------- #
def open_text(path: str, mode: str = "rt"):
    """Open a plain or gzip-compressed text file by extension."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fastq(path: str):
    """Yield (name, sequence, quality) from a FASTQ(.gz), streaming.

    Reads are never all held in memory: the caller decides what to keep. The
    name is the first whitespace-delimited token without the leading '@'.
    """
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual and not seq:
                return
            yield header[1:].split()[0], seq, qual


@contextmanager
def fastq_writer(path: str, level: int = 4):
    """Context manager yielding a text handle that writes gzip-compressed FASTQ.

    Uses the external `gzip`/`pigz` when available because Python's gzip
    module is several times slower, and falls back to it otherwise.
    """
    proc = None
    raw = None
    try:
        if str(path).endswith(".gz"):
            tool = None
            for candidate in ("pigz", "gzip"):
                if subprocess.run(["which", candidate], capture_output=True).returncode == 0:
                    tool = candidate
                    break
            if tool:
                raw = open(path, "wb")
                proc = subprocess.Popen([tool, f"-{level}", "-c"], stdin=subprocess.PIPE, stdout=raw)
                handle = io.TextIOWrapper(proc.stdin, encoding="ascii", write_through=False)
            else:
                handle = gzip.open(path, "wt", compresslevel=level)
        else:
            handle = open(path, "wt")
        yield handle
    finally:
        try:
            handle.close()
        finally:
            if proc is not None:
                rc = proc.wait()
                raw.close()
                if rc != 0:
                    raise RuntimeError(f"compressor exited with status {rc} while writing {path}")


def write_record(handle, name: str, seq: np.ndarray, qual: np.ndarray) -> None:
    """Write one FASTQ record. `qual` holds Phred values (not ASCII)."""
    handle.write("@")
    handle.write(name)
    handle.write("\n")
    handle.write(seq.tobytes().decode("ascii"))
    handle.write("\n+\n")
    handle.write((qual.astype(np.uint8) + 33).tobytes().decode("ascii"))
    handle.write("\n")
