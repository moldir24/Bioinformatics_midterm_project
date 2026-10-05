"""Alignments -> one row per read in a columnar table (task 6).

SAM/BAM is a good exchange format and a poor analysis format: answering "what
fraction of reads has MAPQ >= 30" means decoding every record again. Each
aligner's output is therefore reduced once to the handful of columns the
benchmark needs and stored as Parquet, which is compressed, typed and loads
column by column.

Only the primary alignment of each read is kept:
  * secondary alignments (flag 0x100) are the aligner's alternative guesses
  * supplementary alignments (flag 0x800) are the other pieces of a split read
The primary record is the aligner's single answer to "where is this read
from", and its MAPQ is the confidence it attaches to that answer.
"""
from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

SCHEMA = pa.schema([
    ("qname", pa.string()), ("mate", pa.int8()), ("mapped", pa.bool_()),
    ("contig", pa.string()), ("start", pa.int64()), ("end", pa.int64()), ("strand", pa.string()),
    ("mapq", pa.int16()), ("nm", pa.int32()), ("read_len", pa.int32()), ("clipped", pa.int32()),
    ("proper", pa.bool_()), ("split", pa.bool_()),
])


def bam_to_table(bam_path: str, out_parquet: str) -> dict:
    """Reduce a BAM/SAM to a Parquet table of primary alignments.

    Coordinates follow pysam, i.e. the BED convention: 0-based start,
    exclusive end, both on the reference.
    """
    cols: dict[str, list] = {name: [] for name in SCHEMA.names}
    seen = set()
    n_records = n_duplicate = 0
    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            n_records += 1
            if read.is_secondary or read.is_supplementary:
                continue
            mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
            key = (read.query_name, mate)
            if key in seen:  # an aligner must not emit two primaries for one read
                n_duplicate += 1
                continue
            seen.add(key)
            mapped = not read.is_unmapped
            cols["qname"].append(read.query_name)
            cols["mate"].append(mate)
            cols["mapped"].append(mapped)
            cols["proper"].append(bool(read.is_proper_pair))
            if mapped:
                cigar = read.cigartuples or []
                clipped = sum(length for op, length in cigar if op in (4, 5))  # soft + hard clips
                cols["contig"].append(read.reference_name)
                cols["start"].append(read.reference_start)
                cols["end"].append(read.reference_end)
                cols["strand"].append("-" if read.is_reverse else "+")
                cols["mapq"].append(read.mapping_quality)
                cols["nm"].append(read.get_tag("NM") if read.has_tag("NM") else -1)
                cols["read_len"].append(read.infer_read_length() or 0)
                cols["clipped"].append(clipped)
                cols["split"].append(read.has_tag("SA"))
            else:
                cols["contig"].append("")
                cols["start"].append(-1)
                cols["end"].append(-1)
                cols["strand"].append("")
                cols["mapq"].append(0)
                cols["nm"].append(-1)
                cols["read_len"].append(read.query_length or 0)
                cols["clipped"].append(0)
                cols["split"].append(False)
    table = pa.table({name: pa.array(cols[name], SCHEMA.field(name).type) for name in SCHEMA.names})
    pq.write_table(table, out_parquet, compression="zstd")
    return {"records": n_records, "primary_reads": len(seen), "duplicate_primaries": n_duplicate}
