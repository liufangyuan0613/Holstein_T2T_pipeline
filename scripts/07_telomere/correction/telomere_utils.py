#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared utilities for telomere length correction.

Coordinate convention:
- BED-like coordinates are treated as 0-based, half-open: [start, end).
- FASTA sequence slicing uses the same convention.
"""

from __future__ import annotations

import csv
import gzip
import re
from typing import Dict, Iterable, Iterator, List, Tuple, Optional


NA_VALUES = {"", "NA", "NaN", "nan", "None", "null"}


def smart_open(path: str, mode: str = "rt"):
    """Open plain text or gzip-compressed files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_int(value, default: Optional[int] = None) -> Optional[int]:
    if value is None:
        return default
    s = str(value).strip()
    if s in NA_VALUES:
        return default
    return int(float(s))


def parse_float(value, default: Optional[float] = None) -> Optional[float]:
    if value is None:
        return default
    s = str(value).strip()
    if s in NA_VALUES:
        return default
    return float(s)


def read_tsv(path: str) -> List[dict]:
    with smart_open(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def write_tsv(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_fai(path: str) -> Dict[str, int]:
    fai: Dict[str, int] = {}
    with smart_open(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            fai[parts[0]] = int(parts[1])
    return fai


def contig_to_telogator_chr(contig: str, contig_prefix: str = "cow_chr") -> str:
    """
    Convert a FASTA contig name to telogator-style chromosome name.

    Examples:
      cow_chr01 -> chr1
      cow_chr1  -> chr1
      cow_chrX  -> chrX
      chr01     -> chr1
      chrX      -> chrX
    """
    x = contig
    if contig_prefix and x.startswith(contig_prefix):
        x = x[len(contig_prefix):]
        return "chr" + normalize_chr_suffix(x)
    if x.startswith("chr"):
        return "chr" + normalize_chr_suffix(x[3:])
    return "chr" + normalize_chr_suffix(x)


def normalize_chr_suffix(s: str) -> str:
    s = str(s)
    if re.fullmatch(r"\d+", s):
        return str(int(s))
    return s


def repeat_to_length(motif: str, length: int) -> str:
    if length <= 0:
        return ""
    n = (length + len(motif) - 1) // len(motif)
    return (motif * n)[:length]


def fasta_iter(path: str) -> Iterator[Tuple[str, str, str]]:
    """
    Yield (name, description, sequence) from a FASTA file.
    name excludes the leading '>' and stops at the first whitespace.
    description is the complete header without '>'.
    """
    name = None
    desc = None
    chunks: List[str] = []
    with smart_open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks).upper()
                desc = line[1:]
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, desc, "".join(chunks).upper()


def write_fasta_record(fh, desc: str, seq: str, width: int = 60) -> None:
    fh.write(f">{desc}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def count_prefix_repeat(seq: str, motif: str) -> int:
    motif = motif.upper()
    seq = seq.upper()
    i = 0
    m = len(motif)
    while i + m <= len(seq) and seq[i:i+m] == motif:
        i += m
    # Allow final partial motif at the boundary.
    rem = min(m - 1, len(seq) - i)
    for k in range(rem, 0, -1):
        if seq[i:i+k] == motif[:k]:
            return i + k
    return i


def count_suffix_repeat(seq: str, motif: str) -> int:
    motif = motif.upper()
    seq = seq.upper()
    rc = motif
    i = len(seq)
    m = len(rc)
    count = 0
    while i - m >= 0 and seq[i-m:i] == rc:
        i -= m
        count += m
    # Allow final partial motif at the boundary on the left side of the suffix run.
    rem = min(m - 1, i)
    for k in range(rem, 0, -1):
        if seq[i-k:i] == rc[-k:]:
            return count + k
    return count


def classify_error(abs_error: Optional[int], pass_bp: int = 6, warn_bp: int = 100) -> str:
    if abs_error is None:
        return "MISSING"
    if abs_error <= pass_bp:
        return "PASS"
    if abs_error <= warn_bp:
        return "WARN"
    return "FAIL"
