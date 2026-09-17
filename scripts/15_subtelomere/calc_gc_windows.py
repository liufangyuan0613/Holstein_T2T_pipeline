#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

def fasta_iter(path):
    """Yield (name, seq) from a FASTA file. name is first token in header."""
    name = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            yield name, "".join(seq_chunks)

def gc_ratio(seq: str, ignore_n: bool = True) -> float:
    """
    GC ratio in a window.
    If ignore_n=True: denominator excludes non-ACGT (e.g., N). If denom==0 -> NaN.
    Else: denominator is window length.
    """
    s = seq.upper()
    g = s.count("G")
    c = s.count("C")
    if ignore_n:
        a = s.count("A")
        t = s.count("T")
        denom = a + t + g + c
        if denom == 0:
            return float("nan")
        return (g + c) / denom
    else:
        if len(s) == 0:
            return float("nan")
        return (g + c) / len(s)

def main():
    ap = argparse.ArgumentParser(description="Calculate GC ratio in fixed windows from a genome FASTA.")
    ap.add_argument("-f", "--fasta", required=True, help="Genome FASTA")
    ap.add_argument("-o", "--out", required=True, help="Output TSV: chr start end gc_ratio")
    ap.add_argument("-w", "--window", type=int, default=100000, help="Window size (default 100000)")
    ap.add_argument("--step", type=int, default=None, help="Step size (default = window, i.e., non-overlapping)")
    ap.add_argument("--ignore-n", action="store_true",
                    help="Ignore Ns/non-ACGT in denominator (recommended).")
    ap.add_argument("--one-based", action="store_true",
                    help="Output coordinates as 1-based inclusive. Default: 0-based half-open [start,end).")
    args = ap.parse_args()

    w = args.window
    step = args.step if args.step is not None else w

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("chr\tstart\tend\tgc_ratio\n")
        for chrom, seq in fasta_iter(args.fasta):
            L = len(seq)
            # iterate windows
            for start in range(0, L, step):
                end = min(start + w, L)
                win = seq[start:end]
                gc = gc_ratio(win, ignore_n=args.ignore_n)

                if args.one_based:
                    # 1-based inclusive coordinates
                    s_out = start + 1
                    e_out = end
                else:
                    # 0-based half-open coordinates
                    s_out = start
                    e_out = end

                # print NaN as NA (R-friendly)
                if gc != gc:  # NaN check
                    out.write(f"{chrom}\t{s_out}\t{e_out}\tNA\n")
                else:
                    out.write(f"{chrom}\t{s_out}\t{e_out}\t{gc:.6f}\n")

if __name__ == "__main__":
    main()
