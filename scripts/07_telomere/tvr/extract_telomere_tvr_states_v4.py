#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_telomere_tvr_states_v4.py

Assembly-based telomere TVR map builder with improved TVR recovery.

Why this version
----------------
If v2 still produces too many "O" bins and misses TVR signal, the usual causes are:
1) exact-match motif dictionary is too sparse
2) fixed 6-mer phasing misses local phase drift / interruptions
3) dominant-state binning lets C swallow mixed bins

This version improves all three points:
- uses overlapping 6-mer scanning (not only one global phase)
- supports fuzzy motif matching by Hamming distance
- uses telomeric-support-aware bin calling:
    * if enough recognized telomeric 6-mers exist in a bin, preserve TVR states
    * reserve O for bins with poor telomeric support
- prioritizes dominant non-C state when non-C burden is sufficiently high

Input
-----
--fasta genome.fa
--bed telomere BED-like file with columns:
    chr
    5_telomere_start  5_telomere_end  5_tele_length
    3_telomere_start  3_telomere_end  3_tele_length

Output
------
<prefix>.state_matrix.tsv
<prefix>.summary.tsv
"""

from __future__ import annotations
import argparse
import csv
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

DEFAULT_MOTIF_TO_STATE = {
    "TTAGGG": "C",   # canonical
    "TCAGGG": "A",
    "TGAGGG": "D",
    "TTGGGG": "F",
    "TTAGGC": "H",
    "TTAGAG": "L",
    "TAAGGG": "N",
    "CTAGGG": "S",
    "GTAGGG": "T",
    "TTAGCG": "V",
}

REQ_COLS = [
    "chr",
    "5_telomere_start", "5_telomere_end", "5_tele_length",
    "3_telomere_start", "3_telomere_end", "3_tele_length",
]

RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def check_samtools() -> None:
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH. Please install samtools.")


def fetch_seq_by_samtools(fasta: str, chrom: str, start_1based: int, length_bp: int) -> str:
    if length_bp <= 0:
        return ""
    end_1based = start_1based + length_bp - 1
    region = f"{chrom}:{start_1based}-{end_1based}"
    res = subprocess.run(["samtools", "faidx", fasta, region], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {region}:\n{res.stderr}")
    lines = [x.strip() for x in res.stdout.splitlines() if x.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"Unexpected faidx output for {region}")
    return "".join(lines[1:]).upper()


def classify_kmer_fuzzy(kmer: str,
                        motif_map: dict[str, str],
                        max_mismatch: int = 1):
    """
    Return (best_motif, best_state, mismatch_count) or None.

    Strategy:
    - consider both kmer and its reverse complement
    - choose the nearest motif by Hamming distance
    - require distance <= max_mismatch
    - if tie among different motifs/states at same distance, return None (ambiguous)
    """
    kmer = kmer.upper()
    candidates = [kmer, reverse_complement(kmer)]

    best = None
    best_dist = None
    best_states = set()

    for q in candidates:
        for motif, st in motif_map.items():
            d = hamming(q, motif)
            if best_dist is None or d < best_dist:
                best_dist = d
                best = (motif, st, d)
                best_states = {st}
            elif d == best_dist:
                best_states.add(st)

    if best_dist is None or best_dist > max_mismatch:
        return None
    if len(best_states) > 1:
        return None
    return best


def build_overlapping_state_series(seq: str,
                                   motif_map: dict[str, str],
                                   max_mismatch: int = 1):
    """
    Return list of tuples:
      (start_bp_0based, end_bp_0based, motif_or_raw, state_or_None, mismatch_or_None)
    using overlapping 6-mers across the oriented telomere sequence.
    """
    out = []
    if len(seq) < 6:
        return out

    for i in range(0, len(seq) - 5):
        k = seq[i:i+6]
        hit = classify_kmer_fuzzy(k, motif_map, max_mismatch=max_mismatch)
        if hit is None:
            out.append((i, i + 6, k, None, None))
        else:
            motif, st, d = hit
            out.append((i, i + 6, motif, st, d))
    return out


def summarize_bin(series, bin_start: int, bin_end: int):
    """
    Summarize overlapping 6-mers whose midpoint falls within this bin.
    Returns:
      total_kmers
      recognized_kmers
      state_counts (recognized only)
    """
    total_kmers = 0
    recognized_kmers = 0
    state_counts = Counter()

    for s, e, motif_or_raw, state, mismatch in series:
        mid = (s + e) / 2.0
        if bin_start <= mid < bin_end:
            total_kmers += 1
            if state is not None:
                recognized_kmers += 1
                state_counts[state] += 1

    return total_kmers, recognized_kmers, state_counts


def call_bin_state(total_kmers: int,
                   recognized_kmers: int,
                   state_counts: Counter,
                   canonical_state: str = "C",
                   empty_state: str = "_",
                   other_state: str = "O",
                   min_kmers_per_bin: int = 3,
                   telomeric_support_threshold: float = 0.20,
                   nonc_threshold: float = 0.15) -> str:
    """
    Bin calling logic:
    1) if too few kmers -> "_"
    2) if no recognized telomeric kmers -> "O"
    3) if telomeric support < threshold -> "O"
    4) if non-C fraction among recognized kmers >= nonc_threshold:
         call dominant non-C state
    5) else call dominant recognized state
    """
    if total_kmers < min_kmers_per_bin:
        return empty_state
    if recognized_kmers == 0:
        return other_state

    telomeric_support = recognized_kmers / total_kmers
    if telomeric_support < telomeric_support_threshold:
        return other_state

    total_rec = sum(state_counts.values())
    nonc_counts = {k: v for k, v in state_counts.items() if k != canonical_state}
    nonc_total = sum(nonc_counts.values())

    if total_rec > 0 and nonc_total > 0 and (nonc_total / total_rec) >= nonc_threshold:
        return sorted(nonc_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]

    return sorted(state_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]


def chr_key(ch: str):
    core = ch.replace("chr", "")
    try:
        return (0, int(core))
    except ValueError:
        if core == "X":
            return (1, 1000)
        if core == "Y":
            return (1, 1001)
        return (2, core)


def parse_args():
    ap = argparse.ArgumentParser(description="Build assembly-based TVR state matrix with improved TVR recovery (v4)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--bin-bp", type=int, default=20)
    ap.add_argument("--max-bp", type=int, default=15000)
    ap.add_argument("--max-mismatch", type=int, default=1,
                    help="Maximum Hamming distance for fuzzy motif match (default: 1)")
    ap.add_argument("--nonc-threshold", type=float, default=0.15,
                    help="If non-C fraction among recognized kmers >= this, call dominant non-C state (default: 0.15)")
    ap.add_argument("--telomeric-support-threshold", type=float, default=0.20,
                    help="Minimum recognized/total 6-mer fraction required to avoid O (default: 0.20)")
    ap.add_argument("--min-kmers-per-bin", type=int, default=3,
                    help="Minimum overlapping 6-mers needed to emit a robust bin state (default: 3)")
    return ap.parse_args()


def main():
    args = parse_args()
    check_samtools()

    motif_map = DEFAULT_MOTIF_TO_STATE.copy()
    rows_state = []
    rows_sum = []

    with open(args.bed, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in REQ_COLS if c not in reader.fieldnames]
        if missing:
            raise RuntimeError(f"BED-like file missing required columns: {missing}")

        for row in reader:
            chrom = row["chr"]

            end_specs = [
                {"end_label": "5'", "start_col": "5_telomere_start", "len_col": "5_tele_length", "need_revcomp": True},
                {"end_label": "3'", "start_col": "3_telomere_start", "len_col": "3_tele_length", "need_revcomp": False},
            ]

            for spec in end_specs:
                start_1based = int(float(row[spec["start_col"]]))
                tel_len = int(float(row[spec["len_col"]]))

                seq = fetch_seq_by_samtools(args.fasta, chrom, start_1based, tel_len)
                if spec["need_revcomp"]:
                    seq = reverse_complement(seq)

                series = build_overlapping_state_series(seq, motif_map, max_mismatch=args.max_mismatch)

                # overall summary from recognized overlapping 6-mers
                overall_counts = Counter([x[3] for x in series if x[3] is not None])

                state_row = {
                    "chrom": chrom,
                    "end_label": spec["end_label"],
                    "chr_end": f"{chrom} {spec['end_label']}",
                    "telomere_length_bp": tel_len,
                    "n_overlapping_6mers": len(series),
                    "bin_bp": args.bin_bp,
                    "max_mismatch": args.max_mismatch,
                    "nonc_threshold": args.nonc_threshold,
                    "telomeric_support_threshold": args.telomeric_support_threshold,
                }

                for bin_start in range(0, args.max_bp, args.bin_bp):
                    if bin_start >= tel_len:
                        state = "_"
                    else:
                        bin_end = min(bin_start + args.bin_bp, args.max_bp)
                        total_kmers, recognized_kmers, state_counts = summarize_bin(series, bin_start, bin_end)
                        state = call_bin_state(
                            total_kmers=total_kmers,
                            recognized_kmers=recognized_kmers,
                            state_counts=state_counts,
                            canonical_state="C",
                            empty_state="_",
                            other_state="O",
                            min_kmers_per_bin=args.min_kmers_per_bin,
                            telomeric_support_threshold=args.telomeric_support_threshold,
                            nonc_threshold=args.nonc_threshold
                        )
                    state_row[f"bin_{bin_start}"] = state

                rows_state.append(state_row)

                summary_row = {
                    "chrom": chrom,
                    "end_label": spec["end_label"],
                    "chr_end": f"{chrom} {spec['end_label']}",
                    "telomere_length_bp": tel_len,
                    "n_overlapping_6mers": len(series),
                    "count_recognized_telomeric_6mers": sum(overall_counts.values()),
                    "fraction_canonical_C": (overall_counts.get("C", 0) / sum(overall_counts.values())) if overall_counts else 0.0,
                    "bin_bp": args.bin_bp,
                    "max_mismatch": args.max_mismatch,
                    "nonc_threshold": args.nonc_threshold,
                    "telomeric_support_threshold": args.telomeric_support_threshold,
                }
                for st in sorted(set(motif_map.values())):
                    summary_row[f"count_state_{st}"] = overall_counts.get(st, 0)

                rows_sum.append(summary_row)

    rows_state = sorted(rows_state, key=lambda r: (chr_key(r["chrom"]), 0 if r["end_label"] == "5'" else 1))
    rows_sum = sorted(rows_sum, key=lambda r: (chr_key(r["chrom"]), 0 if r["end_label"] == "5'" else 1))

    state_out = Path(f"{args.out_prefix}.state_matrix.tsv")
    sum_out = Path(f"{args.out_prefix}.summary.tsv")

    with open(state_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_state[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_state)

    with open(sum_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_sum[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_sum)

    sys.stderr.write(f"[OK] Wrote: {state_out}\n")
    sys.stderr.write(f"[OK] Wrote: {sum_out}\n")
    sys.stderr.write(
        "[INFO] v4 used overlapping 6-mer scanning + fuzzy motif matching + telomeric-support-aware bin calling\n"
    )


if __name__ == "__main__":
    main()
