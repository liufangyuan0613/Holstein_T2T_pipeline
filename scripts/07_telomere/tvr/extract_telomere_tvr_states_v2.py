#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
"""
extract_telomere_tvr_states_v2.py
Assembly-based telomere TVR map builder with improved non-C sensitive bin calling.
"""
import argparse, csv, shutil, subprocess, sys
from collections import Counter
from pathlib import Path

DEFAULT_MOTIF_TO_STATE = {
    "TTAGGG": "C",
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
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
REQ_COLS = [
    "chr",
    "5_telomere_start", "5_telomere_end", "5_tele_length",
    "3_telomere_start", "3_telomere_end", "3_tele_length",
]

def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]

def check_samtools():
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH.")

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

def normalize_kmer_to_g_rich(kmer: str, motif_map: dict[str, str]):
    kmer = kmer.upper()
    if kmer in motif_map:
        return kmer
    rc = reverse_complement(kmer)
    if rc in motif_map:
        return rc
    return None

def choose_best_phase(seq: str, motif_map: dict[str, str]) -> int:
    best_offset = 0
    best_score = (-1, -1)
    for offset in range(6):
        recognized = 0
        total = 0
        for i in range(offset, len(seq) - 5, 6):
            total += 1
            if normalize_kmer_to_g_rich(seq[i:i+6], motif_map) is not None:
                recognized += 1
        score = (recognized, total)
        if score > best_score:
            best_score = score
            best_offset = offset
    return best_offset

def build_phased_state_series(seq: str, motif_map: dict[str, str], other_state: str = "O"):
    if len(seq) < 6:
        return [], 0
    offset = choose_best_phase(seq, motif_map)
    out = []
    for i in range(offset, len(seq) - 5, 6):
        k = seq[i:i+6]
        motif = normalize_kmer_to_g_rich(k, motif_map)
        if motif is None:
            state = other_state
            motif = k
        else:
            state = motif_map[motif]
        out.append((i, i + 6, motif, state))
    return out, offset

def summarize_bin(series, bin_start: int, bin_end: int):
    hits = []
    for s, e, motif, state in series:
        mid = (s + e) / 2.0
        if bin_start <= mid < bin_end:
            hits.append(state)
    return Counter(hits)

def call_bin_state_sensitive(counts: Counter, canonical_state: str = "C", empty_state: str = "_",
                             nonc_threshold: float = 0.15, min_kmers_per_bin: int = 1) -> str:
    total = sum(counts.values())
    if total < min_kmers_per_bin:
        return empty_state
    nonc_counts = {k: v for k, v in counts.items() if k != canonical_state}
    nonc_total = sum(nonc_counts.values())
    if nonc_total > 0 and (nonc_total / total) >= nonc_threshold:
        return sorted(nonc_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]
    return sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0][0]

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
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--bin-bp", type=int, default=20)
    ap.add_argument("--max-bp", type=int, default=15000)
    ap.add_argument("--other-state", default="O")
    ap.add_argument("--nonc-threshold", type=float, default=0.15)
    ap.add_argument("--min-kmers-per-bin", type=int, default=1)
    return ap.parse_args()

def main():
    args = parse_args()
    check_samtools()
    motif_map = DEFAULT_MOTIF_TO_STATE.copy()
    rows_state, rows_sum = [], []
    with open(args.bed, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in REQ_COLS if c not in reader.fieldnames]
        if missing:
            raise RuntimeError(f"BED-like file missing required columns: {missing}")
        for row in reader:
            chrom = row["chr"]
            specs = [
                {"end_label": "5'", "start_col": "5_telomere_start", "len_col": "5_tele_length", "need_revcomp": True},
                {"end_label": "3'", "start_col": "3_telomere_start", "len_col": "3_tele_length", "need_revcomp": False},
            ]
            for spec in specs:
                start_1based = int(float(row[spec["start_col"]]))
                tel_len = int(float(row[spec["len_col"]]))
                seq = fetch_seq_by_samtools(args.fasta, chrom, start_1based, tel_len)
                if spec["need_revcomp"]:
                    seq = reverse_complement(seq)
                series, phase_offset = build_phased_state_series(seq, motif_map, other_state=args.other_state)
                counts_all = Counter([x[3] for x in series]) if series else Counter()
                state_row = {
                    "chrom": chrom,
                    "end_label": spec["end_label"],
                    "chr_end": f"{chrom} {spec['end_label']}",
                    "telomere_length_bp": tel_len,
                    "phase_offset_bp": phase_offset,
                    "n_phased_6mers": len(series),
                }
                for bin_start in range(0, args.max_bp, args.bin_bp):
                    if bin_start >= tel_len:
                        state = "_"
                    else:
                        bin_end = min(bin_start + args.bin_bp, args.max_bp)
                        counts = summarize_bin(series, bin_start, bin_end)
                        state = call_bin_state_sensitive(
                            counts,
                            canonical_state="C",
                            empty_state="_",
                            nonc_threshold=args.nonc_threshold,
                            min_kmers_per_bin=args.min_kmers_per_bin
                        )
                    state_row[f"bin_{bin_start}"] = state
                rows_state.append(state_row)
                summary_row = {
                    "chrom": chrom,
                    "end_label": spec["end_label"],
                    "chr_end": f"{chrom} {spec['end_label']}",
                    "telomere_length_bp": tel_len,
                    "phase_offset_bp": phase_offset,
                    "n_phased_6mers": len(series),
                    "fraction_canonical_C": (counts_all.get("C", 0) / len(series)) if series else 0.0,
                    "fraction_other_noncanonical": (counts_all.get(args.other_state, 0) / len(series)) if series else 0.0,
                    "nonc_threshold": args.nonc_threshold,
                    "bin_bp": args.bin_bp,
                }
                for st in sorted(set(motif_map.values()) | {args.other_state}):
                    summary_row[f"count_state_{st}"] = counts_all.get(st, 0)
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
    sys.stderr.write(f"[INFO] non-C sensitive caller used with threshold={args.nonc_threshold}, bin_bp={args.bin_bp}\n")

if __name__ == "__main__":
    main()
