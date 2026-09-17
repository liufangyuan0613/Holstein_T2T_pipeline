#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse

from telomere_utils import (
    read_tsv, write_tsv, fasta_iter, parse_int,
    count_prefix_repeat, count_suffix_repeat, classify_error
)


def main():
    ap = argparse.ArgumentParser(description="Evaluate exact terminal motif lengths in corrected FASTA.")
    ap.add_argument("--corrected-fa", required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--out", default="corrected_telomere_exact_motif_evaluation.tsv")
    ap.add_argument("--pass-bp", type=int, default=6)
    ap.add_argument("--warn-bp", type=int, default=100)
    args = ap.parse_args()

    plan = [r for r in read_tsv(args.plan) if r.get("status") == "OK"]
    plan_by_contig = {}
    for r in plan:
        plan_by_contig.setdefault(r["contig"], []).append(r)

    rows = []
    observed_contigs = set()

    for name, desc, seq in fasta_iter(args.corrected_fa):
        if name not in plan_by_contig:
            continue
        observed_contigs.add(name)
        for r in plan_by_contig[name]:
            motif = r["motif"].upper()
            target = parse_int(r["target_telomere_length"])
            if r["end_side"] == "left":
                obs = count_prefix_repeat(seq, motif)
            else:
                obs = count_suffix_repeat(seq, motif)
            err = obs - target
            abs_err = abs(err)
            rows.append({
                "contig": name,
                "arm": r["arm"],
                "chr_arm": r["chr_arm"],
                "end_side": r["end_side"],
                "motif": motif,
                "target_telomere_length": target,
                "observed_terminal_exact_motif_length": obs,
                "error_bp": err,
                "abs_error_bp": abs_err,
                "status": classify_error(abs_err, args.pass_bp, args.warn_bp),
                "replace_mode": r.get("replace_mode", ""),
            })

    for contig in sorted(set(plan_by_contig) - observed_contigs):
        for r in plan_by_contig[contig]:
            rows.append({
                "contig": contig,
                "arm": r["arm"],
                "chr_arm": r["chr_arm"],
                "end_side": r["end_side"],
                "motif": r["motif"],
                "target_telomere_length": r["target_telomere_length"],
                "observed_terminal_exact_motif_length": "",
                "error_bp": "",
                "abs_error_bp": "",
                "status": "MISSING_CONTIG_IN_FASTA",
                "replace_mode": r.get("replace_mode", ""),
            })

    fields = [
        "contig", "arm", "chr_arm", "end_side", "motif",
        "target_telomere_length", "observed_terminal_exact_motif_length",
        "error_bp", "abs_error_bp", "status", "replace_mode"
    ]
    write_tsv(args.out, rows, fields)

    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    print(f"Wrote motif evaluation: {args.out}")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
