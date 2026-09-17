#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
from collections import Counter

from telomere_utils import (
    read_fai, read_tsv, write_tsv, parse_int, parse_float,
    contig_to_telogator_chr
)


def main():
    ap = argparse.ArgumentParser(
        description="Check consistency among FASTA .fai, current telomere BED, and telogator2 arm_stats.tsv."
    )
    ap.add_argument("--fai", required=True, help="FASTA index, e.g. cow_t2t.fa.fai")
    ap.add_argument("--tele-bed", required=True, help="Current telomere BED-like table.")
    ap.add_argument("--arm-stats", required=True, help="telogator2 summarized arm_stats.tsv.")
    ap.add_argument("--out-prefix", default="telomere_input_qc", help="Output prefix.")
    ap.add_argument("--target-column", default="plot_value_bp", help="Target telomere length column in arm_stats.tsv.")
    ap.add_argument("--contig-prefix", default="cow_chr", help="Prefix used in FASTA contig names.")
    ap.add_argument("--terminal-slack-bp", type=int, default=1000,
                    help="Slack for reporting whether detected telomere block is near the terminal boundary.")
    args = ap.parse_args()

    fai = read_fai(args.fai)
    tele = read_tsv(args.tele_bed)
    arms = read_tsv(args.arm_stats)

    target_by_arm = {}
    for r in arms:
        arm_id = r.get("chr_arm") or ((r.get("chr") or "") + (r.get("arm") or ""))
        target = parse_int(r.get(args.target_column))
        if arm_id and target is not None:
            target_by_arm[arm_id] = r

    tele_by_chr = {r["chr"]: r for r in tele}

    coord_rows = []
    for contig, chrom_len in fai.items():
        r = tele_by_chr.get(contig)
        row = {
            "contig": contig,
            "fai_length": chrom_len,
            "tele_bed_found": "yes" if r else "no",
            "status": "",
            "notes": "",
        }
        issues = []
        if r is None:
            issues.append("MISSING_IN_TELE_BED")
        else:
            vals = {}
            for col in [
                "5_telomere_start", "5_telomere_end", "5_tele_length",
                "3_telomere_start", "3_telomere_end", "3_tele_length"
            ]:
                vals[col] = parse_int(r.get(col))
                row[col] = vals[col]

            if not (0 <= vals["5_telomere_start"] <= vals["5_telomere_end"] <= chrom_len):
                issues.append("INVALID_5_COORD")
            if not (0 <= vals["3_telomere_start"] <= vals["3_telomere_end"] <= chrom_len):
                issues.append("INVALID_3_COORD")
            if vals["5_telomere_end"] - vals["5_telomere_start"] != vals["5_tele_length"]:
                issues.append("MISMATCH_5_LENGTH")
            if vals["3_telomere_end"] - vals["3_telomere_start"] != vals["3_tele_length"]:
                issues.append("MISMATCH_3_LENGTH")

            left_gap = vals["5_telomere_start"]
            right_gap = chrom_len - vals["3_telomere_end"]
            row["left_non_telomeric_prefix_bp"] = left_gap
            row["right_non_telomeric_suffix_bp"] = right_gap

            if left_gap > args.terminal_slack_bp:
                issues.append("LEFT_TELOMERE_NOT_NEAR_CONTIG_START")
            if right_gap > args.terminal_slack_bp:
                issues.append("RIGHT_TELOMERE_NOT_NEAR_CONTIG_END")

        row["status"] = "PASS" if not issues else ";".join(issues)
        coord_rows.append(row)

    # telogator target coverage
    missing_rows = []
    target_rows = []
    for contig in fai:
        tel_chr = contig_to_telogator_chr(contig, args.contig_prefix)
        for arm in ["p", "q"]:
            chr_arm = f"{tel_chr}{arm}"
            arm_row = target_by_arm.get(chr_arm)
            out = {
                "contig": contig,
                "telogator_chr": tel_chr,
                "arm": arm,
                "chr_arm": chr_arm,
                "target_column": args.target_column,
                "has_target": "yes" if arm_row else "no",
                "target_length_bp": parse_int(arm_row.get(args.target_column)) if arm_row else "",
            }
            if not arm_row:
                missing_rows.append(out)
            target_rows.append(out)

    coord_fields = [
        "contig", "fai_length", "tele_bed_found",
        "5_telomere_start", "5_telomere_end", "5_tele_length",
        "3_telomere_start", "3_telomere_end", "3_tele_length",
        "left_non_telomeric_prefix_bp", "right_non_telomeric_suffix_bp",
        "status", "notes"
    ]
    write_tsv(args.out_prefix + ".coordinate_qc.tsv", coord_rows, coord_fields)

    target_fields = ["contig", "telogator_chr", "arm", "chr_arm", "target_column", "has_target", "target_length_bp"]
    write_tsv(args.out_prefix + ".target_coverage.tsv", target_rows, target_fields)
    write_tsv(args.out_prefix + ".missing_targets.tsv", missing_rows, target_fields)

    status_counter = Counter()
    for row in coord_rows:
        for s in row["status"].split(";"):
            status_counter[s] += 1

    print("Input QC summary")
    print(f"  contigs in fai:            {len(fai)}")
    print(f"  rows in tele-bed:          {len(tele)}")
    print(f"  target arms in arm_stats:  {len(target_by_arm)}")
    print(f"  expected arms from fai:    {len(fai) * 2}")
    print(f"  missing target arms:       {len(missing_rows)}")
    print("  coordinate status counts:")
    for key, val in sorted(status_counter.items()):
        print(f"    {key}: {val}")
    print(f"\nWrote: {args.out_prefix}.coordinate_qc.tsv")
    print(f"Wrote: {args.out_prefix}.target_coverage.tsv")
    print(f"Wrote: {args.out_prefix}.missing_targets.tsv")


if __name__ == "__main__":
    main()
