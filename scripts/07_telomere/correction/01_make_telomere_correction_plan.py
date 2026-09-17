#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse

from telomere_utils import (
    read_fai, read_tsv, write_tsv, parse_int, parse_float,
    contig_to_telogator_chr
)


def make_plan(args):
    fai = read_fai(args.fai)
    tele = read_tsv(args.tele_bed)
    arms = read_tsv(args.arm_stats)

    tele_by_contig = {r["chr"]: r for r in tele}

    arm_by_id = {}
    for r in arms:
        arm_id = r.get("chr_arm") or ((r.get("chr") or "") + (r.get("arm") or ""))
        if arm_id:
            arm_by_id[arm_id] = r

    rows = []
    for contig, chrom_len in fai.items():
        tel_chr = contig_to_telogator_chr(contig, args.contig_prefix)
        tele_row = tele_by_contig.get(contig)

        for arm in ["p", "q"]:
            end_side = "left" if arm == "p" else "right"
            chr_arm = f"{tel_chr}{arm}"
            arm_row = arm_by_id.get(chr_arm)
            target = parse_int(arm_row.get(args.target_column)) if arm_row else None

            row = {
                "contig": contig,
                "telogator_chr": tel_chr,
                "arm": arm,
                "chr_arm": chr_arm,
                "end_side": end_side,
                "fai_length": chrom_len,
                "target_column": args.target_column,
                "target_telomere_length": target if target is not None else "",
                "motif": args.left_motif if end_side == "left" else args.right_motif,
                "replace_mode": args.left_replace_mode if end_side == "left" else args.right_replace_mode,
                "status": "OK",
                "notes": "",
            }

            if tele_row is None:
                row["status"] = "SKIP_NO_TELE_BED"
                row["notes"] = "contig missing from tele.bed"
                rows.append(row)
                continue

            if target is None:
                row["status"] = "SKIP_NO_TARGET"
                row["notes"] = f"missing {chr_arm} in arm_stats.tsv or target column"
                rows.append(row)
                continue

            # Current detected telomere interval from tele.bed
            if end_side == "left":
                tel_start = parse_int(tele_row["5_telomere_start"])
                tel_end = parse_int(tele_row["5_telomere_end"])
                current_len = parse_int(tele_row["5_tele_length"])
                mode = args.left_replace_mode
                if mode == "terminal":
                    repl_start, repl_end = 0, tel_end
                elif mode == "bed":
                    repl_start, repl_end = tel_start, tel_end
                else:
                    raise ValueError(f"Unsupported left replace mode: {mode}")
            else:
                tel_start = parse_int(tele_row["3_telomere_start"])
                tel_end = parse_int(tele_row["3_telomere_end"])
                current_len = parse_int(tele_row["3_tele_length"])
                mode = args.right_replace_mode
                if mode == "terminal":
                    repl_start, repl_end = tel_start, chrom_len
                elif mode == "bed":
                    repl_start, repl_end = tel_start, tel_end
                else:
                    raise ValueError(f"Unsupported right replace mode: {mode}")

            row.update({
                "current_telomere_start": tel_start,
                "current_telomere_end": tel_end,
                "current_telomere_length": current_len,
                "replacement_start": repl_start,
                "replacement_end": repl_end,
                "replaced_region_length": repl_end - repl_start,
                "delta_vs_detected_telomere_length": target - current_len,
                "expected_contig_length_delta": target - (repl_end - repl_start),
            })

            issues = []
            if not (0 <= tel_start <= tel_end <= chrom_len):
                issues.append("BAD_TELOMERE_COORD")
            if tel_end - tel_start != current_len:
                issues.append("CURRENT_LENGTH_MISMATCH")
            if not (0 <= repl_start <= repl_end <= chrom_len):
                issues.append("BAD_REPLACEMENT_COORD")
            if target < 0:
                issues.append("NEGATIVE_TARGET")
            if repl_end - repl_start < 0:
                issues.append("NEGATIVE_REPLACED_LENGTH")
            if end_side == "left" and mode == "terminal" and tel_start > 0:
                row["notes"] = f"terminal mode removes {tel_start} bp before detected 5-prime telomere"
            if end_side == "right" and mode == "terminal" and tel_end < chrom_len:
                suffix = chrom_len - tel_end
                row["notes"] = f"terminal mode removes {suffix} bp after detected 3-prime telomere"

            # Add telogator2 QC fields if present.
            if arm_row:
                for col in [
                    "n_alleles", "n_unique_anchor_alleles", "mean_tl_bp", "median_tl_bp",
                    "min_tl_bp", "max_tl_bp", "representative_allele_id",
                    "representative_tl_bp", "representative_anchor_position",
                    "representative_ref_samp", "representative_support_reads",
                    "representative_median_mapq", "has_ambiguous_allele", "plot_stat"
                ]:
                    row[col] = arm_row.get(col, "")

            if issues:
                row["status"] = "SKIP_" + ";".join(issues)
                row["notes"] = (row.get("notes", "") + "; " if row.get("notes") else "") + ";".join(issues)

            rows.append(row)

    return rows


def main():
    ap = argparse.ArgumentParser(description="Create a telomere correction plan from .fai, tele.bed, and arm_stats.tsv.")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--tele-bed", required=True)
    ap.add_argument("--arm-stats", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target-column", default="plot_value_bp",
                    help="Target length column in arm_stats.tsv. Default: plot_value_bp")
    ap.add_argument("--contig-prefix", default="cow_chr")
    ap.add_argument("--left-motif", default="CCCTAA",
                    help="Motif used to rewrite the left/p telomere.")
    ap.add_argument("--right-motif", default="TTAGGG",
                    help="Motif used to rewrite the right/q telomere.")
    ap.add_argument("--left-replace-mode", choices=["terminal", "bed"], default="terminal",
                    help="terminal: replace 0..5_telomere_end; bed: replace 5_telomere_start..5_telomere_end")
    ap.add_argument("--right-replace-mode", choices=["terminal", "bed"], default="terminal",
                    help="terminal: replace 3_telomere_start..contig_end; bed: replace 3_telomere_start..3_telomere_end")
    args = ap.parse_args()

    rows = make_plan(args)
    fields = [
        "contig", "telogator_chr", "arm", "chr_arm", "end_side", "fai_length",
        "current_telomere_start", "current_telomere_end", "current_telomere_length",
        "replacement_start", "replacement_end", "replaced_region_length",
        "target_column", "target_telomere_length",
        "delta_vs_detected_telomere_length", "expected_contig_length_delta",
        "motif", "replace_mode", "status", "notes",
        "n_alleles", "n_unique_anchor_alleles", "mean_tl_bp", "median_tl_bp",
        "min_tl_bp", "max_tl_bp", "representative_allele_id",
        "representative_tl_bp", "representative_anchor_position",
        "representative_ref_samp", "representative_support_reads",
        "representative_median_mapq", "has_ambiguous_allele", "plot_stat"
    ]
    write_tsv(args.out, rows, fields)

    ok = sum(1 for r in rows if r["status"] == "OK")
    skipped = len(rows) - ok
    print(f"Wrote correction plan: {args.out}")
    print(f"  rows:    {len(rows)}")
    print(f"  OK:      {ok}")
    print(f"  skipped: {skipped}")


if __name__ == "__main__":
    main()
