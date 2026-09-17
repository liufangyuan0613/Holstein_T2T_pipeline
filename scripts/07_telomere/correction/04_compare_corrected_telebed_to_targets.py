#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse

from telomere_utils import read_tsv, write_tsv, parse_int, classify_error


def main():
    ap = argparse.ArgumentParser(
        description="Compare a newly generated corrected.tele.bed against telogator2 target lengths in the correction plan."
    )
    ap.add_argument("--plan", required=True)
    ap.add_argument("--corrected-tele-bed", required=True)
    ap.add_argument("--out", default="corrected_telebed_vs_telogator_targets.tsv")
    ap.add_argument("--pass-bp", type=int, default=6)
    ap.add_argument("--warn-bp", type=int, default=100)
    args = ap.parse_args()

    plan = [r for r in read_tsv(args.plan) if r.get("status") == "OK"]
    tele = {r["chr"]: r for r in read_tsv(args.corrected_tele_bed)}

    rows = []
    for r in plan:
        contig = r["contig"]
        target = parse_int(r["target_telomere_length"])
        t = tele.get(contig)

        if t is None:
            obs = None
            status = "MISSING_CONTIG_IN_CORRECTED_TELEBED"
        else:
            obs = parse_int(t["5_tele_length"] if r["end_side"] == "left" else t["3_tele_length"])
            status = classify_error(abs(obs - target), args.pass_bp, args.warn_bp) if obs is not None else "MISSING_LENGTH"

        rows.append({
            "contig": contig,
            "arm": r["arm"],
            "chr_arm": r["chr_arm"],
            "end_side": r["end_side"],
            "target_telomere_length": target,
            "corrected_telebed_length": obs if obs is not None else "",
            "error_bp": (obs - target) if obs is not None else "",
            "abs_error_bp": abs(obs - target) if obs is not None else "",
            "status": status,
        })

    fields = [
        "contig", "arm", "chr_arm", "end_side", "target_telomere_length",
        "corrected_telebed_length", "error_bp", "abs_error_bp", "status"
    ]
    write_tsv(args.out, rows, fields)

    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    print(f"Wrote comparison: {args.out}")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
