#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from collections import defaultdict

from telomere_utils import (
    read_tsv, write_tsv, fasta_iter, write_fasta_record,
    parse_int, repeat_to_length
)


def main():
    ap = argparse.ArgumentParser(description="Rewrite FASTA terminal telomere lengths according to a correction plan.")
    ap.add_argument("--input-fa", required=True, help="Input genome FASTA.")
    ap.add_argument("--plan", required=True, help="telomere_correction_plan.tsv from 01_make_telomere_correction_plan.py")
    ap.add_argument("--output-fa", required=True, help="Corrected FASTA.")
    ap.add_argument("--summary", default="telomere_correction_apply_summary.tsv")
    ap.add_argument("--wrap", type=int, default=60)
    args = ap.parse_args()

    plan = read_tsv(args.plan)
    replacements_by_contig = defaultdict(list)
    for r in plan:
        if r.get("status") != "OK":
            continue
        contig = r["contig"]
        start = parse_int(r["replacement_start"])
        end = parse_int(r["replacement_end"])
        target_len = parse_int(r["target_telomere_length"])
        motif = r["motif"].upper()
        new_seq = repeat_to_length(motif, target_len)
        replacements_by_contig[contig].append({
            "arm": r["arm"],
            "end_side": r["end_side"],
            "start": start,
            "end": end,
            "target_len": target_len,
            "motif": motif,
            "new_seq": new_seq,
            "expected_delta": target_len - (end - start),
        })

    summary_rows = []
    seen = set()

    with open(args.output_fa, "w") as out:
        for name, desc, seq in fasta_iter(args.input_fa):
            old_len = len(seq)
            reps = replacements_by_contig.get(name, [])
            seen.add(name)

            status = "OK"
            notes = []
            # Check overlap in original coordinates.
            sorted_by_start = sorted(reps, key=lambda x: x["start"])
            for i in range(1, len(sorted_by_start)):
                if sorted_by_start[i]["start"] < sorted_by_start[i-1]["end"]:
                    status = "ERROR"
                    notes.append("overlapping replacement intervals")

            for rep in reps:
                if not (0 <= rep["start"] <= rep["end"] <= old_len):
                    status = "ERROR"
                    notes.append(
                        f"bad interval for {rep['arm']}: {rep['start']}-{rep['end']} outside 0-{old_len}"
                    )

            if status == "OK":
                # Apply from right to left, preserving original coordinates.
                new_seq = seq
                for rep in sorted(reps, key=lambda x: x["start"], reverse=True):
                    new_seq = new_seq[:rep["start"]] + rep["new_seq"] + new_seq[rep["end"]:]
                write_fasta_record(out, desc, new_seq, args.wrap)
                new_len = len(new_seq)
            else:
                # Do not modify the contig if intervals are unsafe.
                write_fasta_record(out, desc, seq, args.wrap)
                new_len = old_len

            summary_rows.append({
                "contig": name,
                "old_length": old_len,
                "new_length": new_len,
                "observed_length_delta": new_len - old_len,
                "expected_length_delta": sum(x["expected_delta"] for x in reps) if status == "OK" else "",
                "n_replacements": len(reps) if status == "OK" else 0,
                "arms_replaced": ",".join(x["arm"] for x in sorted_by_start) if status == "OK" else "",
                "status": status,
                "notes": "; ".join(notes),
            })

    missing_contigs = sorted(set(replacements_by_contig) - seen)
    for contig in missing_contigs:
        summary_rows.append({
            "contig": contig,
            "old_length": "",
            "new_length": "",
            "observed_length_delta": "",
            "expected_length_delta": "",
            "n_replacements": 0,
            "arms_replaced": "",
            "status": "ERROR_MISSING_IN_FASTA",
            "notes": "contig present in plan but missing from input FASTA",
        })

    fields = [
        "contig", "old_length", "new_length", "observed_length_delta",
        "expected_length_delta", "n_replacements", "arms_replaced", "status", "notes"
    ]
    write_tsv(args.summary, summary_rows, fields)

    ok = sum(1 for r in summary_rows if r["status"] == "OK")
    err = sum(1 for r in summary_rows if r["status"] != "OK")
    print(f"Wrote corrected FASTA: {args.output_fa}")
    print(f"Wrote apply summary:   {args.summary}")
    print(f"  contigs OK:    {ok}")
    print(f"  contigs error: {err}")


if __name__ == "__main__":
    main()
