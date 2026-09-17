#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
import statistics
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Calculate GC content in subtelomeric 500 kb regions using fixed windows "
            "(default 10 kb). Supports 5p/3p orientation relative to telomere."
        )
    )
    p.add_argument("-b", "--bed", required=True, help="Subtelomere BED file (4 columns recommended)")
    p.add_argument("-f", "--fasta", required=True, help="FASTA of subtelomere sequences extracted from BED")
    p.add_argument("-o", "--outprefix", required=True, help="Output prefix")
    p.add_argument("-w", "--window", type=int, default=10000, help="Window size, default: 10000")
    p.add_argument(
        "--partial",
        action="store_true",
        help="Keep last partial window if sequence length is not an exact multiple of window size",
    )
    return p.parse_args()


def read_fasta(path):
    seqs = {}
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks).upper()
    return seqs


def load_bed(path):
    records = []
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(f"BED line {i} has <3 columns: {line.strip()}")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            arm = infer_arm(name, start, end)
            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "name": name,
                    "arm": arm,
                    "coord_id": f"{chrom}:{start}-{end}",
                    "length": end - start,
                }
            )
    return records


def infer_arm(name, start, end):
    # Priority: infer from region name
    if re.search(r"(?:^|[_\-.])5p(?:[_\-.]|$)", name):
        return "5p"
    if re.search(r"(?:^|[_\-.])3p(?:[_\-.]|$)", name):
        return "3p"
    # Fallback: unknown orientation
    return "NA"


def normalize_header(h):
    return h.split("(")[0]


def match_sequence(record, fasta_dict):
    # Typical exact IDs
    candidates = [
        record["name"],
        record["coord_id"],
        normalize_header(record["name"]),
        normalize_header(record["coord_id"]),
    ]
    for key in candidates:
        if key in fasta_dict:
            return key, fasta_dict[key]

    # bedtools getfasta may create headers like chr:start-end or chr:start-end(+)
    for h, seq in fasta_dict.items():
        hn = normalize_header(h)
        if hn == record["coord_id"] or hn == record["name"]:
            return h, seq

    # relaxed matching: unique header containing both name or coordinates
    relaxed = []
    for h, seq in fasta_dict.items():
        hn = normalize_header(h)
        if record["name"] in hn or record["coord_id"] in hn:
            relaxed.append((h, seq))
    if len(relaxed) == 1:
        return relaxed[0]

    raise KeyError(
        f"Cannot match BED entry '{record['name']}' ({record['coord_id']}) to FASTA header. "
        f"Please ensure FASTA headers are region names or chrom:start-end."
    )


def gc_percent(seq):
    seq = seq.upper()
    atgc = sum(seq.count(x) for x in "ATGC")
    if atgc == 0:
        return None, 0
    gc = seq.count("G") + seq.count("C")
    return (gc / atgc) * 100.0, atgc


def region_windows(record, seq, window_size=10000, keep_partial=False):
    L = len(seq)
    full_n = L // window_size
    remainder = L % window_size
    n_windows = full_n + (1 if (keep_partial and remainder > 0) else 0)

    rows = []
    for i in range(n_windows):
        if record["arm"] == "3p":
            # For 3p, window 1 is the telomere-proximal LAST 10 kb of the genomic interval.
            sub_end = L - i * window_size
            sub_start = max(0, sub_end - window_size)
            window_seq = seq[sub_start:sub_end]
            g_start = record["start"] + sub_start
            g_end = record["start"] + sub_end
        else:
            # 5p or unknown: window 1 is from region start toward centromere.
            sub_start = i * window_size
            sub_end = min(L, sub_start + window_size)
            if (sub_end - sub_start) < window_size and not keep_partial:
                continue
            window_seq = seq[sub_start:sub_end]
            g_start = record["start"] + sub_start
            g_end = record["start"] + sub_end

        gc, valid_bases = gc_percent(window_seq)
        rows.append(
            {
                "region_name": record["name"],
                "chrom": record["chrom"],
                "arm": record["arm"],
                "region_start": record["start"],
                "region_end": record["end"],
                "region_length_bed": record["length"],
                "region_length_fasta": L,
                "window_id_from_telomere": i + 1,
                "window_size": len(window_seq),
                "distance_from_telomere_start": i * window_size,
                "distance_from_telomere_end": i * window_size + len(window_seq),
                "window_start_genomic": g_start,
                "window_end_genomic": g_end,
                "gc_percent": gc,
                "valid_bases": valid_bases,
            }
        )
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            rr = dict(r)
            if rr.get("gc_percent") is not None:
                rr["gc_percent"] = f"{rr['gc_percent']:.4f}"
            w.writerow(rr)


def summarize_by_distance(rows):
    bucket = defaultdict(list)
    for r in rows:
        if r["gc_percent"] is not None:
            bucket[(r["window_id_from_telomere"], r["distance_from_telomere_start"], r["distance_from_telomere_end"])].append(r["gc_percent"])
    summary = []
    for key in sorted(bucket):
        vals = bucket[key]
        wid, d1, d2 = key
        summary.append(
            {
                "window_id_from_telomere": wid,
                "distance_from_telomere_start": d1,
                "distance_from_telomere_end": d2,
                "n_regions": len(vals),
                "mean_gc_percent": round(statistics.mean(vals), 4),
                "median_gc_percent": round(statistics.median(vals), 4),
                "std_gc_percent": round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0,
                "min_gc_percent": round(min(vals), 4),
                "max_gc_percent": round(max(vals), 4),
            }
        )
    return summary


def try_plot(summary_rows, outprefix):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    if not summary_rows:
        return False

    x = [(r["distance_from_telomere_start"] + r["distance_from_telomere_end"]) / 2 / 1000 for r in summary_rows]
    y = [r["mean_gc_percent"] for r in summary_rows]
    yerr = [r["std_gc_percent"] for r in summary_rows]

    plt.figure(figsize=(10, 5))
    plt.plot(x, y, marker="o", linewidth=1)
    lower = [max(0, a - b) for a, b in zip(y, yerr)]
    upper = [a + b for a, b in zip(y, yerr)]
    plt.fill_between(x, lower, upper, alpha=0.2)
    plt.xlabel("Distance from telomere (kb)")
    plt.ylabel("Mean GC (%)")
    plt.title("GC profile across subtelomeric 500 kb regions")
    plt.tight_layout()
    plt.savefig(f"{outprefix}.gc_profile.png", dpi=300)
    plt.close()
    return True


def main():
    args = parse_args()

    bed_records = load_bed(args.bed)
    fasta_dict = read_fasta(args.fasta)

    all_rows = []
    matched = 0
    for rec in bed_records:
        header, seq = match_sequence(rec, fasta_dict)
        rec["matched_header"] = header
        rows = region_windows(rec, seq, window_size=args.window, keep_partial=args.partial)
        all_rows.extend(rows)
        matched += 1

    if not all_rows:
        raise SystemExit("No windows generated. Check FASTA/BED inputs.")

    detail_fields = [
        "region_name", "chrom", "arm", "region_start", "region_end",
        "region_length_bed", "region_length_fasta", "window_id_from_telomere",
        "window_size", "distance_from_telomere_start", "distance_from_telomere_end",
        "window_start_genomic", "window_end_genomic", "gc_percent", "valid_bases"
    ]
    write_tsv(f"{args.outprefix}.windows.tsv", all_rows, detail_fields)

    summary_rows = summarize_by_distance(all_rows)
    summary_fields = [
        "window_id_from_telomere", "distance_from_telomere_start", "distance_from_telomere_end",
        "n_regions", "mean_gc_percent", "median_gc_percent", "std_gc_percent",
        "min_gc_percent", "max_gc_percent"
    ]
    write_tsv(f"{args.outprefix}.summary.tsv", summary_rows, summary_fields)

    plotted = try_plot(summary_rows, args.outprefix)

    n_regions = len({r['region_name'] for r in all_rows})
    n_windows = len(all_rows)
    print(f"Matched regions: {matched}")
    print(f"Regions with windows: {n_regions}")
    print(f"Total windows: {n_windows}")
    print(f"Detailed output : {args.outprefix}.windows.tsv")
    print(f"Summary output  : {args.outprefix}.summary.tsv")
    if plotted:
        print(f"Plot            : {args.outprefix}.gc_profile.png")
    else:
        print("Plot            : skipped (matplotlib unavailable)")


if __name__ == "__main__":
    main()
