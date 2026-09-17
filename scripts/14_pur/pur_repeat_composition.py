#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict, OrderedDict

# ---------------- basic interval utils (0-based half-open) ----------------
def merge_intervals(ivls):
    """Merge overlapping/adjacent intervals. ivls: list of (s,e)"""
    if not ivls:
        return []
    ivls = sorted(ivls)
    out = []
    cs, ce = ivls[0]
    for s, e in ivls[1:]:
        if s <= ce:  # overlap or adjacent
            if e > ce:
                ce = e
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out

def intersect_intervals(A, B):
    """Return intersection intervals of merged A and merged B."""
    if not A or not B:
        return []
    i = j = 0
    out = []
    while i < len(A) and j < len(B):
        a1, a2 = A[i]
        b1, b2 = B[j]
        s = max(a1, b1)
        e = min(a2, b2)
        if s < e:
            out.append((s, e))
        if a2 <= b2:
            i += 1
        else:
            j += 1
    return out

def subtract_intervals(A, B):
    """
    Return A - B. A and B must be merged & sorted.
    """
    if not A:
        return []
    if not B:
        return A[:]
    out = []
    j = 0
    for a1, a2 in A:
        cur = a1
        while j < len(B) and B[j][1] <= a1:
            j += 1
        k = j
        while k < len(B) and B[k][0] < a2:
            b1, b2 = B[k]
            if b1 > cur:
                out.append((cur, min(b1, a2)))
            cur = max(cur, b2)
            if cur >= a2:
                break
            k += 1
        if cur < a2:
            out.append((cur, a2))
    return out

def total_len(ivls):
    return sum(e - s for s, e in ivls)

# ---------------- parsers ----------------
def read_chrlen(path):
    # expects 2 columns with header: Chromosome Length (tab/space)
    chroms = OrderedDict()
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline()
        if not header:
            raise ValueError("Empty chr length file")
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chroms[parts[0]] = int(parts[1])
    return chroms

def read_bed(path):
    """Read BED: chr start end (0-based). return dict chr->list[(s,e)] merged later."""
    d = defaultdict(list)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            if e > s:
                d[chrom].append((s, e))
    return d

def read_repeats_gff(path):
    """
    Read RepeatMasker-like GFF:
    col1 chr, col4 start(1-based), col5 end(1-based), col3 type (e.g. LINE/Penelope, Satellite, DNA/Ginger)
    return dict chr->list[(s0,e0,major_class,full_type)]
    """
    d = defaultdict(list)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            ftype = parts[2]  # e.g. LINE/Penelope
            try:
                s1 = int(parts[3])
                e1 = int(parts[4])
            except ValueError:
                continue
            s0 = s1 - 1
            e0 = e1
            if e0 <= s0:
                continue
            major = ftype.split("/")[0] if "/" in ftype else ftype
            d[chrom].append((s0, e0, major, ftype))
    return d

# ---------------- classification ----------------
def classify_repeat_major(major):
    """
    Map major class -> one of: LINE, LTR, Satellite, SINE, Others_repeat
    """
    m = major.upper()
    if m == "LINE":
        return "LINE"
    if m == "LTR":
        return "LTR"
    if m == "SATELLITE":
        return "Satellite"
    if m == "SINE":
        return "SINE"
    return "Others_repeat"

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(
        description="Compute composition of repeat types within PUR regions (exclusive by precedence)."
    )
    ap.add_argument("--pur", required=True, help="PUR regions BED (0-based): chr start end")
    ap.add_argument("--repeats_gff", required=True, help="RepeatMasker repeats.gff (GFF with type in col3)")
    ap.add_argument("--centromere_bed", required=True, help="Centromere BED (0-based): chr start end")
    ap.add_argument("--chrlen", required=True, help="Chromosome length file (2 cols: Chromosome Length)")
    ap.add_argument("--rdna_bed", default=None, help="Optional rDNA BED (0-based). If not given, rDNA=0.")
    ap.add_argument("-o", "--out", required=True, help="Output TSV summary")
    ap.add_argument("--per_chr_out", default=None, help="Optional per-chromosome TSV output")
    args = ap.parse_args()

    chroms = read_chrlen(args.chrlen)

    pur = read_bed(args.pur)
    cen = read_bed(args.centromere_bed)
    rdna = read_bed(args.rdna_bed) if args.rdna_bed else defaultdict(list)
    reps = read_repeats_gff(args.repeats_gff)

    # merge PUR/centromere/rDNA per chr
    for c in list(pur.keys()):
        pur[c] = merge_intervals(pur[c])
    for c in list(cen.keys()):
        cen[c] = merge_intervals(cen[c])
    for c in list(rdna.keys()):
        rdna[c] = merge_intervals(rdna[c])

    categories = ["rDNA", "centromere", "LINE", "LTR", "Satellite", "SINE", "Others_repeat", "Unannotated_in_PUR"]
    overall_bp = dict((k, 0) for k in categories)
    overall_pur_bp = 0

    per_chr_rows = []

    # process chromosomes in chrlen order
    for chrom in chroms.keys():
        pur_iv = pur.get(chrom, [])
        if not pur_iv:
            continue

        pur_len = total_len(pur_iv)
        overall_pur_bp += pur_len

        # start with remaining PUR
        rem = pur_iv

        # --- rDNA ---
        rdna_iv = rdna.get(chrom, [])
        hit = merge_intervals(intersect_intervals(rem, rdna_iv))
        rdna_bp = total_len(hit)
        rem = merge_intervals(subtract_intervals(rem, hit))

        # --- centromere ---
        cen_iv = cen.get(chrom, [])
        hit = merge_intervals(intersect_intervals(rem, cen_iv))
        cen_bp = total_len(hit)
        rem = merge_intervals(subtract_intervals(rem, hit))

        # --- repeats by class (exclusive by precedence LINE > LTR > Satellite > SINE > Others_repeat) ---
        # collect repeat intervals per mapped class, but only those that overlap current rem to reduce cost
        rep_class_iv = defaultdict(list)
        for s0, e0, major, ftype in reps.get(chrom, []):
            cls = classify_repeat_major(major)
            rep_class_iv[cls].append((s0, e0))
        for cls in rep_class_iv:
            rep_class_iv[cls] = merge_intervals(rep_class_iv[cls])

        def take_cls(rem_iv, cls_iv):
            hit = merge_intervals(intersect_intervals(rem_iv, cls_iv))
            bp = total_len(hit)
            rem2 = merge_intervals(subtract_intervals(rem_iv, hit))
            return bp, rem2

        line_bp, rem = take_cls(rem, rep_class_iv.get("LINE", []))
        ltr_bp, rem = take_cls(rem, rep_class_iv.get("LTR", []))
        sat_bp, rem = take_cls(rem, rep_class_iv.get("Satellite", []))
        sine_bp, rem = take_cls(rem, rep_class_iv.get("SINE", []))
        oth_bp, rem = take_cls(rem, rep_class_iv.get("Others_repeat", []))

        un_bp = total_len(rem)  # remaining bases in PUR not assigned

        # accumulate overall
        overall_bp["rDNA"] += rdna_bp
        overall_bp["centromere"] += cen_bp
        overall_bp["LINE"] += line_bp
        overall_bp["LTR"] += ltr_bp
        overall_bp["Satellite"] += sat_bp
        overall_bp["SINE"] += sine_bp
        overall_bp["Others_repeat"] += oth_bp
        overall_bp["Unannotated_in_PUR"] += un_bp

        if args.per_chr_out:
            per_chr_rows.append((chrom, pur_len, rdna_bp, cen_bp, line_bp, ltr_bp, sat_bp, sine_bp, oth_bp, un_bp))

    # write overall summary
    with open(args.out, "w", encoding="utf-8") as out:
        out.write("Category\tbp\tMb\tpct_of_PUR\n")
        for cat in categories:
            bp = overall_bp[cat]
            mb = bp / 1e6
            pct = (bp / overall_pur_bp * 100.0) if overall_pur_bp > 0 else 0.0
            out.write(f"{cat}\t{bp}\t{mb:.6f}\t{pct:.4f}\n")
        out.write(f"PUR_total\t{overall_pur_bp}\t{overall_pur_bp/1e6:.6f}\t100.0000\n")

    # write per chromosome table if requested
    if args.per_chr_out:
        with open(args.per_chr_out, "w", encoding="utf-8") as out:
            out.write("chr\tPUR_bp\trDNA_bp\tcentromere_bp\tLINE_bp\tLTR_bp\tSatellite_bp\tSINE_bp\tOthers_repeat_bp\tUnannotated_bp\n")
            for r in per_chr_rows:
                out.write("\t".join(map(str, r)) + "\n")

if __name__ == "__main__":
    main()
