#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import re
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="按固定窗口统计FASTA每条染色体全长范围内窗口序列的最高ACTG碱基比例和GC比例"
    )
    parser.add_argument("-f", "--fasta", required=True, help="参考基因组FASTA文件（需有.fai索引）")
    parser.add_argument("-o", "--output", required=True, help="输出TSV文件")
    parser.add_argument("-w", "--window", type=int, default=100000, help="窗口大小，默认100000")
    parser.add_argument("--chrom", default=None,
                        help="只分析指定染色体（正则表达式），如 '^chr[0-9XY]+$' 或 '^1$' ")
    parser.add_argument("--exclude-chrom", default=None,
                        help="排除指定染色体（正则表达式），如 'chrM|MT|chrUn' ")
    parser.add_argument("--drop-last-short-window", action="store_true",
                        help="丢弃每条染色体末端不足window大小的短窗口（默认保留）")
    return parser.parse_args()

def safe_ratio(num, den):
    return 0.0 if den == 0 else num / den

def split_windows(start, end, window_size, drop_last_short=False):
    """将[start, end)切分为固定窗口；最后一个窗口可不足window_size（可选丢弃）"""
    s = start
    while s < end:
        e = min(s + window_size, end)
        if drop_last_short and (e - s) < window_size:
            break
        yield s, e
        s = e

def count_bases(seq):
    """统计ACTG，忽略N及其他字符"""
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    actg = a + c + g + t
    gc = g + c
    max_base_count = max(a, c, g, t) if actg > 0 else 0
    return a, c, g, t, actg, gc, max_base_count

def pick_max_base(a, c, g, t, actg):
    """若并列，按A/C/G/T顺序返回第一个"""
    if actg == 0:
        return "NA"
    base_counts = [("A", a), ("C", c), ("G", g), ("T", t)]
    base_counts.sort(key=lambda x: x[1], reverse=True)
    # 处理并列：按A/C/G/T优先
    maxv = base_counts[0][1]
    for b in ["A", "C", "G", "T"]:
        if {"A": a, "C": c, "G": g, "T": t}[b] == maxv:
            return b
    return base_counts[0][0]

def main():
    args = parse_args()

    if not os.path.exists(args.fasta):
        sys.exit(f"[ERROR] FASTA文件不存在: {args.fasta}")
    if not os.path.exists(args.fasta + ".fai"):
        sys.exit(f"[ERROR] 未找到FASTA索引文件: {args.fasta}.fai\n请先运行: samtools faidx {args.fasta}")

    chrom_re = re.compile(args.chrom) if args.chrom else None
    excl_re = re.compile(args.exclude_chrom) if args.exclude_chrom else None

    fa = pysam.FastaFile(args.fasta)

    with open(args.output, "w") as out:
        out.write("\t".join([
            "chrom",
            "chr_len",
            "window_start",
            "window_end",
            "window_len",
            "A", "C", "G", "T",
            "ACTG_count",
            "max_base",
            "max_base_count",
            "max_base_ratio",
            "GC_count",
            "GC_ratio"
        ]) + "\n")

        for chrom in fa.references:
            if chrom_re and (not chrom_re.search(chrom)):
                continue
            if excl_re and excl_re.search(chrom):
                continue

            chr_len = fa.get_reference_length(chrom)
            for wstart, wend in split_windows(
                0, chr_len, args.window, drop_last_short=args.drop_last_short_window
            ):
                try:
                    seq = fa.fetch(chrom, wstart, wend)
                except Exception as e:
                    print(f"[WARN] 提取序列失败 {chrom}:{wstart}-{wend}，跳过。错误: {e}", file=sys.stderr)
                    continue

                a, c, g, t, actg, gc, max_base_count = count_bases(seq)
                max_base = pick_max_base(a, c, g, t, actg)
                max_base_ratio = safe_ratio(max_base_count, actg)
                gc_ratio = safe_ratio(gc, actg)

                out.write("\t".join([
                    chrom,
                    str(chr_len),
                    str(wstart),
                    str(wend),
                    str(wend - wstart),
                    str(a), str(c), str(g), str(t),
                    str(actg),
                    max_base,
                    str(max_base_count),
                    f"{max_base_ratio:.6f}",
                    str(gc),
                    f"{gc_ratio:.6f}",
                ]) + "\n")

    fa.close()
    print(f"[INFO] 完成，结果已输出到: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
