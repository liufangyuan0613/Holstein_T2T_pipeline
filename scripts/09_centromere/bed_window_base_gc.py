#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="按固定窗口统计BED区域序列的最高ACTG碱基比例和GC比例"
    )
    parser.add_argument("-b", "--bed", required=True, help="输入BED文件")
    parser.add_argument("-f", "--fasta", required=True, help="参考基因组FASTA文件（需有.fai索引）")
    parser.add_argument("-o", "--output", required=True, help="输出TSV文件")
    parser.add_argument("-w", "--window", type=int, default=100000, help="窗口大小，默认100000")
    parser.add_argument("--skip-header", action="store_true", help="BED文件无表头（默认按无表头处理，此参数仅兼容保留）")
    return parser.parse_args()

def safe_ratio(num, den):
    return 0.0 if den == 0 else num / den

def iter_bed_regions(bed_file):
    """
    读取BED前三列：chrom, start, end
    自动跳过空行、注释行、格式错误行
    """
    with open(bed_file, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                # 尝试按任意空白分隔
                fields = line.split()
                if len(fields) < 3:
                    print(f"[WARN] 跳过第{line_no}行（列数不足3）: {line}", file=sys.stderr)
                    continue

            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                # 如果第一行是表头等情况，跳过
                print(f"[WARN] 跳过第{line_no}行（start/end不是整数）: {line}", file=sys.stderr)
                continue

            if start < 0 or end <= start:
                print(f"[WARN] 跳过第{line_no}行（坐标非法）: {line}", file=sys.stderr)
                continue

            extra = fields[3:] if len(fields) > 3 else []
            yield line_no, chrom, start, end, extra

def split_windows(start, end, window_size):
    """将[start, end)切分为固定窗口，最后一个窗口可不足window_size"""
    s = start
    while s < end:
        e = min(s + window_size, end)
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

def main():
    args = parse_args()

    # 检查fasta和索引
    if not os.path.exists(args.fasta):
        sys.exit(f"[ERROR] FASTA文件不存在: {args.fasta}")
    if not os.path.exists(args.fasta + ".fai"):
        sys.exit(f"[ERROR] 未找到FASTA索引文件: {args.fasta}.fai\n请先运行: samtools faidx {args.fasta}")

    fa = pysam.FastaFile(args.fasta)
    ref_names = set(fa.references)

    with open(args.output, "w") as out:
        # 输出表头
        out.write("\t".join([
            "bed_line",
            "chrom",
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

        for line_no, chrom, region_start, region_end, extra in iter_bed_regions(args.bed):
            if chrom not in ref_names:
                print(f"[WARN] BED第{line_no}行染色体不在FASTA中，跳过: {chrom}", file=sys.stderr)
                continue

            for wstart, wend in split_windows(region_start, region_end, args.window):
                try:
                    seq = fa.fetch(chrom, wstart, wend)
                except Exception as e:
                    print(f"[WARN] 提取序列失败（第{line_no}行 {chrom}:{wstart}-{wend}），跳过。错误: {e}", file=sys.stderr)
                    continue

                a, c, g, t, actg, gc, max_base_count = count_bases(seq)

                # 找出最高的碱基（若并列，按A/C/G/T顺序返回第一个）
                base_counts = {"A": a, "C": c, "G": g, "T": t}
                max_base = max(base_counts, key=base_counts.get) if actg > 0 else "NA"

                max_base_ratio = safe_ratio(max_base_count, actg)
                gc_ratio = safe_ratio(gc, actg)

                out.write("\t".join([
                    str(line_no),
                    chrom,
                    str(wstart),
                    str(wend),
                    str(wend - wstart),
                    str(a), str(c), str(g), str(t),
                    str(actg),
                    max_base,
                    str(max_base_count),
                    f"{max_base_ratio:.6f}",
                    str(gc),
                    f"{gc_ratio:.6f}"
                ]) + "\n")

    fa.close()
    print(f"[INFO] 完成，结果已输出到: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
