#!/bin/bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome.fa> <telomere_table.tsv> <output_prefix>"
    echo "Example: bash $0 cow_t2t_v2.fa tele.bed cow_subtelomere"
    exit 1
fi

FA="$1"
TEL="$2"
PREFIX="$3"

OUTBED="${PREFIX}.bed"
OUTFA="${PREFIX}.fa"
GENOMESIZE="${PREFIX}.genome.size"

# 检查输入文件
if [ ! -f "$FA" ]; then
    echo "Error: fasta file not found -> $FA"
    exit 1
fi

if [ ! -f "$TEL" ]; then
    echo "Error: telomere table not found -> $TEL"
    exit 1
fi

# 检查依赖
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found"; exit 1; }

# 如果没有 fai，就建立索引
if [ ! -f "${FA}.fai" ]; then
    samtools faidx "$FA"
fi

# 生成染色体长度文件
cut -f1,2 "${FA}.fai" > "$GENOMESIZE"

# 根据 5' 和 3' 端坐标生成 subtelomere 500kb BED
awk -v OFS="\t" '
BEGIN {
    WIN=500000
}
NR==FNR {
    chrlen[$1]=$2
    next
}
NR==1 {
    next
}
{
    chr=$1
    t5s=$2
    t5e=$3
    t3s=$4
    t3e=$5

    gsub(/[[:space:]]/, "", chr)

    gsub(/[[:space:]]/, "", t5s)
    gsub(/[[:space:]]/, "", t5e)
    gsub(/[[:space:]]/, "", t3s)
    gsub(/[[:space:]]/, "", t3e)

    gsub(/"/, "", t5s)
    gsub(/"/, "", t5e)
    gsub(/"/, "", t3s)
    gsub(/"/, "", t3e)

    gsub(/,/, "", t5s)
    gsub(/,/, "", t5e)
    gsub(/,/, "", t3s)
    gsub(/,/, "", t3e)

    if (!(chr in chrlen)) {
        print "Warning: " chr " not found in genome size file" > "/dev/stderr"
        next
    }

    len = chrlen[chr]

    t5s += 0
    t5e += 0
    t3s += 0
    t3e += 0

    # 5p subtelomere: telomere end -> inward 500 kb
    sub5_start = t5e
    sub5_end   = t5e + WIN
    if (sub5_end > len) sub5_end = len

    if (sub5_start < sub5_end) {
        print chr, sub5_start, sub5_end, chr"_5p_subtelomere_500kb"
    }

    # 3p subtelomere: 500 kb before telomere start
    sub3_start = t3s - WIN
    if (sub3_start < 0) sub3_start = 0
    sub3_end = t3s

    if (sub3_start < sub3_end) {
        print chr, sub3_start, sub3_end, chr"_3p_subtelomere_500kb"
    }
}
' "$GENOMESIZE" "$TEL" > "$OUTBED"

# 提取序列
bedtools getfasta -fi "$FA" -bed "$OUTBED" -name -fo "$OUTFA"

echo "Done."
echo "Genome size file : $GENOMESIZE"
echo "Subtelomere BED  : $OUTBED"
echo "Subtelomere FASTA: $OUTFA"
