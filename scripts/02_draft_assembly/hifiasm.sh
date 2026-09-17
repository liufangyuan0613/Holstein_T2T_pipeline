#!/usr/bin/env bash
# 02_draft_assembly/hifiasm.sh — draft (i): HiFi-only hifiasm
# Result for Holstein A-1: 3,315.78 Mb, 281 contigs, contig N50 90.22 Mb
set -euo pipefail

HIFI=${1:-hifi.fastq.gz}
OUT=${2:-A-1.hifiasm}

hifiasm -o "${OUT}" -t 64 "${HIFI}"
# primary contigs: ${OUT}.bp.p_ctg.gfa -> fasta
awk '/^S/{print ">"$2"\n"$3}' "${OUT}.bp.p_ctg.gfa" > "${OUT}.bp.p_ctg.fasta"
seqkit stats "${OUT}.bp.p_ctg.fasta" > "${OUT}.bp.p_ctg.fasta.n50"
