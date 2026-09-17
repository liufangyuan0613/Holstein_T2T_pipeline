#!/usr/bin/env bash
# 02_draft_assembly/hifiasm_ul.sh — draft (ii): HiFiasmMix (HiFi + ONT >50 kb, --ul)
# Result for Holstein A-1: 3,318.89 Mb, 179 contigs, contig N50 101.93 Mb  [SELECTED]
set -euo pipefail

HIFI=${1:-hifi.fastq.gz}
ONT50K=${2:-ont.50k.fastq}
OUT=${3:-A-1.hifiasm_ul}

hifiasm -o "${OUT}" -t 64 --ul "${ONT50K}" "${HIFI}"
awk '/^S/{print ">"$2"\n"$3}' "${OUT}.bp.p_ctg.gfa" > "${OUT}.bp.p_ctg.fasta"
seqkit stats "${OUT}.bp.p_ctg.fasta" > "${OUT}.bp.p_ctg.fasta.n50"
