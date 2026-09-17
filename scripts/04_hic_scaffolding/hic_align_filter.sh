#!/usr/bin/env bash
# 04_hic_scaffolding/hic_align_filter.sh — Hi-C read alignment & filtering
set -euo pipefail

FA=${1:?decontaminated fasta required}
R1=${2:-hic_R1.clean.fq.gz}
R2=${3:-hic_R2.clean.fq.gz}
THREADS=${THREADS:-64}

bwa index "${FA}"
bwa mem -5SP -t "${THREADS}" "${FA}" "${R1}" "${R2}" \
  | samtools view -bS - > hic.bam

# retain pairs with mapping quality >= 1 and edit distance <= 3
bam-filter 1 --nm 3 hic.bam > hic.filtered.bam
samtools sort -@ "${THREADS}" -n -o hic.filtered.namesort.bam hic.filtered.bam
