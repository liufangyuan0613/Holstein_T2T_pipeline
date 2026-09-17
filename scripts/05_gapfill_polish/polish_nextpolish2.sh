#!/usr/bin/env bash
# 05_gapfill_polish/polish_nextpolish2.sh — repeat-aware polishing
# Result for Holstein A-1: 3,012.58 Mb, contig N50 116.94 Mb.
set -euo pipefail

FA=${1:?gapfilled fasta required}
SHORT1=${2:-short_R1.clean.fq.gz}
SHORT2=${3:-short_R2.clean.fq.gz}
HIFI=${4:-hifi.fastq.gz}
THREADS=${THREADS:-64}

# 1. two yak k-mer databases from quality-filtered short reads
yak count -b37 -k21 -o short.k21.yak <(cat "${SHORT1}" "${SHORT2}") <(cat "${SHORT1}" "${SHORT2}")
yak count -b37 -k31 -o short.k31.yak <(cat "${SHORT1}" "${SHORT2}") <(cat "${SHORT1}" "${SHORT2}")

# 2. HiFi alignments
minimap2 --secondary=no -ax map-hifi -t "${THREADS}" "${FA}" "${HIFI}" \
  | samtools sort -@ "${THREADS}" -o hifi2asm.bam
samtools index hifi2asm.bam

# 3. NextPolish2 default parameters (repeat-aware)
NextPolish2 -t "${THREADS}" -g "${FA}" -l hifi2asm.bam -k short.k21.yak -K short.k31.yak \
  -o A-1.polished.fa
seqkit stats A-1.polished.fa
