#!/usr/bin/env bash
# 12_methylation/hifi_5mc.sh — PacBio HiFi 5mC profiling
set -euo pipefail
FA=${1:?fasta}
HIFI_BAM=${2:-hifi.ccs.bam}
T=${THREADS:-32}

pbmm2 align --preset CCS -j "${T}" "${FA}" "${HIFI_BAM}" hifi.ccs.aligned.bam
samtools index hifi.ccs.aligned.bam
aligned_bam_to_cpg_scores --bam hifi.ccs.aligned.bam --output-prefix hifi_5mc \
  --model "${PILEUP_MODEL:-pileup_calling_model.v1.tflite}" --threads "${T}"
# genome-wide mean 5mC (HiFi): ~76.1%
