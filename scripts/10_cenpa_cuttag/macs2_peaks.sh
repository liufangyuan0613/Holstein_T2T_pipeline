#!/usr/bin/env bash
# 10_cenpa_cuttag/macs2_peaks.sh — CENP-A peak calling
# Result: 6,091 / 3,930 peaks in the two replicates (q < 0.05).
set -euo pipefail
BAM=${1:?mapq10 bam}
OUT=${2:-cut_1.macs2}

macs2 callpeak -t "${BAM}" -f BAMPE -g 3.01e9 -n "${OUT}" --outdir macs2_out -q 0.05
