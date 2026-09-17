#!/usr/bin/env bash
# 10_cenpa_cuttag/deeptools_tracks.sh — signal tracks normalized to input coverage,
# scaled by R/T x 10^6 (R = reads in peaks, T = total reads), IgG as background.
set -euo pipefail
BAM=${1:?mapq10 bam}
IGG=${2:?igg mapq10 bam}
OUT=${3:-cut_1}

bamCompare -b1 "${BAM}" -b2 "${IGG}" --operation ratio --scaleFactorsMethod None \
  --normalizeUsing RPKM -o "${OUT}.vsIgG.bw" -p "${THREADS:-16}"
bamCoverage -b "${BAM}" --normalizeUsing RPKM -o "${OUT}.rpkm.bw" -p "${THREADS:-16}"
