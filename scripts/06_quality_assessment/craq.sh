#!/usr/bin/env bash
# 06_quality_assessment/craq.sh — single-base & structural accuracy (reference-free)
set -euo pipefail
FA=${1:?fasta required}
SHORTBAM=${2:-short2asm.bam}
LONGBAM=${3:-long2asm.bam}

# CRAQ default parameters; clipping-based regional errors (CRE -> R-AQI)
# and large structural errors (CSE -> S-AQI).
runCRAQ.sh -g "${FA}" -sr "${SHORTBAM}" -lr "${LONGBAM}" -m 500000 -x 3
# Holstein_T2T final: covered rate 99.53%, low-confidence rate 1.00e-4,
# R-AQI 99.40, S-AQI 99.22.
