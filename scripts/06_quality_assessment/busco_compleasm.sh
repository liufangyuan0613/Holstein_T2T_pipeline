#!/usr/bin/env bash
# 06_quality_assessment/busco_compleasm.sh — gene-space completeness (mammalia_odb10)
set -euo pipefail
FA=${1:?fasta required}
compleasm run -a "${FA}" -o compleasm_out -l mammalia_odb10 -t "${THREADS:-64}"
# Holstein_T2T final: 9,218/9,226 complete (99.92%; 98.8% single-copy, 1.12% duplicated),
# 0.04% fragmented, 0.04% missing.
