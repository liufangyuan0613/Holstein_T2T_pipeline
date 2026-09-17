#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_telomere_correction.sh cow_t2t.fa cow_t2t.fa.fai tele.bed arm_stats.tsv
#
# Inputs:
#   1. genome FASTA
#   2. genome FASTA index (.fai)
#   3. current telomere BED generated from the same FASTA
#   4. telogator2 arm_stats.tsv

FA=${1:?input fasta required}
FAI=${2:?input fai required}
TELE_BED=${3:?tele.bed required}
ARM_STATS=${4:?arm_stats.tsv required}

PREFIX="cow_t2t.telomere_corrected"

python 00_check_telomere_inputs.py \
  --fai "$FAI" \
  --tele-bed "$TELE_BED" \
  --arm-stats "$ARM_STATS" \
  --out-prefix "${PREFIX}.input_qc" \
  --target-column plot_value_bp

python 01_make_telomere_correction_plan.py \
  --fai "$FAI" \
  --tele-bed "$TELE_BED" \
  --arm-stats "$ARM_STATS" \
  --out "${PREFIX}.plan.tsv" \
  --target-column plot_value_bp \
  --left-replace-mode terminal \
  --right-replace-mode terminal

python 02_correct_telomeres.py \
  --input-fa "$FA" \
  --plan "${PREFIX}.plan.tsv" \
  --output-fa "${PREFIX}.fa" \
  --summary "${PREFIX}.apply_summary.tsv"

samtools faidx "${PREFIX}.fa"

python 03_evaluate_telomere_correction.py \
  --corrected-fa "${PREFIX}.fa" \
  --plan "${PREFIX}.plan.tsv" \
  --out "${PREFIX}.exact_motif_evaluation.tsv"

cat <<EOF

Primary correction is complete.

Recommended external evaluation:
  1. Re-run the same telomere scanner that produced the original tele.bed:
       your_telomere_scanner ${PREFIX}.fa > ${PREFIX}.tele.bed

  2. Compare scanner-detected lengths with telogator2 targets:
       python 04_compare_corrected_telebed_to_targets.py \\
         --plan ${PREFIX}.plan.tsv \\
         --corrected-tele-bed ${PREFIX}.tele.bed \\
         --out ${PREFIX}.telebed_vs_telogator_targets.tsv

EOF
