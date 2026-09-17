#!/usr/bin/env bash
# 09_centromere/compare_centromeres.sh — homologous centromere comparison
# (Holstein_T2T vs Wagyu UOA_Wagyu_1; manuscript Fig. S5–S7)
set -euo pipefail
COWFA=${1:?holstein t2t fasta}
REFFA=${2:-UOA_Wagyu_1.fa}
CENTBED=${3:-cow_T2T_centromere.bed}
T=${THREADS:-64}

# 1. extract homologous centromeric intervals and align (minimap2 asm5 --eqx)
mkdir -p cen_align
while read -r chr s e rest; do
  seqkit subseq --chr "${chr}" -r "${s}-${e}" "${COWFA}" > "cen_align/cow_${chr}.cen.fa"
  seqkit subseq --chr "${chr}" -r "1-999999999" "${REFFA}" > /dev/null 2>&1 || true
done < "${CENTBED}"

minimap2 -ax asm5 --eqx -t "${T}" "${COWFA}" "${REFFA}" > cow_vs_wagyu.paf
minimap2 -ax asm5 --eqx -t "${T}" "${REFFA}" "${COWFA}" > wagyu_vs_cow.paf

# 2. dotplots for homologous sequence blocks / repeat-array continuity
#    (paf2dotplot.r from dotPlotly, or moddotplot — see cow_trash/moddotplot.sh)

# 3. indel classification from alignment structure (insertion-dominant /
#    deletion-dominant / balanced) and higher-order-repeat organization summaries
#    (manuscript Fig. S7A–F); see 15_subtelomere/bed_window_base_gc.py helpers.
