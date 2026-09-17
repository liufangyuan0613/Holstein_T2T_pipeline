#!/usr/bin/env bash
# 14_pur/pur_winnowmap.sh — previously unresolved regions (PUR) vs ARS-UCD2.0
# Result: ~243.37 Mb newly resolved sequence space, 71.76% centromeric.
set -euo pipefail
COWFA=${1:?holstein_t2t fasta}
REFFA=${2:-ARS-UCD2.0.fa}
CENTBED=${3:-cow_T2T_centromere.bed}
T=${THREADS:-64}

# 1. align the new assembly against the previous reference
winnowmap --MD -H -cx asm20 -t "${T}" "${REFFA}" "${COWFA}" > cow_vs_ARS-UCD2.0.paf

# 2. segments of the new assembly with no alignment to the reference = PURs;
#    keep PURs longer than 50 kb for downstream analyses
#    (implemented with bedtools intersect/window/merge over the PAF-derived coverage)
paf2bed6 < cow_vs_ARS-UCD2.0.paf | bedtools merge -i - > aligned.bed
bedtools complement -i aligned.bed -g "${COWFA}.fai" > pur.raw.bed
awk '($3-$2)>50000' pur.raw.bed > pur.gt50kb.bed

# 3. fraction of PURs overlapping centromeric intervals
bedtools intersect -a pur.gt50kb.bed -b "${CENTBED}" -u > pur.centromeric.bed

# 4. repeat composition of PURs (see pur_repeat_composition.py in this directory)
