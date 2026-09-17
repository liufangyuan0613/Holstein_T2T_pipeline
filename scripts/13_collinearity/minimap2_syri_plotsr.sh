#!/usr/bin/env bash
# 13_collinearity/minimap2_syri_plotsr.sh — chromosome-scale collinearity vs references
set -euo pipefail
COWFA=${1:?holstein_t2t fasta}
REFFA=${2:-ARS-UCD2.0.fa}
T=${THREADS:-64}

minimap2 --secondary=no -ax asm5 --eqx -t "${T}" "${REFFA}" "${COWFA}" > ref_vs_cow.paf
minimap2 --secondary=no -ax asm5 --eqx -t "${T}" "${COWFA}" "${REFFA}" > cow_vs_ref.paf
syri -c ref_vs_cow.sam -r "${REFFA}" -q "${COWFA}" -F S -k -s
plotsr --sr syri.out --genomes genomes.txt -o cow_vs_ref.synteny.png
# chromosome-level synteny overview with NGenomeSyn:
#   perl GetTwoGenomeSyn.pl -InGenomeA cow.fa -InGenomeB ref.fa -MappingBin minimap2 \
#        -OutPrefix cow_vs_ref -TopBinNum 30 -BinMapX 2000 -BinMapY 2000
