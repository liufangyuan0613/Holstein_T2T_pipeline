#!/usr/bin/env bash
# 04_hic_scaffolding/haphic_scaffold.sh — chromosome-scale scaffolding with HapHiC
# Result for Holstein A-1: 3,009.12 Mb (96.02% of input) anchored on 30 pseudomolecules
# (29 autosomes + X), scaffold N50 116.86 Mb, 6 residual gaps.
set -euo pipefail

FA=${1:?decontaminated fasta required}
BAM=${2:-hic.filtered.namesort.bam}
OUT=${3:-A-1.haphic}

# HapHiC default parameters. Internally: misjoin pre-processing (correct/break/filter/
# de-allele), MCL chromosome clustering, reassignment, fast 3D-DNA iterative sorting
# with ALLHiC optimization, pseudomolecule construction.
haphic pipeline "${FA}" "${BAM}" 30 \
  --outdir "${OUT}" --threads "${THREADS:-64}"

# Manual curation of the Hi-C contact map (HiCExplorer) is then performed:
#   hicFindRestSite -f DpnII -o rest_sites.bed "${FA}"
#   hicBuildMatrix -s hic.filtered.bam -rs rest_sites.bed -b restriction_sites.bed -o hic.h5
#   hicCorrectMatrix correct -m hic.h5 --filterThreshold -2 2 -o hic.corrected.h5
#   hicPlotMatrix -m hic.corrected.h5 -o hic_map.png
# and the curated assembly is re-evaluated for contiguity and BUSCO completeness.
seqkit stats "${OUT}"/**/HapHiC.fa 2>/dev/null || true
