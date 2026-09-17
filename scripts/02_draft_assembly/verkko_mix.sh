#!/usr/bin/env bash
# 02_draft_assembly/verkko_mix.sh — draft (iii): VerkkoMix
# Result for Holstein A-1: 6,046.55 Mb, 2,340 contigs, contig N50 10.97 Mb
set -euo pipefail

HIFI=${1:-hifi.fastq.gz}
ONT50K=${2:-ont.50k.fastq}
OUTDIR=${3:-A-1.verkko}

verkko -d "${OUTDIR}" --hifi "${HIFI}" --nano "${ONT50K}" --no-cleanup
seqkit stats "${OUTDIR}/assembly.fasta" > "${OUTDIR}/assembly.fasta.n50"
