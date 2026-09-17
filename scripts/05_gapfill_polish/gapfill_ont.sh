#!/usr/bin/env bash
# 05_gapfill_polish/gapfill_ont.sh — close residual gaps with ONT local assembly
# Result for Holstein A-1: all 6 gaps closed -> 3,012.77 Mb, gap-free.
set -euo pipefail

FA=${1:?scaffolded fasta required}
ONT=${2:-ont.fastq}
THREADS=${THREADS:-64}

# 1. locate gaps
seqkit locate -i -r -p "N+" "${FA}" > gaps.bed

# 2. align ONT reads; keep candidates that are unmapped, poorly aligned
#    (identity < 80% or aligned read coverage < 80%), or uniquely span both gap flanks
minimap2 --secondary=no -ax map-ont -t "${THREADS}" "${FA}" "${ONT}" > ont2gap.paf

# 3. per gap: collect candidate reads from flanking regions, local iterative assembly
#    (hifiasm on the read subset), align local contigs back with
#    minimap2 -ax asm5 --eqx, and accept a replacement only when the contig bridges
#    the gap with reliable overlap on BOTH flanks.
#    (implemented as a vendor/local utility; see docs/pipeline_overview.md#05)
python3 "$(dirname "$0")/gapfill_local_assembly.py" \
  --gaps gaps.bed --paf ont2gap.paf --reads "${ONT}" --genome "${FA}" \
  --out A-1.gapfilled.fa --threads "${THREADS}"

seqkit stats A-1.gapfilled.fa
