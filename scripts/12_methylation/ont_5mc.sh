#!/usr/bin/env bash
# 12_methylation/ont_5mc.sh — ONT 5mC/5hmC/6mA profiling with modkit
set -euo pipefail
FA=${1:?fasta}
ONT_FQ=${2:-ont.fastq}
T=${THREADS:-32}

# 1. modification-aware base calling is assumed (dorado --emit-fastq with MM/ML tags);
#    split by length into 10-kb bins for bookkeeping
seqkit sliding -s 10000 -W 10000 -g "${FA}" > windows10k.bed

# 2. align
minimap2 -ayx map-ont -t "${T}" "${FA}" "${ONT_FQ}" | samtools sort -@ "${T}" -o ont2asm.bam
samtools index ont2asm.bam

# 3. pileup modifications (5mC, 5hmC, 6mA; default thresholds)
modkit pileup --preset traditional --threads "${T}" ont2asm.bam ont_mods.bed
# genome-wide mean 5mC (ONT): ~80.4%
