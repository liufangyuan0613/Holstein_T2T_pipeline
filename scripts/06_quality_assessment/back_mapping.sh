#!/usr/bin/env bash
# 06_quality_assessment/back_mapping.sh — read back-mapping & coverage
set -euo pipefail
FA=${1:?fasta required}
SHORT1=${2:-short_R1.clean.fq.gz}
SHORT2=${3:-short_R2.clean.fq.gz}
HIFI=${4:-hifi.fastq.gz}
ONT=${5:-ont.fastq}
T=${THREADS:-64}

minimap2 --secondary=no -ax sr       -t $T "${FA}" "${SHORT1}" "${SHORT2}" | samtools sort -@ $T -o short2asm.bam
minimap2 --secondary=no -ax map-hifi -t $T "${FA}" "${HIFI}"               | samtools sort -@ $T -o hifi2asm.bam
minimap2 --secondary=no -ax map-ont  -t $T "${FA}" "${ONT}"                | samtools sort -@ $T -o ont2asm.bam
for b in short2asm hifi2asm ont2asm; do
  samtools index $b.bam
  samtools stats $b.bam > $b.stats.txt
done
# Holstein_T2T: short 96.86%, HiFi 100%, ONT 93.60% (99.76% excluding ONT reads <1 kb)

# per-base coverage in 100-kb windows (PanDepth), plotted with karyoploteR
PanDepth -t $T -w 100000 -o short.100kb.depth short2asm.bam
PanDepth -t $T -w 100000 -o hifi.100kb.depth  hifi2asm.bam
PanDepth -t $T -w 100000 -o ont.100kb.depth   ont2asm.bam
