#!/usr/bin/env bash
# 10_cenpa_cuttag/trimmomatic_bwa_align.sh — CUT&Tag processing (this study's own data)
set -euo pipefail
FA=${1:?holstein_t2t fasta}
R1=${2:?clean R1 fq}
R2=${3:?clean R2 fq}
SAMPLE=${4:-cut_1}
T=${THREADS:-32}

# 1. trim adapters (Trimmomatic default)
trimmomatic PE -threads "${T}" "${R1}" "${R2}" \
  "${SAMPLE}_1P.fq.gz" "${SAMPLE}_1U.fq.gz" "${SAMPLE}_2P.fq.gz" "${SAMPLE}_2U.fq.gz" \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:36

# 2. align with BWA (default parameters)
bwa index "${FA}"
bwa mem -t "${T}" "${FA}" "${SAMPLE}_1P.fq.gz" "${SAMPLE}_2P.fq.gz" | samtools view -bS - > "${SAMPLE}.bam"
samtools sort -@ "${T}" -o "${SAMPLE}.sorted.bam" "${SAMPLE}.bam"
samtools index "${SAMPLE}.sorted.bam"

# 3. unique, MAPQ >= 10 reads
samtools view -b -q 10 -F 260 "${SAMPLE}.sorted.bam" > "${SAMPLE}.mapq10.bam"
samtools index "${SAMPLE}.mapq10.bam"
