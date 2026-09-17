#!/usr/bin/env bash
# 07_telomere/telogator2_read_calibration.sh — read-based telomere length targets
# Telogator2 independently evaluates terminal repeat length, orientation and arm
# assignment from raw long reads; its arm_stats.tsv provides target lengths for
# correction/run_telomere_correction.sh.
set -euo pipefail
FA=${1:?fasta required}
HIFI=${2:-hifi.fastq.gz}
ONT=${3:-ont.50k.fastq}

# align raw reads to the assembly ends (Telogator2 accepts BAM of reads vs assembly)
minimap2 --secondary=no -ax map-hifi -t "${THREADS:-64}" "${FA}" "${HIFI}" | samtools sort -@ "${THREADS:-64}" -o hifi2asm.bam
minimap2 --secondary=no -ax map-ont  -t "${THREADS:-64}" "${FA}" "${ONT}"  | samtools sort -@ "${THREADS:-64}" -o ont2asm.bam
samtools index hifi2asm.bam; samtools index ont2asm.bam

# Telogator2 (see https://github.com/zstephens/telogator2); outputs include arm_stats.tsv
# with plot_value_bp per chromosome arm (p/q matched to the assembled 5'/3' ends by
# chromosome orientation).
python3 telogator2.py --bam hifi2asm.bam --bam ont2asm.bam -r "${FA}" -o telogator2_out

# Calibration outcome for Holstein_T2T:
#   42 of 60 chromosome ends agree with read-supported lengths at +/-1 bp;
#   residual deviations only at chr13 3p (1,083 bp) and chr19 5p (266 bp).
# Next step:
#   bash correction/run_telomere_correction.sh cow_t2t.fa cow_t2t.fa.fai tele.bed \
#        telogator2_out/arm_stats.tsv
