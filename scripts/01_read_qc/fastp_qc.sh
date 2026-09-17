#!/usr/bin/env bash
# 01_read_qc/fastp_qc.sh — read QC for the three data types (manuscript parameters)
set -euo pipefail

# --- MGI short-read PE150 (309.60 Gb raw) ---
fastp -i short_R1.fq.gz -I short_R2.fq.gz -o short_R1.clean.fq.gz -O short_R2.clean.fq.gz \
  -n 0 -f 5 -F 5 -t 5 -T 5 -q 20 \
  -h fastp_short.html -j fastp_short.json
fastqc short_R1.clean.fq.gz short_R2.clean.fq.gz

# --- Hi-C PE150 (407.22 Gb raw) ---
# fastp -z 4 (PE mode QC); pairs discarded when containing adapters, >5 N bases,
# or >40% bases with Q<15 in either read.
fastp -i hic_R1.fq.gz -I hic_R2.fq.gz -o hic_R1.clean.fq.gz -O hic_R2.clean.fq.gz \
  -z 4 -h fastp_hic.html -j fastp_hic.json

# --- PacBio HiFi: SMRT Link v25.1, default parameters (Q20 pass) ---
# CCS itself: min passes = 3, min predicted accuracy RQ = 0.99 (run on instrument SMRT Link).

# --- ONT ultra-long ---
# Base-call: dorado v0.7.1 sup mode, model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
dorado basecaller --emit-pod5 --emit-fastq \
  sup,dna_r10.4.1_e8.2_400bps_sup@v5.0.0 raw_pod5/ > ont.fastq
# keep pass reads only (sup mode, QV >= 7), then the >50 kb subset used for hybrid assembly:
seqkit seq --min-len 50000 ont.fastq > ont.50k.fastq
