#!/usr/bin/env bash
# 06_quality_assessment/merqury_qv.sh — reference-free consensus quality
set -euo pipefail
FA=${1:?fasta required}
SHORT1=${2:-short_R1.clean.fq.gz}
SHORT2=${3:-short_R2.clean.fq.gz}
HIFI=${4:-hifi.fastq.gz}

# 21-mer databases from short reads and HiFi, merfin-filtered (k=21)
meryl k=21 count output short.k21.meryl <(cat "${SHORT1}" "${SHORT2}")
meryl k=21 count output hifi.k21.meryl "${HIFI}"
merfin -filter -k 21 --db short.k21.meryl --db hifi.k21.meryl --out k21.filtered
merqury k21.filtered.meryl "${FA}" merqury_out
# Holstein_T2T final: QV 72.54
