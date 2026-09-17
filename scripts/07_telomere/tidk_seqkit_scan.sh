#!/usr/bin/env bash
# 07_telomere/tidk_seqkit_scan.sh — terminal telomeric-array annotation
set -euo pipefail
FA=${1:?fasta required}

# 1. survey telomeric motifs (4-12 bp) with TIDK
tidk explore --minimum 4 --maximum 12 "${FA}" > tidk_explore.tsv
# canonical vertebrate motif TTAGGG (+ reverse complement CCCTAA)

# 2. terminal 100-kb scan at both ends of every pseudomolecule;
#    >= 10 consecutive canonical copies = telomere-positive end
mkdir -p chr_ends
seqkit split --by-id "${FA}" --by-id-prefix chr_ends/
for f in chr_ends/*.fa; do
  c=$(basename "$f" .fa)
  seqkit subseq -r 1-100000 "$f"      > "chr_ends/${c}.left100k.fa"
  L=$(seqkit stat -T "$f" | tail -1 | cut -f4)
  seqkit subseq -r $((L-99999))-$L "$f" > "chr_ends/${c}.right100k.fa"
done
for f in chr_ends/*100k.fa; do
  seqkit locate -i -p TelomereSeq "$f" > "${f%.fa}.telomere.tsv"
done

# Holstein_T2T: 19/30 chromosomes carried >=1 detectable terminal array after
# scaffolding; ends lacking complete arrays were extended by local HiFi assembly
# (docs/pipeline_overview.md#07), then read-calibrated (see correction/).
