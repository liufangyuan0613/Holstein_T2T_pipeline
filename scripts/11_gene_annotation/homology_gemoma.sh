#!/usr/bin/env bash
# 11_gene_annotation/homology_gemoma.sh — homology-based gene prediction (GeMoMa)
# Reference protein sets: Bison, Bos grunniens, Bos indicus, Bos javanicus, Bos mutus.
# Result: 29,740 homology-supported models.
set -euo pipefail
FA=${1:?fasta}
T=${THREADS:-32}

mkdir -p gemoma
for REF in Bison.grunniens Bos.indicus Bos.javanicus Bos.mutus Bos.bison; do
  # extract CDS + intron-position info from each reference
  java -jar GeMoMa.jar CLI=Extractor \
      -a "${REF}.gff" -g "${REF}.fa" -out "${REF}.gemoma.cds"
  # map to the target genome
  java -jar GeMoMa.jar CLI=GeMoMaPipeline \
      -t "${FA}" -a "${REF}.gff" -g "${REF}.fa" \
      -outdir gemoma/${REF}_vs_target -threads "${T}"
done
# combine the per-reference predictions (GAF filter default)
java -jar GeMoMa.jar CLI=GAF -p gemoma/ -outdir gemoma/combined
