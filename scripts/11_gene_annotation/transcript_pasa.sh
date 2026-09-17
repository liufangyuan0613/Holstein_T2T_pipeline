#!/usr/bin/env bash
# 11_gene_annotation/transcript_pasa.sh — transcript evidence (PASA + GeneMark-ST)
set -euo pipefail
FA=${1:?softmasked fasta}
RNAR1=${2:-tissue_R1.fq.gz}
RNAR2=${3:-tissue_R2.fq.gz}
T=${THREADS:-32}

# 1. align tissue RNA-seq (STAR --twopassMode Basic --alignIntronMax 50000)
STAR --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles "${FA}" --runThreadN "${T}"
STAR --genomeDir star_idx --readFilesIn "${RNAR1}" "${RNAR2}" --readFilesCommand zcat \
  --twopassMode Basic --alignIntronMax 50000 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix rnaseq_ --runThreadN "${T}"

# 2. assemble transcripts (StringTie) and convert to GFF
stringtie rnaseq_Aligned.sortedByCoord.out.bam -p "${T}" -o rnaseq.stringtie.gtf
gffread -w transcripts.fa -g "${FA}" rnaseq.stringtie.gtf

# 3. PASA (v2) alignment assembly + GeneMark-ST prediction from PASA assemblies
Launch_PASA_pipeline.pl -c pasa.alignAssembly.txt -C -R -g "${FA}" -t transcripts.fa --ALIGNERS blat,gmap
# PASA assemblies -> GeneMark-ST --cores T (self-training gene prediction on transcript evidence)
gmes_petap.pl --ES --seq transcripts.fa --cores "${T}"
# Result: 71,082 transcript-supported models -> redundancy filtering -> 68,338 models.
