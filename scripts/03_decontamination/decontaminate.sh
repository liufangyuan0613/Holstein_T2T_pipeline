#!/usr/bin/env bash
# 03_decontamination/decontaminate.sh — contaminant screening before chromosome construction
# Usage: bash decontaminate.sh <draft.fa> [nt_db]
# Logic (manuscript): soft-mask -> split >1 Mb sequences into 50-kb bins -> BLASTN vs NT ->
# classify <=1 Mb contigs by best hit, longer contigs by majority of bins ->
# remove non-metazoan / mitochondrial / plastid contigs.
# Result for Holstein A-1: 3,256.47 Mb (contig N50 101.93 Mb), BUSCO 99.94%.
set -euo pipefail

FA=${1:?draft fasta required}
NT=${2:-/db/ncbi/nt/nt}
THREADS=${THREADS:-64}

# 1. soft-mask repeats (custom library 'lib')
RepeatMasker -nolow -no_is -gff -norna -engine abblast -lib lib -pa "${THREADS}" "${FA}"

# 2. split sequences >1 Mb into non-overlapping 50-kb bins for classification
seqkit split --by-size 1000000 --by-size-prefix part_ "${FA}"
for f in part_*.fa; do
  seqkit sliding -s 50000 -W 50000 "${f}" > "${f%.fa}.bins.fa"
done
cat part_*.bins.fa "${FA}" | seqkit sort -N > for_blast.fa

# 3. BLASTN vs NCBI NT (default parameters)
blastn -query for_blast.fa -db "${NT}" -outfmt 6 -num_threads "${THREADS}" \
  -out A-1_blast_out.summary

# 4. taxonomic assignment:
#    <=1 Mb contigs -> best hit; >1 Mb -> majority of bins (see vendor A-1_blast_out.summary.merge
#    with -identity 70.0 -coverage 1.0). Remove non-metazoan, mitochondrial, plastid contigs:
#    A-1_contig_contamination.txt lists the excluded ids.
seqkit grep -v -f A-1_contig_contamination.txt "${FA}" > A-1.decontam.fa
seqkit stats A-1.decontam.fa
compleasm run -a A-1.decontam.fa -o compleasm_decontam -l mammalia_odb10 -t "${THREADS}"
