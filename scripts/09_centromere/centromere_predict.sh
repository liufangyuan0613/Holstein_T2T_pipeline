#!/usr/bin/env bash
# 09_centromere/centromere_predict.sh — centromere coordinate prediction
# Final centromere coordinates of all 30 chromosomes = SRF/BLASTN-supported intervals
# integrated with quarTeT CentroMiner predictions (Supplementary Table S3).
set -euo pipefail
FA=${1:?fasta required}
TEGFF=${2:-repeats.TE.gff}
TRGFF=${3:-repeats.TR.gff}
T=${THREADS:-32}

# 1. high-repeat-density intervals -> KMC k-mer db -> SRF candidate satellite units
kmc -k21 -t"${T}" -fm "${FA}" kmc_db tmp
SRF -l 6 kmc_db > srf_units.fa

# 2. map candidate units back to the genome -> repeat-enriched candidate intervals
makeblastdb -in "${FA}" -dbtype nucl
blastn -query srf_units.fa -db "${FA}" -outfmt 6 -num_threads "${T}" -out srf_units.blastn.tsv

# 3. quarTeT CentroMiner with TRF + RepeatMasker annotations
quarTeT CentroMiner -t "${T}" -i "${FA}" --TE "${TEGFF}" --TR "${TRGFF}"

# 4. integrate SRF/BLASTN intervals + CentroMiner predictions into final coordinates
#    (study-specific curation; final table = Supplementary Table S3 in the manuscript).
