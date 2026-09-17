#!/usr/bin/env bash
# 11_gene_annotation/abinitio_augustus.sh — ab initio prediction (AUGUSTUS v3.5.0)
set -euo pipefail
FA=${1:?softmasked fasta}
HINTS_GFF=${2:-pasa_train.3000.gff}

# 1. train AUGUSTUS on 3,000 high-confidence transcript-supported genes
#    (autoAugTrain pipeline; training set derived from PASA assemblies)
new_species.pl --species=holstein_t2t
etraining --species=holstein_t2t "${HINTS_GFF}"

# 2. run AUGUSTUS on the masked genome
augustus --species=holstein_t2t "${FA}" > augustus.holstein.gff
