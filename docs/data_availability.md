# Data availability

## This study
| Resource | Accession |
|---|---|
| Raw sequencing data (PacBio HiFi, ONT ultra-long, MGI short-read, Hi-C), donor A-1 | NCBI BioProject PRJNA1413493 |
| Holstein_T2T genome assembly | deposited under PRJNA1413493 (GenBank/CNGB accession to be added on publication) |
| CENP-A CUT&Tag (cut_1, cut_2 + IgG controls) | SRR40514003, SRR40514004, SRR40514005, SRR40514006 |

## Re-used public / previously generated data
| Resource | Accession / note |
|---|---|
| Population WGS cohort (9,591 Holstein cows) | not publicly available; access via lead contact (Xihe Li, lixh@imu.edu.cn) |
| Rumen RNA-seq (317 samples) | BioProject PRJNA894354 |
| Public cattle CENP-A CUT&Tag (Pineda et al., 2025) | SRR28462358, SRR28462359 |
| CattleGTEx multi-tissue cis-eQTL/sQTL summary data | https://cgtex.roslin.ed.ac.uk/ |
| FarmGTEx TWAS prediction models (ARS-UCD1.2/bosTau9) | https://twas.farmgtex.org/ |
| Reference genomes | ARS-UCD2.0, ARS-UCD1.2, Wagyu UOA_Wagyu_1 |

## Code
This repository. Analysis-specific note: the telomere-correction suite in
`scripts/07_telomere/correction/` is self-contained Python (stdlib only) and runs
without installation.
