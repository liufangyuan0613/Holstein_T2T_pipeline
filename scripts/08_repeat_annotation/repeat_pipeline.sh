#!/usr/bin/env bash
# 08_repeat_annotation/repeat_pipeline.sh — integrated repeat annotation
# Result for Holstein_T2T: ~1.55 Gb (51.36%) repetitive sequence; LINEs 25.75%,
# SINEs 4.69%, LTRs 5.88%, satellite 9.42% of the genome.
set -euo pipefail
FA=${1:?fasta required}
T=${THREADS:-64}

# 1. tandem repeats (TRF)
trf "${FA}" 2 7 7 80 10 50 500 -f -d -h -l 6

# 2. interspersed repeats / TEs (RepeatMasker, custom library 'lib')
RepeatMasker -nolow -no_is -gff -norna -engine abblast -lib lib -pa "${T}" "${FA}"

# 3. satellite-rich regions: k-mer spectra + SRF unit reconstruction
kmc -k21 -t"${T}" -fm "${FA}" kmc_db tmp
SRF -l 6 kmc_db > srf_satellites.fa

# 4. repeat density in fixed windows along chromosomes (see
#    15_subtelomere/calc_gc_windows.py for the window pattern), integrated with
#    centromeric and subtelomeric annotations for downstream architecture analyses.
