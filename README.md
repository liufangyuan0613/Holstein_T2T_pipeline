# Holstein_T2T_pipeline

Analysis pipeline for the telomere-to-telomere (T2T) genome assembly of Holstein cattle
(*Bos taurus*) and its chromosome-end analyses, as described in:

> *A telomere-to-telomere genome of Holstein cattle resolves chromosome-end architecture
> and the genetic basis of telomere length* (manuscript; corresponding author: Xihe Li,
> Inner Mongolia University).

The final assembly spans **3.01 Gb in 30 gap-free chromosomes** (29 autosomes + X),
contig N50 **116.94 Mb**, Merqury QV **72.54**, BUSCO completeness **99.92 %**
(mammalia_odb10), and is annotated with **21,200 protein-coding genes**.

## Pipeline overview

| Stage | Directory | Key tools |
|---|---|---|
| 01 Read QC | `scripts/01_read_qc/` | fastp, FastQC, SMRT Link, Dorado, SeqKit |
| 02 Draft assemblies (×3) & selection | `scripts/02_draft_assembly/` | hifiasm, hifiasm --ul, Verkko |
| 03 Decontamination | `scripts/03_decontamination/` | RepeatMasker, BLASTN vs NCBI NT |
| 04 Hi-C scaffolding | `scripts/04_hic_scaffolding/` | BWA, bam-filter, HapHiC |
| 05 Gap filling & polishing | `scripts/05_gapfill_polish/` | minimap2, local assembly, yak, NextPolish2 |
| 06 Quality assessment | `scripts/06_quality_assessment/` | compleasm, Merqury/merfin, CRAQ, SAMtools, PanDepth |
| 07 Telomere scan & end calibration | `scripts/07_telomere/` | TIDK, SeqKit, Telogator2 + in-house correction suite |
| 08 Repeat annotation | `scripts/08_repeat_annotation/` | TRF, RepeatMasker, KMC, SRF |
| 09 Centromere prediction & comparison | `scripts/09_centromere/` | SRF, BLASTN, quarTeT CentroMiner, minimap2 |
| 10 CENP-A CUT&Tag | `scripts/10_cenpa_cuttag/` | Trimmomatic, BWA/bowtie2, MACS2, deepTools |
| 11 Gene annotation | `scripts/11_gene_annotation/` | STAR, StringTie, PASA, GeMoMa, AUGUSTUS, EvidenceModeler, TransposonPSI |
| 12 CpG methylation | `scripts/12_methylation/` | pbmm2, pb-CpG-tools, minimap2, modkit |
| 13 Collinearity | `scripts/13_collinearity/` | minimap2, SyRI, plotsr, NGenomeSyn |
| 14 Previously unresolved regions (PUR) | `scripts/14_pur/` | Winnowmap, BEDtools |
| 15 Subtelomere composition | `scripts/15_subtelomere/` | in-house GC/CpG window tools |

Stages 01–06 and 08, 11 reproduce the assembly/annotation workflow with the exact
parameters reported in the manuscript's STAR Methods. Stages 07, 09, 10, 13–15 contain
the in-house downstream analyses (telomere calibration, centromere comparison, CENP-A
CUT&Tag validation, collinearity, PUR and subtelomere profiling), including the original
script sets used in the study (kept verbatim under `scripts/07_telomere/correction/`,
`scripts/07_telomere/tvr/`, `scripts/09_centromere/`, `scripts/14_pur/`,
`scripts/15_subtelomere/` and `scripts/tools/`).

## Quick start

```bash
# 1. QC of the three data types
bash scripts/01_read_qc/fastp_qc.sh

# 2. Draft assemblies (run all three, then select; see scripts/02_draft_assembly/README.md)
bash scripts/02_draft_assembly/hifiasm.sh
bash scripts/02_draft_assembly/hifiasm_ul.sh
bash scripts/02_draft_assembly/verkko_mix.sh

# 3. Decontamination -> 4. Hi-C scaffolding -> 5. Gap fill & polish
bash scripts/03_decontamination/decontaminate.sh A-1.draft.fa
bash scripts/04_hic_scaffolding/hic_align_filter.sh A-1.decontam.fa
bash scripts/04_hic_scaffolding/haphic_scaffold.sh A-1.decontam.fa
bash scripts/05_gapfill_polish/gapfill_ont.sh A-1.haphic.fa
bash scripts/05_gapfill_polish/polish_nextpolish2.sh A-1.gapfilled.fa

# 6. Quality assessment
bash scripts/06_quality_assessment/busco_compleasm.sh A-1.polished.fa
bash scripts/06_quality_assessment/merqury_qv.sh    A-1.polished.fa

# 7. Telomere scan, read-based calibration and end correction
bash scripts/07_telomere/tidk_seqkit_scan.sh A-1.polished.fa
bash scripts/07_telomere/telogator2_read_calibration.sh
bash scripts/07_telomere/correction/run_telomere_correction.sh \
     cow_t2t.fa cow_t2t.fa.fai tele.bed arm_stats.tsv
```

Every script is self-documented in its header (inputs, outputs, parameters).

## Data availability

- Raw sequencing data (PacBio HiFi, ONT ultra-long, short-read, Hi-C) and the
  Holstein_T2T assembly: NCBI BioProject **PRJNA1413493**
- CENP-A CUT&Tag (this study): SRR40514003, SRR40514004, SRR40514005, SRR40514006
- Public cattle CENP-A CUT&Tag re-used here: SRR28462358, SRR28462359
- Rumen RNA-seq (317 samples): BioProject **PRJNA894354**

See `docs/data_availability.md` for details.

## Documentation

- `docs/pipeline_overview.md` — stage-by-stage inputs/outputs and logic
- `docs/software_versions.md` — pinned tool versions and citations
- `docs/data_availability.md` — accession inventory

## License

MIT (see `LICENSE`). If you use this pipeline or the telomere-correction suite,
please cite the manuscript above.
