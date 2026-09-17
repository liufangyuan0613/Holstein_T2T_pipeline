# Pipeline overview — Holstein_T2T

Inputs: PacBio HiFi (~71×), ONT ultra-long (~62×, N50 90 kb), MGI short-read PE150 (~103×),
Hi-C (~135×), all from peripheral blood of one seven-year-old Holstein cow (donor A-1).

## 01 Read QC (`scripts/01_read_qc/`)
- Short-read & Hi-C: `fastp` (short reads `-n 0 -f 5 -F 5 -t 5 -T 5 -q 20`; Hi-C `-z 4` plus
  adapter/N/quality-pair rules), FastQC inspection.
- HiFi: SMRT Link v25.1 default QC (Q20).
- ONT: Dorado v0.7.1 sup base-calling; pass reads only; hybrid-assembly subset
  `seqkit seq --min-len 50000`.

## 02 Draft assemblies (`scripts/02_draft_assembly/`)
Three drafts are built and compared head-to-head (contig stats, collinearity dot plots vs
ARS-UCD2.0 and between drafts, compleasm BUSCO, Merqury QV):

| Draft | Size | Contigs | Contig N50 |
|---|---|---|---|
| hifiasm (HiFi only) | 3,315.78 Mb | 281 | 90.22 Mb |
| HiFiasmMix (HiFi + ONT>50 kb) | 3,318.89 Mb | 179 | 101.93 Mb |
| VerkkoMix | 6,046.55 Mb | 2,340 | 10.97 Mb |

HiFiasmMix was selected (best continuity/completeness balance; VerkkoMix is inflated and
fragmented). See `scripts/02_draft_assembly/README.md`.

## 03 Decontamination (`scripts/03_decontamination/`)
RepeatMasker soft-mask → split >1 Mb sequences into 50-kb bins → BLASTN vs NCBI NT →
classify ≤1 Mb contigs by best hit, longer ones by majority of bins → remove non-metazoan,
mitochondrial and plastid contigs. Result: 3,256.47 Mb, contig N50 101.93 Mb,
BUSCO 99.94 %.

## 04 Hi-C scaffolding (`scripts/04_hic_scaffolding/`)
`bwa mem -5SP` + `bam-filter 1 --nm 3` → HapHiC v1.0.5 default (misjoin correction,
MCL clustering, reassignment, fast-3D-DNA/ALLHiC ordering) → manual Hi-C map curation.
Result: 3,009.12 Mb (96.02 %) anchored on 30 pseudomolecules, scaffold N50 116.86 Mb,
6 residual gaps.

## 05 Gap filling & polishing (`scripts/05_gapfill_polish/`)
ONT reads uniquely spanning gap flanks (or unmapped/poorly aligned) → local iterative
assembly → gap replacement when both flanks are reliably bridged → all 6 gaps closed
(3,012.77 Mb, gap-free) → NextPolish2 repeat-aware polishing with two yak k-mer databases
(-b37 -k21 / -b37 -k31) + HiFi alignments → 3,012.58 Mb, contig N50 116.94 Mb.

## 06 Quality assessment (`scripts/06_quality_assessment/`)
- compleasm v0.2.5 (mammalia_odb10): 9,218/9,226 complete (99.92 %).
- Merqury v1.3 + merfin-filtered 21-mers (short-read + HiFi): QV 72.54.
- Back-mapping: short reads 96.86 %, HiFi 100 %, ONT 93.60 % (99.76 % after excluding
  <1 kb reads); per-base coverage in 100-kb windows (PanDepth).
- CRAQ: covered rate 99.53 %, low-confidence rate 1.00e-4, R-AQI 99.40, S-AQI 99.22.

## 07 Telomere (`scripts/07_telomere/`)
1. TIDK explore (motifs 4–12 bp) + SeqKit locate of TTAGGG/CCCTAA at chromosome termini
   (≥10 consecutive copies = telomere-positive end; 19/30 chromosomes after scaffolding).
2. Ends lacking complete arrays are extended by local HiFi assembly.
3. Telogator2 on raw HiFi + ONT reads provides arm-resolved target telomere lengths;
   the in-house suite `correction/` rewrites terminal arrays to read-supported lengths
   (p=CCCTAA on 5′/left, q=TTAGGG on 3′/right), then evaluates correction
   (00_check → 01_plan → 02_correct → samtools faidx → 03_evaluate → 04_compare).
   42 of 60 chromosome ends agree with read-based estimates at ±1 bp (residuals only at
   chr13-3p, 1,083 bp and chr19-5p, 266 bp).
4. `tvr/`: TVR state extraction, state-matrix plotting and multi-sample TVR analyses.

## 08 Repeat annotation (`scripts/08_repeat_annotation/`)
TRF `2 7 7 80 10 50 500 -f -d -h -l 6` + RepeatMasker (custom library,
`-nolow -no_is -gff -norna -engine abblast -lib`) + KMC k-mer spectra + SRF `-l 6`
satellite-unit reconstruction. Result: ~1.55 Gb (51.36 %) repetitive sequence.

## 09 Centromere (`scripts/09_centromere/`)
SRF/BLASTN-supported candidate intervals + quarTeT CentroMiner
(`-t 10 -i Genome --TE TEgff`) integrated into final centromere coordinates for all 30
chromosomes (Supplementary Table S3). Homologous comparison vs the Wagyu reference
(UOA_Wagyu_1): minimap2 `-ax asm5 --eqx` + indel classification + higher-order-repeat
organization summaries.

## 10 CENP-A CUT&Tag (`scripts/10_cenpa_cuttag/`)
Two biological replicates + matched IgG controls (Jiayin Bio, NovaSeq PE):
Trimmomatic → BWA (own data) / bowtie2 very-sensitive (public SRR28462358/9) → MAPQ ≥ 10
unique reads → MACS2 (q < 0.05; 6,091 / 3,930 peaks) → deepTools R/T × 10^6 tracks with
IgG as background. CENP-A signal validates centromere coordinates genome-wide.

## 11 Gene annotation (`scripts/11_gene_annotation/`)
Transcript evidence (STAR + StringTie → PASA + GeneMark-ST; 71,082 → 68,338 models),
homology (GeMoMa vs Bison, B. grunniens, B. indicus, B. javanicus, B. mutus; 29,740
models), ab initio (AUGUSTUS trained on 3,000 transcript-supported genes) →
EvidenceModeler integration → TransposonPSI filtering + frame/stop QC →
**21,200 protein-coding genes** (BUSCO-protein 98.28 %). Functional annotation:
NR/Swiss-Prot/KEGG/KOG/GO + InterProScan (96.03 % annotated).

## 12 CpG methylation (`scripts/12_methylation/`)
HiFi: pbmm2 `--preset CCS` + pb-CpG-tools. ONT: modification-aware FASTQ → minimap2
`-ayx map-ont` → modkit (5mC/5hmC/6mA). Genome-wide mean 5mC ≈ 76.1 % (HiFi) / 80.4 % (ONT).

## 13 Collinearity (`scripts/13_collinearity/`)
minimap2 `-ax asm5 --eqx` vs previous references → SyRI `-F S` → plotsr; chromosome-level
synteny with NGenomeSyn (see also `scripts/09_centromere/` for centromere-focused
alignments and dot plots).

## 14 PUR (`scripts/14_pur/`)
Winnowmap `-cx asm20 -H --MD` of Holstein_T2T vs ARS-UCD2.0 → BEDtools
(intersect/window/merge) → PURs > 50 kb (~243.37 Mb newly resolved sequence space;
71.76 % centromeric).

## 15 Subtelomere (`scripts/15_subtelomere/`)
500-kb subtelomeric intervals internal to terminal arrays; sliding-window GC content,
CpG density and methylation; chromosome-end feature matrix for clustering/ordination
(see manuscript Fig. 2G–I, Fig. S9/S10).
