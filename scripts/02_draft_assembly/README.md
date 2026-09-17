# Draft assembly comparison & selection

Three drafts (hifiasm / HiFiasmMix / VerkkoMix) are compared on:

1. **Contig statistics** — seqkit stats (size, count, N50). See `02_Assembly/04_compare`
   in the vendor delivery (`A-1.N50.tsv`).
2. **Collinearity dot plots** — minimap2 v2.26 + paf2dotplot.r, each draft vs ARS-UCD2.0
   and drafts against each other (refHifiasmul_qryHifiasm / refHifiasmul_qryVerkko /
   refHifiasm_qryVerkko).
3. **BUSCO completeness** — compleasm v0.2.5, mammalia_odb10 (`A-1.BUSCO.tsv`).
4. **Consensus quality** — Merqury v1.3 with merfin-filtered 21-mers from short reads
   and HiFi reads (`A-1.QV.tsv`).

Outcome for Holstein A-1:

| Draft | Size | Contigs | N50 | Verdict |
|---|---|---|---|---|
| hifiasm (HiFi only) | 3,315.78 Mb | 281 | 90.22 Mb | fewer errors, more fragmented |
| **HiFiasmMix** | **3,318.89 Mb** | **179** | **101.93 Mb** | **selected** |
| VerkkoMix | 6,046.55 Mb | 2,340 | 10.97 Mb | inflated size, fragmented |

HiFiasmMix gives the best balance of continuity and completeness and is carried into
`03_decontamination/`.
