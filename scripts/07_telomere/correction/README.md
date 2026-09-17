# Telomere length correction scripts v2

This script set corrects terminal telomere lengths in a genome FASTA by combining:

1. `cow_t2t.fa.fai`: chromosome/contig lengths.
2. `tele.bed`: current terminal telomere positions and lengths in the FASTA.
3. `arm_stats.tsv`: telogator2-derived target telomere lengths.

## Assumptions

- `tele.bed` uses 0-based half-open intervals: `[start, end)`.
- `p` arm corresponds to the FASTA left/start side.
- `q` arm corresponds to the FASTA right/end side.
- The left/p terminal telomere motif is rewritten as `CCCTAA`.
- The right/q terminal telomere motif is rewritten as `TTAGGG`.
- The default target length is `arm_stats.tsv:plot_value_bp`.

## Important correction modes

The default plan uses `terminal` mode:

- left/p: replace `[0, 5_telomere_end)` with `CCCTAA` repeated to the target length.
- right/q: replace `[3_telomere_start, contig_length)` with `TTAGGG` repeated to the target length.

This removes short non-telomeric sequence outside the detected telomeric block at the chromosome ends. This is usually the desired behavior for terminal telomere correction.

If you only want to replace the detected telomeric repeat block and preserve sequence outside it, use `bed` mode:

```bash
python 01_make_telomere_correction_plan.py \
  --fai cow_t2t.fa.fai \
  --tele-bed tele.bed \
  --arm-stats arm_stats.tsv \
  --out cow_t2t.bed_mode.plan.tsv \
  --left-replace-mode bed \
  --right-replace-mode bed
```

## Standard workflow

```bash
bash run_telomere_correction.sh cow_t2t.fa cow_t2t.fa.fai tele.bed arm_stats.tsv
```

The main outputs are:

- `cow_t2t.telomere_corrected.plan.tsv`
- `cow_t2t.telomere_corrected.fa`
- `cow_t2t.telomere_corrected.apply_summary.tsv`
- `cow_t2t.telomere_corrected.exact_motif_evaluation.tsv`

## Post-correction evaluation

After correction, re-run the same tool that generated the original `tele.bed`:

```bash
your_telomere_scanner cow_t2t.telomere_corrected.fa > cow_t2t.telomere_corrected.tele.bed
```

Then compare detected telomere lengths with telogator2 targets:

```bash
python 04_compare_corrected_telebed_to_targets.py \
  --plan cow_t2t.telomere_corrected.plan.tsv \
  --corrected-tele-bed cow_t2t.telomere_corrected.tele.bed \
  --out cow_t2t.telomere_corrected.telebed_vs_telogator_targets.tsv
```

## Notes on skipped arms

The script skips chromosome arms that are missing from `arm_stats.tsv`. In the current dataset, those arms are retained unchanged in the FASTA.
