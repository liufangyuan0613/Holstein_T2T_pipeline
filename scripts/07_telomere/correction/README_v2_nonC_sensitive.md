# v2: non-C sensitive bin calling

这版优化的核心不是继续减小 bin，而是改变 bin 内状态判定：

- 原版：直接取 dominant state，容易被大量 C 吞掉
- v2：如果一个 bin 内 non-C 比例 >= --nonc-threshold，
      就优先输出该 bin 中最常见的 non-C state

推荐先试：

python extract_telomere_tvr_states_v2.py \
  --fasta genome.fa \
  --bed cow_t2t.telomere_corrected.tele.bed \
  -o assembly_telomere_v2 \
  --bin-bp 20 \
  --max-bp 15000 \
  --nonc-threshold 0.15 \
  --min-kmers-per-bin 1

如果你想让 non-C 更容易被保留，可试：
--nonc-threshold 0.10

然后在 RStudio 里打开：
plot_telomere_state_matrix_RStudio_v2.R

修改：
state_matrix_file <- "assembly_telomere_v2.state_matrix.tsv"
output_prefix <- "assembly_telomere_v2_TVR_map"

再直接 Source。
