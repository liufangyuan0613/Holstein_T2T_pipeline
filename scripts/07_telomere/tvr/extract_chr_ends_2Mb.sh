#!/usr/bin/env bash
set -euo pipefail

# 用法:
#   bash extract_chr_ends_2Mb.sh genome.fa out_prefix [window_bp]
# 示例:
#   bash extract_chr_ends_2Mb.sh genome.fa chr_ends 2000000
#
# 输出:
#   out_prefix.5p_2Mb.fa   (所有chr的5'端窗口合并在一个fa)
#   out_prefix.3p_2Mb.fa   (所有chr的3'端窗口合并在一个fa)

FA=${1:? "ERROR: need genome fasta, e.g. genome.fa"}
OUT_PREFIX=${2:? "ERROR: need out prefix, e.g. chr_ends"}
WIN=${3:-2000000}

# 1) 检查 samtools
command -v samtools >/dev/null 2>&1 || {
  echo "ERROR: samtools not found in PATH" >&2
  exit 1
}

# 2) 建立 faidx 索引（若已存在则跳过）
if [[ ! -f "${FA}.fai" ]]; then
  samtools faidx "$FA"
fi

OUT5="${OUT_PREFIX}.5p_2Mb.fa"
OUT3="${OUT_PREFIX}.3p_2Mb.fa"
: > "$OUT5"
: > "$OUT3"

# 3) 遍历 .fai：第1列=序列名，第2列=长度
# 注意：不强制只选 "chr" 开头，默认对所有条目都提取
while IFS=$'\t' read -r CHR LEN _; do
  # 5' 端：1..min(WIN, LEN)
  END5=$WIN
  if (( END5 > LEN )); then END5=$LEN; fi
  REGION5="${CHR}:1-${END5}"

  # 3' 端：max(1, LEN-WIN+1)..LEN
  START3=$(( LEN - WIN + 1 ))
  if (( START3 < 1 )); then START3=1; fi
  REGION3="${CHR}:${START3}-${LEN}"

  # 提取并写入输出
  # samtools faidx 输出header形如：>chr:start-end
  samtools faidx "$FA" "$REGION5" >> "$OUT5"
  samtools faidx "$FA" "$REGION3" >> "$OUT3"

done < "${FA}.fai"

echo "DONE"
echo "  5' ends  -> $OUT5"
echo "  3' ends  -> $OUT3"
