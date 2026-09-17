#!/usr/bin/env bash
set -euo pipefail

# 用法:
#   bash keep_only_ChromosomeN.sh in.fa out.fa
#
# 规则:
#   - 只处理以 > 开头的header行
#   - 从header中提取 "Chromosome 1" / "Chromosome    12" 里的数字
#   - 输出header为：>Chromosome<数字>
#   - 如果某条header里找不到Chromosome数字，则原样保留（可按需改成报错）

IN_FA=${1:? "need input fasta"}
OUT_FA=${2:? "need output fasta"}

awk '
  /^>/{
    if (match($0, /Chromosome[[:space:]]*([0-9]+)/, a)) {
      print ">Chromosome" a[1]
    } else {
      # 找不到就原样输出；如需强制报错可改成：exit 1
      print $0
    }
    next
  }
  {print}
' "$IN_FA" > "$OUT_FA"
