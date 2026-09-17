#!/usr/bin/env bash
set -euo pipefail

# 用法:
#   bash rename_fa_headers_to_cow_chr.sh in.fa out.fa
#
# 将类似 ">::cow_chr1:21589-13912593" 的header改为 ">cow_chr1"

IN_FA=${1:? "need input fasta"}
OUT_FA=${2:? "need output fasta"}

awk '
  BEGIN{OFS=""}
  /^>/{
    # 去掉开头的 > 以及紧随其后的任意冒号(:)
    line=$0
    sub(/^>:+/, "", line)

    # 去掉坐标部分 :start-end 以及其后所有内容
    sub(/:[0-9]+-[0-9]+.*/, "", line)

    # 如果还残留空格，只取第一段（可选，保险）
    split(line, a, /[[:space:]]+/)

    print ">" a[1]
    next
  }
  {print}
' "$IN_FA" > "$OUT_FA"
