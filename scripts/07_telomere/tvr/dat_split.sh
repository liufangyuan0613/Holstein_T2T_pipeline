#!/bin/bash

# 检查输入参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.dat output.txt"
    exit 1
fi

# 输入的TRF .dat文件和输出文件
INPUT_DAT_FILE="$1"
OUTPUT_TXT_FILE="$2"

# 检查.dat文件是否存在
if [ ! -f "$INPUT_DAT_FILE" ]; then
    echo "File not found: $INPUT_DAT_FILE"
    exit 1
fi

# 打印列名
echo "Sequence ID\tStart Position\tEnd Position\tRepeat Size\tCopy Number\tPercent Match\tPercent Indels\tAlignment Score\tSequence" > "$OUTPUT_TXT_FILE"

# 使用awk提取信息并追加到输出文件
awk 'BEGIN { OFS = "\t" }
     /^Sequence: / { seq=$2 }
     /^[0-9]+/ {
       print seq, $1, $2, $3, $4, $5, $6, $7, $14
     }' "$INPUT_DAT_FILE" >> "$OUTPUT_TXT_FILE"
