#!/usr/bin/env bash

set -u
set -o pipefail

# =========================================================
# 用法:
#   bash extract_chr_ends_1Mb_trf_final.sh genome.fa out_prefix [window_bp]
#
# 示例:
#   bash extract_chr_ends_1Mb_trf_final.sh cow.fa cow_chr 1000000
#
# 主要输出:
#   cow_chr.5p_1Mb.fa
#   cow_chr.3p_1Mb.fa
#   cow_chr.5p_1Mb.fa.2.7.7.80.10.50.100.dat
#   cow_chr.3p_1Mb.fa.2.7.7.80.10.50.100.dat
#   cow_chr.5p_1Mb.trf.bed
#   cow_chr.3p_1Mb.trf.bed
#   cow_chr.5p_1Mb.telomere_like.bed
#   cow_chr.3p_1Mb.telomere_like.bed
#   cow_chr.1Mb.telomere_like.sorted.bed
#   cow_chr.1Mb.telomere_like.merged.bed
# =========================================================

FA=${1:? "ERROR: need genome fasta, e.g. genome.fa"}
OUT_PREFIX=${2:? "ERROR: need out prefix, e.g. chr_ends"}
WIN=${3:-1000000}

# -------------------------
# TRF 参数
# -------------------------
TRF_MATCH=2
TRF_MISMATCH=7
TRF_DELTA=7
TRF_PM=80
TRF_PI=10
TRF_MINSCORE=50
TRF_MAXPERIOD=100

# -------------------------
# 输出文件名
# -------------------------
OUT5="${OUT_PREFIX}.5p_1Mb.fa"
OUT3="${OUT_PREFIX}.3p_1Mb.fa"

DAT5="${OUT5}.${TRF_MATCH}.${TRF_MISMATCH}.${TRF_DELTA}.${TRF_PM}.${TRF_PI}.${TRF_MINSCORE}.${TRF_MAXPERIOD}.dat"
DAT3="${OUT3}.${TRF_MATCH}.${TRF_MISMATCH}.${TRF_DELTA}.${TRF_PM}.${TRF_PI}.${TRF_MINSCORE}.${TRF_MAXPERIOD}.dat"

BED5="${OUT_PREFIX}.5p_1Mb.trf.bed"
BED3="${OUT_PREFIX}.3p_1Mb.trf.bed"

FILTER5="${OUT_PREFIX}.5p_1Mb.telomere_like.bed"
FILTER3="${OUT_PREFIX}.3p_1Mb.telomere_like.bed"

MERGED_SORTED="${OUT_PREFIX}.1Mb.telomere_like.sorted.bed"
MERGED_FINAL="${OUT_PREFIX}.1Mb.telomere_like.merged.bed"

LOGDIR="${OUT_PREFIX}.trf_logs"
mkdir -p "$LOGDIR"

# -------------------------
# 打印参数
# -------------------------
echo "[INFO] Input fasta      : $FA"
echo "[INFO] Output prefix    : $OUT_PREFIX"
echo "[INFO] Window size      : $WIN"
echo "[INFO] Log directory    : $LOGDIR"

# -------------------------
# 检查依赖
# -------------------------
command -v samtools >/dev/null 2>&1 || {
    echo "[ERROR] samtools not found in PATH" >&2
    exit 1
}

command -v trf >/dev/null 2>&1 || {
    echo "[ERROR] trf not found in PATH" >&2
    exit 1
}

command -v bedtools >/dev/null 2>&1 || {
    echo "[ERROR] bedtools not found in PATH" >&2
    exit 1
}

# -------------------------
# 建立 faidx 索引
# -------------------------
if [[ ! -f "${FA}.fai" ]]; then
    echo "[INFO] Building FASTA index..."
    samtools faidx "$FA" || {
        echo "[ERROR] samtools faidx failed" >&2
        exit 1
    }
fi

# -------------------------
# 提取染色体两端 1 Mb
# -------------------------
extract_chr_ends() {
    echo "[INFO] Extracting chromosome ends..."

    : > "$OUT5"
    : > "$OUT3"

    while IFS=$'\t' read -r CHR LEN REST; do
        [[ -z "${CHR:-}" ]] && continue
        [[ -z "${LEN:-}" ]] && continue

        # 5' 端: 1 到 min(WIN, LEN)
        END5=$WIN
        if (( END5 > LEN )); then
            END5=$LEN
        fi
        REGION5="${CHR}:1-${END5}"

        # 3' 端: max(1, LEN-WIN+1) 到 LEN
        START3=$(( LEN - WIN + 1 ))
        if (( START3 < 1 )); then
            START3=1
        fi
        REGION3="${CHR}:${START3}-${LEN}"

        samtools faidx "$FA" "$REGION5" >> "$OUT5" || {
            echo "[ERROR] Failed extracting $REGION5" >&2
            exit 1
        }

        samtools faidx "$FA" "$REGION3" >> "$OUT3" || {
            echo "[ERROR] Failed extracting $REGION3" >&2
            exit 1
        }

    done < "${FA}.fai"

    [[ -s "$OUT5" ]] || { echo "[ERROR] $OUT5 is empty" >&2; exit 1; }
    [[ -s "$OUT3" ]] || { echo "[ERROR] $OUT3 is empty" >&2; exit 1; }

    echo "[INFO] Sequence extraction finished:"
    echo "       $OUT5"
    echo "       $OUT3"
}

# -------------------------
# 运行单个 TRF
# -------------------------
run_trf_one() {
    local fasta="$1"
    local base
    base=$(basename "$fasta")

    local stdout_log="${LOGDIR}/${base}.stdout.log"
    local stderr_log="${LOGDIR}/${base}.stderr.log"

    local dat_file="${fasta}.${TRF_MATCH}.${TRF_MISMATCH}.${TRF_DELTA}.${TRF_PM}.${TRF_PI}.${TRF_MINSCORE}.${TRF_MAXPERIOD}.dat"
    local mask_file="${fasta}.${TRF_MATCH}.${TRF_MISMATCH}.${TRF_DELTA}.${TRF_PM}.${TRF_PI}.${TRF_MINSCORE}.${TRF_MAXPERIOD}.mask"

    echo "[INFO] Running TRF on $fasta"
    echo "[INFO]   stdout -> $stdout_log"
    echo "[INFO]   stderr -> $stderr_log"

    rm -f "$dat_file" "$mask_file"

    trf "$fasta" \
        "$TRF_MATCH" "$TRF_MISMATCH" "$TRF_DELTA" \
        "$TRF_PM" "$TRF_PI" "$TRF_MINSCORE" "$TRF_MAXPERIOD" \
        -d -h -m >"$stdout_log" 2>"$stderr_log"

    local trf_exit=$?

    if [[ $trf_exit -ne 0 ]]; then
        echo "[WARN] TRF exited non-zero on $fasta (exit=$trf_exit)" >&2
    fi

    sync

    if [[ ! -s "$dat_file" ]]; then
        echo "[WARN] dat file not found after first run: $dat_file" >&2
        echo "[INFO] Re-running TRF once for $fasta" >&2

        trf "$fasta" \
            "$TRF_MATCH" "$TRF_MISMATCH" "$TRF_DELTA" \
            "$TRF_PM" "$TRF_PI" "$TRF_MINSCORE" "$TRF_MAXPERIOD" \
            -d -h -m >>"$stdout_log" 2>>"$stderr_log"

        trf_exit=$?
        sync
    fi

    if [[ -s "$dat_file" ]]; then
        echo "[INFO] TRF success: $dat_file"
    else
        echo "[ERROR] TRF failed to generate dat file for $fasta" >&2
        echo "[ERROR] Please check log: $stderr_log" >&2
        return 1
    fi

    if [[ -s "$mask_file" ]]; then
        echo "[INFO] MASK file generated: $mask_file"
    else
        echo "[WARN] MASK file missing: $mask_file" >&2
    fi

    return 0
}

# -------------------------
# dat 转 bed
# 修复点：
#   TRF 的坐标是相对于提取出来的子序列
#   这里解析 Sequence: chr:start-end
#   把局部坐标换算成染色体真实坐标
#
# 输出列:
# 1 chrom
# 2 start (1-based, real genomic coordinate)
# 3 end   (1-based, real genomic coordinate)
# 4 repeat_size
# 5 copy_number
# 6 percent_match
# 7 percent_indels
# 8 alignment_score
# 9 repeat_sequence
# -------------------------
dat_to_bed() {
    local input_dat="$1"
    local output_bed="$2"

    [[ -f "$input_dat" ]] || {
        echo "[ERROR] dat file not found: $input_dat" >&2
        return 1
    }

    echo "[INFO] Converting dat to real-coordinate bed: $input_dat -> $output_bed"

    awk '
    BEGIN {
        OFS = "\t"
        chrom = ""
        region_start = 1
        region_end = 1
    }

    /^Sequence: / {
        fullseq = $2

        chrom = fullseq
        sub(/:.*/, "", chrom)

        coord = fullseq
        sub(/^[^:]*:/, "", coord)

        split(coord, a, "-")
        region_start = a[1] + 0
        region_end   = a[2] + 0

        if (region_start < 1) region_start = 1
        if (region_end < region_start) region_end = region_start
        next
    }

    /^[0-9]+/ {
        local_start = $1 + 0
        local_end   = $2 + 0

        real_start = region_start + local_start - 1
        real_end   = region_start + local_end   - 1

        if (real_start < 1) real_start = 1
        if (real_end < real_start) real_end = real_start

        print chrom, real_start, real_end, $3, $4, $5, $6, $7, $14
    }' "$input_dat" > "$output_bed"

    if [[ -s "$output_bed" ]]; then
        echo "[INFO] BED generated: $output_bed"
    else
        echo "[WARN] BED file is empty: $output_bed" >&2
        : > "$output_bed"
    fi
}

# -------------------------
# 筛选端粒样 motif
# 条件:
#   第4列 repeat_size == 6
#   第5列 copy_number > 100
# 输出带表头
# -------------------------
filter_telomere_bed() {
    local input_bed="$1"
    local output_bed="$2"

    [[ -f "$input_bed" ]] || {
        echo "[ERROR] input bed not found: $input_bed" >&2
        return 1
    }

    echo "[INFO] Filtering telomere-like motifs: $input_bed -> $output_bed"
    echo "[INFO]   Criteria: repeat_size == 6 AND copy_number > 100"

    {
        echo -e "chrom\tstart\tend\trepeat_size\tcopy_number\tpercent_match\tpercent_indels\talignment_score\trepeat_sequence"
        awk 'BEGIN{OFS="\t"}
             $4 == 6 && $5 > 100 {
                 if ($2 < 1) $2 = 1
                 if ($3 < $2) $3 = $2
                 print $0
             }' "$input_bed"
    } > "$output_bed"

    if [[ $(wc -l < "$output_bed") -gt 1 ]]; then
        echo "[INFO] Filtered BED generated: $output_bed"
    else
        echo "[WARN] No records passed filter in: $input_bed" >&2
    fi
}

# -------------------------
# 合并筛选后的 bed
# sorted 输出保留全部列并添加表头
# merged 输出仅保留 3 列并添加表头
# -------------------------
merge_telomere_bed() {
    local input1="$1"
    local input2="$2"
    local sorted_out="$3"
    local merged_out="$4"

    local tmp_sorted="${sorted_out}.tmp"
    local tmp_merge_input="${merged_out}.merge_input.tmp"

    echo "[INFO] Sorting filtered BED files -> $sorted_out"

    {
        tail -n +2 "$input1" 2>/dev/null
        tail -n +2 "$input2" 2>/dev/null
    } | awk 'BEGIN{OFS="\t"} NF > 0 {
            if ($2 < 1) $2 = 1
            if ($3 < $2) $3 = $2
            print $0
        }' | sort -k1,1 -k2,2n > "$tmp_sorted"

    {
        echo -e "chrom\tstart\tend\trepeat_size\tcopy_number\tpercent_match\tpercent_indels\talignment_score\trepeat_sequence"
        cat "$tmp_sorted"
    } > "$sorted_out"

    echo "[INFO] Merging overlapping intervals -> $merged_out"

    awk 'BEGIN{OFS="\t"} NF > 0 { print $1, $2, $3 }' "$tmp_sorted" > "$tmp_merge_input"

    {
        echo -e "chrom\tstart\tend"
        bedtools merge -i "$tmp_merge_input"
    } > "$merged_out"

    rm -f "$tmp_sorted" "$tmp_merge_input"

    if [[ $(wc -l < "$merged_out") -gt 1 ]]; then
        echo "[INFO] Final merged BED generated: $merged_out"
    else
        echo "[WARN] Final merged BED is empty: $merged_out" >&2
    fi
}

# =========================================================
# 主流程
# =========================================================

extract_chr_ends

echo "[INFO] Starting TRF analyses..."
fail_n=0

run_trf_one "$OUT5" || fail_n=$((fail_n + 1))
run_trf_one "$OUT3" || fail_n=$((fail_n + 1))

echo "[INFO] TRF analyses finished."

if [[ $fail_n -ne 0 ]]; then
    echo "[ERROR] $fail_n TRF task(s) failed." >&2
    exit 1
fi

dat_to_bed "$DAT5" "$BED5" || exit 1
dat_to_bed "$DAT3" "$BED3" || exit 1

filter_telomere_bed "$BED5" "$FILTER5" || exit 1
filter_telomere_bed "$BED3" "$FILTER3" || exit 1

merge_telomere_bed "$FILTER5" "$FILTER3" "$MERGED_SORTED" "$MERGED_FINAL" || exit 1

echo "[INFO] All tasks completed successfully."
echo "[INFO] Output files:"
echo "       $OUT5"
echo "       $OUT3"
echo "       $DAT5"
echo "       $DAT3"
echo "       $BED5"
echo "       $BED3"
echo "       $FILTER5"
echo "       $FILTER3"
echo "       $MERGED_SORTED"
echo "       $MERGED_FINAL"
