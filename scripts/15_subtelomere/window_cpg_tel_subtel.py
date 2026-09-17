#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from bisect import bisect_left


def norm_chr(chrom: str) -> str:
    """Normalize chromosome naming so Chr01 -> cow_chr01."""
    chrom = chrom.strip()
    if chrom.startswith("cow_chr"):
        return chrom
    if chrom.startswith("Chr"):
        return "cow_chr" + chrom[3:]
    if chrom.startswith("chr"):
        return "cow_chr" + chrom[3:]
    return chrom


def read_cpg_bed(path):
    """
    Read CpG BED with columns:
    chr, start, end, coverage, methylated_count
    Returns dict[chrom] -> sorted list of (start, end, coverage, meth_count, meth_ratio)
    """
    by_chr = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom = norm_chr(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            cov = float(parts[3])
            meth = float(parts[4])
            ratio = meth / cov if cov > 0 else math.nan
            by_chr[chrom].append((start, end, cov, meth, ratio))
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: x[0])
    return by_chr


def read_telomere_bed(path):
    """
    Read wide telomere BED with header:
    chr 5_telomere_start 5_telomere_end 5_tele_length 3_telomere_start 3_telomere_end 3_tele_length
    Returns list of region dicts.
    """
    regions = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = norm_chr(row['chr'])
            s5 = int(row['5_telomere_start'])
            e5 = int(row['5_telomere_end'])
            s3 = int(row['3_telomere_start'])
            e3 = int(row['3_telomere_end'])
            regions.append({
                'chrom': chrom,
                'start': s5,
                'end': e5,
                'arm': '5p',
                'region_type': 'telomere',
                'region_name': f'{chrom}_5p_telomere'
            })
            regions.append({
                'chrom': chrom,
                'start': s3,
                'end': e3,
                'arm': '3p',
                'region_type': 'telomere',
                'region_name': f'{chrom}_3p_telomere'
            })
    return regions


def infer_arm(region_name: str, start: int, end: int) -> str:
    if '_5p_' in region_name or region_name.endswith('_5p') or '_5p' in region_name:
        return '5p'
    if '_3p_' in region_name or region_name.endswith('_3p') or '_3p' in region_name:
        return '3p'
    # fallback not ideal, but keep explicit
    raise ValueError(f'Cannot infer arm (5p/3p) from region name: {region_name}')


def read_subtel_bed(path):
    """Read standard 4-col BED: chrom start end name."""
    regions = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                raise ValueError('subtel BED must have at least 4 columns: chrom start end region_name')
            chrom = norm_chr(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            arm = infer_arm(name, start, end)
            regions.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'arm': arm,
                'region_type': 'subtelomere',
                'region_name': name
            })
    return regions


def build_start_index(records):
    starts = [r[0] for r in records]
    return starts


def window_stats(records, starts, region, window_size):
    chrom = region['chrom']
    rstart = region['start']
    rend = region['end']
    arm = region['arm']
    rtype = region['region_type']
    rname = region['region_name']

    out = []
    region_len = rend - rstart
    if region_len <= 0:
        return out

    # For forward genomic iteration. Orientation is added separately.
    n_windows = math.ceil(region_len / window_size)

    # Candidate records start before region end.
    left = bisect_left(starts, rstart - 1)
    cand = []
    for rec in records[left:]:
        if rec[0] >= rend:
            break
        # keep if overlaps region at all
        if rec[1] > rstart and rec[0] < rend:
            cand.append(rec)

    # CpGs are single-base intervals here, so start-based inclusion is effectively enough.
    for i in range(n_windows):
        wstart = rstart + i * window_size
        wend = min(rstart + (i + 1) * window_size, rend)
        vals = [x for x in cand if x[0] >= wstart and x[0] < wend]
        n_cpg = len(vals)
        cov_sum = sum(v[2] for v in vals)
        meth_sum = sum(v[3] for v in vals)
        per_site = [v[4] for v in vals if not math.isnan(v[4])]
        mean_site_ratio = sum(per_site) / len(per_site) if per_site else math.nan
        weighted_ratio = meth_sum / cov_sum if cov_sum > 0 else math.nan

        # Oriented coordinate: distance from chromosome end/telomere-adjacent side inward.
        if arm == '5p':
            dist_start = wstart - rstart
            dist_end = wend - rstart
            oriented_window_idx = i + 1
        else:  # 3p, reverse orientation so closest to chromosome end becomes window 1
            dist_start = rend - wend
            dist_end = rend - wstart
            oriented_window_idx = n_windows - i

        out.append({
            'region_type': rtype,
            'region_name': rname,
            'chrom': chrom,
            'arm': arm,
            'region_start': rstart,
            'region_end': rend,
            'region_length': region_len,
            'window_size': window_size,
            'window_index_genomic': i + 1,
            'window_index_from_terminal': oriented_window_idx,
            'window_start_genomic': wstart,
            'window_end_genomic': wend,
            'distance_from_terminal_start': dist_start,
            'distance_from_terminal_end': dist_end,
            'n_cpg_sites': n_cpg,
            'coverage_sum': round(cov_sum, 6),
            'methylated_sum': round(meth_sum, 6),
            'mean_cpg_level_unweighted': round(mean_site_ratio, 6) if not math.isnan(mean_site_ratio) else 'NA',
            'mean_cpg_level_weighted': round(weighted_ratio, 6) if not math.isnan(weighted_ratio) else 'NA'
        })
    return out


def write_tsv(rows, path, fieldnames):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize(rows):
    grp = defaultdict(list)
    for row in rows:
        key = (row['region_type'], row['arm'], row['window_index_from_terminal'], row['distance_from_terminal_start'], row['distance_from_terminal_end'])
        grp[key].append(row)
    out = []
    for key, vals in sorted(grp.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        def clean_numeric(col):
            arr = []
            for v in vals:
                x = v[col]
                if x != 'NA':
                    arr.append(float(x))
            return arr
        uw = clean_numeric('mean_cpg_level_unweighted')
        ww = clean_numeric('mean_cpg_level_weighted')
        nsites = [int(v['n_cpg_sites']) for v in vals]
        covs = [float(v['coverage_sum']) for v in vals]
        meths = [float(v['methylated_sum']) for v in vals]

        region_type, arm, win_idx, dstart, dend = key
        out.append({
            'region_type': region_type,
            'arm': arm,
            'window_index_from_terminal': win_idx,
            'distance_from_terminal_start': dstart,
            'distance_from_terminal_end': dend,
            'n_regions': len(vals),
            'mean_n_cpg_sites': round(sum(nsites) / len(nsites), 4),
            'sum_coverage': round(sum(covs), 6),
            'sum_methylated': round(sum(meths), 6),
            'mean_cpg_level_unweighted_mean': round(sum(uw) / len(uw), 6) if uw else 'NA',
            'mean_cpg_level_weighted_mean': round(sum(ww) / len(ww), 6) if ww else 'NA',
            'pooled_weighted_cpg_level': round(sum(meths) / sum(covs), 6) if sum(covs) > 0 else 'NA'
        })
    return out


def main():
    ap = argparse.ArgumentParser(
        description='Compute windowed CpG methylation levels across telomere and subtelomere regions.'
    )
    ap.add_argument('-c', '--cpg-bed', required=True, help='CpG BED: chr start end coverage methylated_count')
    ap.add_argument('-t', '--telomere-bed', required=True, help='Wide telomere BED/TSV from tele.bed')
    ap.add_argument('-s', '--subtelomere-bed', required=True, help='Subtelomere BED (4 columns)')
    ap.add_argument('-w', '--window-size', type=int, default=10000, help='Window size in bp [default: 10000]')
    ap.add_argument('-o', '--out-prefix', required=True, help='Output prefix')
    args = ap.parse_args()

    cpg = read_cpg_bed(args.cpg_bed)
    tel_regions = read_telomere_bed(args.telomere_bed)
    subtel_regions = read_subtel_bed(args.subtelomere_bed)
    regions = tel_regions + subtel_regions

    all_rows = []
    for region in regions:
        chrom = region['chrom']
        records = cpg.get(chrom, [])
        starts = build_start_index(records) if records else []
        all_rows.extend(window_stats(records, starts, region, args.window_size))

    detail_fields = [
        'region_type','region_name','chrom','arm','region_start','region_end','region_length',
        'window_size','window_index_genomic','window_index_from_terminal',
        'window_start_genomic','window_end_genomic',
        'distance_from_terminal_start','distance_from_terminal_end',
        'n_cpg_sites','coverage_sum','methylated_sum',
        'mean_cpg_level_unweighted','mean_cpg_level_weighted'
    ]
    write_tsv(all_rows, args.out_prefix + '.windowed_cpg.tsv', detail_fields)

    summary_rows = summarize(all_rows)
    summary_fields = [
        'region_type','arm','window_index_from_terminal','distance_from_terminal_start','distance_from_terminal_end',
        'n_regions','mean_n_cpg_sites','sum_coverage','sum_methylated',
        'mean_cpg_level_unweighted_mean','mean_cpg_level_weighted_mean','pooled_weighted_cpg_level'
    ]
    write_tsv(summary_rows, args.out_prefix + '.windowed_cpg.summary.tsv', summary_fields)

    print(f'Wrote {args.out_prefix}.windowed_cpg.tsv')
    print(f'Wrote {args.out_prefix}.windowed_cpg.summary.tsv')


if __name__ == '__main__':
    main()
