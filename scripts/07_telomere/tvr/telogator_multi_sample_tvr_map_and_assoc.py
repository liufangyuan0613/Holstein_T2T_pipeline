#!/usr/bin/env python3
import argparse
import collections
import glob
import math
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from scipy.stats import spearmanr

PALETTE = {
    'C':'#2f36ff','V':'#5fd0ff','H':'#f3a7ff','F':'#49c34a','D':'#ff3b30',
    'L':'#f1ea2f','S':'#27d7d7','A':'#ff9f1a','T':'#b276ff','N':'#c9b4ff',
    'P':'#ff66c4','E':'#ff7f7f','M':'#a8ff60','Q':'#7dc65a','W':'#9acd32',
    'R':'#f28e2b','K':'#8c564b','G':'#98df8a','Y':'#c5b0d5','_':'#d9d9d9'
}

KEY_FEATURES = [
    'tvr_len', 'tvr_symbol_entropy', 'tvr_nonC_fraction',
    'tvr_transition_rate', 'support_n', 'mapq_median'
]


def chr_sort_key(arm_label):
    m = re.match(r'chr(\d+|X|Y)([pqU]?)$', str(arm_label))
    if not m:
        return (9999, 9, str(arm_label))
    chrom, arm = m.group(1), m.group(2)
    c = 1000 if chrom == 'X' else 1001 if chrom == 'Y' else int(chrom)
    a = {'p': 0, 'q': 1, 'U': 2, '': 3}.get(arm, 9)
    return (c, a, str(arm_label))


def derive_sample_id_from_path(fp: str) -> str:
    name = Path(fp).name
    if name == 'tlens_by_allele.tsv':
        return Path(fp).resolve().parent.name
    m = re.match(r'(.+?)_tlens_by_allele\.tsv$', name)
    if m:
        return m.group(1)
    return Path(fp).stem


def parse_list_field(val):
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s:
        return []
    return [x for x in s.split(',') if x != '']


def safe_float_list(val):
    out = []
    for x in parse_list_field(val):
        try:
            out.append(float(x))
        except Exception:
            pass
    return out


def load_manifest(args):
    rows = []
    if args.manifest:
        man = pd.read_csv(args.manifest, sep='\t')
        required = {'sample_id', 'tlens'}
        if not required.issubset(man.columns):
            raise ValueError('Manifest must contain columns: sample_id and tlens')
        for _, r in man.iterrows():
            rows.append({
                'sample_id': str(r['sample_id']),
                'tlens': str(r['tlens']),
                'final_zip': None if 'final_zip' not in man.columns or pd.isna(r.get('final_zip')) else str(r['final_zip'])
            })
        return rows

    if args.inputs:
        for fp in args.inputs:
            rows.append({'sample_id': derive_sample_id_from_path(fp), 'tlens': fp, 'final_zip': None})
        return rows

    if args.scan_root:
        fps = sorted(glob.glob(os.path.join(args.scan_root, '**', '*tlens_by_allele.tsv'), recursive=True))
        for fp in fps:
            rows.append({'sample_id': derive_sample_id_from_path(fp), 'tlens': fp, 'final_zip': None})
        return rows

    raise ValueError('Provide --manifest, -i/--inputs, or --scan-root')


def load_consensus(zip_fp):
    cons = {}
    if not zip_fp:
        return cons
    try:
        with zipfile.ZipFile(zip_fp) as zf:
            for name in zf.namelist():
                if name.startswith('fa/consensus_') and name.endswith('.fa'):
                    idx = int(re.search(r'(\d+)\.fa$', name).group(1))
                    seq = ''.join(zf.read(name).decode().splitlines()[1:])
                    cons[idx] = seq
    except Exception:
        return {}
    return cons


def shannon_entropy(seq: str) -> float:
    if not seq:
        return np.nan
    cnt = collections.Counter(seq)
    total = len(seq)
    return -sum((c/total) * math.log2(c/total) for c in cnt.values() if c > 0)


def transition_rate(seq: str) -> float:
    if not seq or len(seq) < 2:
        return np.nan
    transitions = sum(1 for i in range(1, len(seq)) if seq[i] != seq[i-1])
    return transitions / (len(seq) - 1)


def symbol_fraction(seq: str, sym: str) -> float:
    if not seq:
        return np.nan
    return seq.count(sym) / len(seq)


def build_allele_table(sample_id, tlens_fp, final_zip=None, min_support=0, min_mapq=0.0):
    df = pd.read_csv(tlens_fp, sep='\t')
    cons = load_consensus(final_zip)
    out_rows = []
    all_symbols = set(['C'])

    for _, r in df.iterrows():
        chr_raw = str(r.get('#chr', ''))
        primary_chr_arm = chr_raw.split(',')[0] if chr_raw else ''
        arm = primary_chr_arm[-1] if primary_chr_arm.endswith(('p', 'q')) else 'U'
        chrom = primary_chr_arm[:-1] if arm in ('p', 'q') else primary_chr_arm
        support_n = len(parse_list_field(r.get('supporting_reads', '')))
        mapqs = safe_float_list(r.get('read_mapq', ''))
        mapq_median = float(np.median(mapqs)) if mapqs else np.nan
        tvr_consensus = '' if pd.isna(r.get('tvr_consensus')) else str(r.get('tvr_consensus'))
        all_symbols.update(list(tvr_consensus))
        full_consensus = cons.get(int(r.get('allele_id')) if not pd.isna(r.get('allele_id')) else -1, '')
        seq_source = 'final_consensus' if full_consensus else 'tvr_consensus'

        row = {
            'sample_id': sample_id,
            'tlens_file': tlens_fp,
            'final_zip': final_zip or '',
            'allele_id': r.get('allele_id'),
            '#chr': chr_raw,
            'primary_chr_arm': primary_chr_arm,
            'chrom': chrom,
            'arm': arm,
            'position': r.get('position'),
            'ref_samp': r.get('ref_samp'),
            'is_ambiguous_chr': ',' in chr_raw,
            'TL_p75': pd.to_numeric(r.get('TL_p75'), errors='coerce'),
            'tvr_len': pd.to_numeric(r.get('tvr_len'), errors='coerce'),
            'support_n': support_n,
            'mapq_median': mapq_median,
            'tvr_consensus': tvr_consensus,
            'tvr_symbol_count': len(tvr_consensus),
            'tvr_symbol_entropy': shannon_entropy(tvr_consensus),
            'tvr_nonC_fraction': 1.0 - symbol_fraction(tvr_consensus, 'C') if tvr_consensus else np.nan,
            'tvr_transition_rate': transition_rate(tvr_consensus),
            'dominant_tvr_symbol': collections.Counter(tvr_consensus).most_common(1)[0][0] if tvr_consensus else '',
            'dominant_nonC_tvr_symbol': max({k:v for k,v in collections.Counter(tvr_consensus).items() if k != 'C'}.items(), key=lambda kv: kv[1])[0] if sum(v for k,v in collections.Counter(tvr_consensus).items() if k != 'C') > 0 else '',
            'full_consensus_len': len(full_consensus),
            'consensus_source': seq_source,
            'pass_min_support': support_n >= min_support,
            'pass_min_mapq': (not np.isnan(mapq_median)) and mapq_median >= min_mapq,
        }
        row['pass_basic_qc'] = bool(row['pass_min_support'] and row['pass_min_mapq'] and pd.notna(row['TL_p75']) and row['TL_p75'] > 0)
        if r.get('allele_id') is not None:
            row['allele_key'] = f"{sample_id}|{int(r['allele_id'])}" if pd.notna(r['allele_id']) else f"{sample_id}|NA"
        else:
            row['allele_key'] = f"{sample_id}|NA"
        out_rows.append(row)

    feat = pd.DataFrame(out_rows)
    for sym in sorted(all_symbols):
        feat[f'tvr_frac_{sym}'] = feat['tvr_consensus'].fillna('').apply(lambda s: symbol_fraction(s, sym))
    return feat


def summarize_sample_qc(alleles):
    grp = alleles.groupby('sample_id')
    qc = grp.agg(
        total_alleles=('allele_id', 'count'),
        nonambiguous_alleles=('is_ambiguous_chr', lambda x: int((~x).sum())),
        pass_basic_qc_alleles=('pass_basic_qc', 'sum'),
        pass_basic_qc_nonambiguous=('pass_basic_qc', lambda x: int(x.sum())),
        unique_chr_arms_total=('primary_chr_arm', 'nunique'),
        unique_chr_arms_nonambiguous=('primary_chr_arm', lambda s: s[~alleles.loc[s.index, 'is_ambiguous_chr']].nunique()),
        unique_p_arms=('arm', lambda s: int((s == 'p').sum())),
        unique_q_arms=('arm', lambda s: int((s == 'q').sum())),
        median_TL_bp=('TL_p75', 'median'),
        median_TVR_bp=('tvr_len', 'median'),
        median_tvr_nonC_fraction=('tvr_nonC_fraction', 'median'),
        median_mapq=('mapq_median', 'median')
    ).reset_index()
    return qc


def make_sample_arm_summary(alleles):
    good = alleles[alleles['pass_basic_qc'] & (~alleles['is_ambiguous_chr']) & alleles['arm'].isin(['p', 'q'])].copy()
    grp = good.groupby(['sample_id', 'primary_chr_arm', 'chrom', 'arm'])
    out = grp.agg(
        allele_count=('allele_id', 'count'),
        median_TL_bp=('TL_p75', 'median'),
        mean_TL_bp=('TL_p75', 'mean'),
        median_TVR_bp=('tvr_len', 'median'),
        mean_TVR_bp=('tvr_len', 'mean'),
        median_tvr_entropy=('tvr_symbol_entropy', 'median'),
        median_tvr_nonC_fraction=('tvr_nonC_fraction', 'median'),
        median_tvr_transition_rate=('tvr_transition_rate', 'median'),
        median_support_reads=('support_n', 'median'),
        median_mapq=('mapq_median', 'median')
    ).reset_index()
    frac_cols = [c for c in good.columns if c.startswith('tvr_frac_')]
    if frac_cols:
        extra = grp[frac_cols].mean().reset_index()
        out = out.merge(extra, on=['sample_id', 'primary_chr_arm', 'chrom', 'arm'], how='left')
    return out.sort_values(['arm', 'primary_chr_arm', 'sample_id'], key=lambda s: s.map(chr_sort_key) if s.name == 'primary_chr_arm' else s)


def make_cohort_arm_summary(sample_arm):
    grp = sample_arm.groupby(['primary_chr_arm', 'chrom', 'arm'])
    out = grp.agg(
        sample_count=('sample_id', 'nunique'),
        row_count=('sample_id', 'count'),
        median_of_sample_median_TL_bp=('median_TL_bp', 'median'),
        mean_of_sample_median_TL_bp=('median_TL_bp', 'mean'),
        median_of_sample_median_TVR_bp=('median_TVR_bp', 'median'),
        mean_of_sample_median_TVR_bp=('median_TVR_bp', 'mean'),
        median_of_sample_tvr_nonC_fraction=('median_tvr_nonC_fraction', 'median'),
        median_of_sample_tvr_entropy=('median_tvr_entropy', 'median')
    ).reset_index()
    frac_cols = [c for c in sample_arm.columns if c.startswith('tvr_frac_')]
    if frac_cols:
        extra = grp[frac_cols].mean().reset_index()
        out = out.merge(extra, on=['primary_chr_arm', 'chrom', 'arm'], how='left')
    return out.sort_values(['arm', 'primary_chr_arm'], key=lambda s: s.map(chr_sort_key) if s.name == 'primary_chr_arm' else s)


def _prep_centered(df, feature, y='TL_p75', group='primary_chr_arm', min_group_n=3):
    tmp = df[[feature, y, group]].copy().dropna()
    if tmp.empty:
        return tmp
    group_sizes = tmp.groupby(group).size()
    valid_groups = group_sizes[group_sizes >= min_group_n].index
    tmp = tmp[tmp[group].isin(valid_groups)].copy()
    if tmp.empty:
        return tmp
    tmp['x_centered'] = tmp[feature] - tmp.groupby(group)[feature].transform('mean')
    tmp['y_centered'] = tmp[y] - tmp.groupby(group)[y].transform('mean')
    return tmp


def association_table(df, feature_cols, y='TL_p75', group='primary_chr_arm', min_group_n=3, level='allele'):
    rows = []
    for feature in feature_cols:
        if feature not in df.columns:
            continue
        tmp = df[[feature, y, group]].copy().dropna()
        n = len(tmp)
        rho, p = (np.nan, np.nan)
        if n >= 3 and tmp[feature].nunique() > 1 and tmp[y].nunique() > 1:
            rho, p = spearmanr(tmp[feature], tmp[y])
        centered = _prep_centered(df, feature, y=y, group=group, min_group_n=min_group_n)
        wrho, wp = (np.nan, np.nan)
        wn = len(centered)
        if wn >= 3 and centered['x_centered'].nunique() > 1 and centered['y_centered'].nunique() > 1:
            wrho, wp = spearmanr(centered['x_centered'], centered['y_centered'])
        rows.append({
            'level': level,
            'feature': feature,
            'n_global': n,
            'spearman_rho_global': rho,
            'p_global': p,
            'n_within_arm': wn,
            'spearman_rho_within_arm': wrho,
            'p_within_arm': wp
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out['abs_rho_within_arm'] = out['spearman_rho_within_arm'].abs()
        out = out.sort_values(['abs_rho_within_arm', 'p_within_arm', 'spearman_rho_global'], ascending=[False, True, False])
    return out


def per_arm_assoc(sample_arm, feature='median_TVR_bp', y='median_TL_bp', min_samples=3):
    rows = []
    for arm, sub in sample_arm.groupby('primary_chr_arm'):
        tmp = sub[[feature, y]].dropna()
        n = len(tmp)
        rho, p = (np.nan, np.nan)
        if n >= min_samples and tmp[feature].nunique() > 1 and tmp[y].nunique() > 1:
            rho, p = spearmanr(tmp[feature], tmp[y])
        rows.append({'primary_chr_arm': arm, 'n_samples': n, 'feature': feature, 'rho': rho, 'p': p})
    return pd.DataFrame(rows).sort_values('primary_chr_arm', key=lambda s: s.map(chr_sort_key))


def generate_inference_text(assoc_allele, assoc_sample_arm, out_fp):
    lines = []
    lines.append('Interpretation notes')
    lines.append('====================')
    lines.append('These associations are correlation-based and do not establish causality.')
    lines.append('Within-arm correlations are usually more informative than global correlations because they reduce the confounding effect of intrinsic chromosome-arm length differences.')
    lines.append('')
    def top_line(df, label):
        if df.empty:
            return f'{label}: no analyzable features.'
        top = df.dropna(subset=['spearman_rho_within_arm']).sort_values(['abs_rho_within_arm', 'p_within_arm'], ascending=[False, True]).head(1)
        if top.empty:
            return f'{label}: no analyzable within-arm associations.'
        r = top.iloc[0]
        direction = 'positive' if r['spearman_rho_within_arm'] > 0 else 'negative'
        return (
            f"{label}: strongest within-arm association = {r['feature']} "
            f"(rho={r['spearman_rho_within_arm']:.3f}, p={r['p_within_arm']:.3g}, {direction}; n={int(r['n_within_arm'])})."
        )
    lines.append(top_line(assoc_allele, 'Allele-level'))
    lines.append(top_line(assoc_sample_arm, 'Sample-arm-level'))
    lines.append('')
    lines.append('Recommended emphasis: prioritize sample-arm-level and within-arm metrics when summarizing the TVR–TL relationship across multiple samples.')
    with open(out_fp, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


def symbol_at_bp(row, center_bp, max_bp):
    tl = row['TL_p75']
    if pd.isna(tl) or tl <= 0 or center_bp > min(tl, max_bp):
        return '_'

    full_seq = row.get('_full_consensus', '')
    if full_seq:
        idx = min(len(full_seq) - 1, int((center_bp / tl) * len(full_seq)))
        return full_seq[idx]

    tvr_seq = row.get('tvr_consensus', '')
    tvr_len = row.get('tvr_len', np.nan)
    if tvr_seq and pd.notna(tvr_len) and tvr_len > 0 and center_bp <= min(tvr_len, tl):
        idx = min(len(tvr_seq) - 1, int((center_bp / tvr_len) * len(tvr_seq)))
        return tvr_seq[idx]

    return 'C'


def draw_variant_map(df, out_png, out_pdf, title, arm, max_bp=15000, bin_bp=10):
    sub = df[(df['arm'] == arm) & (~df['is_ambiguous_chr']) & (df['pass_basic_qc'])].copy()
    sub = sub.sort_values(['primary_chr_arm', 'sample_id', 'TL_p75', 'allele_id'], ascending=[True, True, False, True], key=lambda s: s.map(chr_sort_key) if s.name == 'primary_chr_arm' else s)
    if sub.empty:
        return False
    n_bins = max_bp // bin_bp
    plot_rows = []
    labels = []
    for _, row in sub.iterrows():
        sym_row = []
        for j in range(n_bins):
            center = (j + 0.5) * bin_bp
            sym_row.append(symbol_at_bp(row, center, max_bp))
        plot_rows.append(sym_row)
        labels.append(f"{row['sample_id']} | ({int(row['allele_id'])}) {row['primary_chr_arm']}")

    rgb = np.zeros((len(plot_rows), n_bins, 3), dtype=np.uint8)
    for i, row in enumerate(plot_rows):
        for j, sym in enumerate(row):
            rgb[i, j] = (np.array(to_rgb(PALETTE.get(sym, '#000000'))) * 255).astype(np.uint8)

    fig_h = max(6, len(plot_rows) * 0.28)
    plt.figure(figsize=(16, fig_h), dpi=200)
    ax = plt.gca()
    ax.imshow(rgb, aspect='auto', interpolation='nearest', extent=[0, max_bp, len(plot_rows), 0])
    ax.set_yticks(np.arange(len(labels)) + 0.5)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('distance from subtelomere boundary (bp)', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.set_xlim(0, max_bp)
    ax.set_xticks(np.arange(0, max_bp + 1, 1000))
    ax.tick_params(axis='x', labelsize=9)
    ax.grid(axis='x', color='white', linestyle='--', linewidth=0.3, alpha=0.6)
    for y in range(len(labels) + 1):
        ax.axhline(y, color='black', linewidth=0.25)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()
    return True


def scatter_plot(df, x, y, out_fp, title, xlabel=None, ylabel=None):
    tmp = df[[x, y]].dropna().copy()
    if tmp.empty:
        return False
    plt.figure(figsize=(6, 5), dpi=200)
    plt.scatter(tmp[x], tmp[y], s=18, alpha=0.7)
    plt.xlabel(xlabel or x)
    plt.ylabel(ylabel or y)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fp, bbox_inches='tight')
    plt.close()
    return True


def run(args):
    os.makedirs(args.outdir, exist_ok=True)
    records = load_manifest(args)
    if not records:
        raise ValueError('No tlens_by_allele.tsv files found.')

    final_map = {}
    if args.final_manifest:
        fm = pd.read_csv(args.final_manifest, sep='\t')
        if not {'sample_id', 'final_zip'}.issubset(fm.columns):
            raise ValueError('Final manifest must contain columns: sample_id and final_zip')
        final_map = dict(zip(fm['sample_id'].astype(str), fm['final_zip'].astype(str)))

    all_tables = []
    for rec in records:
        sample_id = rec['sample_id']
        tlens = rec['tlens']
        final_zip = rec.get('final_zip') or final_map.get(sample_id)
        try:
            tbl = build_allele_table(sample_id, tlens, final_zip=final_zip, min_support=args.min_support, min_mapq=args.min_mapq)
            cons = load_consensus(final_zip)
            tbl['_full_consensus'] = tbl['allele_id'].apply(lambda aid: cons.get(int(aid), '') if pd.notna(aid) else '')
            all_tables.append(tbl)
        except Exception as e:
            fail_row = pd.DataFrame([{
                'sample_id': sample_id,
                'tlens_file': tlens,
                'final_zip': final_zip or '',
                'processing_error': str(e)
            }])
            fail_row.to_csv(os.path.join(args.outdir, f'{sample_id}.processing_error.tsv'), sep='\t', index=False)

    if not all_tables:
        raise ValueError('No sample could be processed successfully.')

    alleles = pd.concat(all_tables, ignore_index=True)
    alleles.sort_values(['arm', 'primary_chr_arm', 'sample_id', 'allele_id'], inplace=True, key=lambda s: s.map(chr_sort_key) if s.name == 'primary_chr_arm' else s)
    alleles.to_csv(os.path.join(args.outdir, 'merged_allele_tvr_features.tsv'), sep='\t', index=False)

    qc = summarize_sample_qc(alleles)
    qc.to_csv(os.path.join(args.outdir, 'sample_qc_summary.tsv'), sep='\t', index=False)

    sample_arm = make_sample_arm_summary(alleles)
    sample_arm.to_csv(os.path.join(args.outdir, 'sample_arm_tvr_summary.tsv'), sep='\t', index=False)

    cohort_arm = make_cohort_arm_summary(sample_arm)
    cohort_arm.to_csv(os.path.join(args.outdir, 'cohort_arm_tvr_summary.tsv'), sep='\t', index=False)

    assoc_features = [c for c in KEY_FEATURES if c in alleles.columns]
    assoc_features += [c for c in alleles.columns if c.startswith('tvr_frac_')]
    good_alleles = alleles[(alleles['pass_basic_qc']) & (~alleles['is_ambiguous_chr']) & (alleles['arm'].isin(['p', 'q']))].copy()
    assoc_allele = association_table(good_alleles, assoc_features, y='TL_p75', group='primary_chr_arm', min_group_n=args.min_group_n, level='allele')
    assoc_allele.to_csv(os.path.join(args.outdir, 'TVR_TL_association_allele_level.tsv'), sep='\t', index=False)

    sample_arm_assoc_features = [
        'median_TVR_bp', 'median_tvr_entropy', 'median_tvr_nonC_fraction',
        'median_tvr_transition_rate', 'median_support_reads', 'median_mapq'
    ] + [c for c in sample_arm.columns if c.startswith('tvr_frac_')]
    assoc_sample_arm = association_table(sample_arm, sample_arm_assoc_features, y='median_TL_bp', group='primary_chr_arm', min_group_n=args.min_group_n, level='sample_arm')
    assoc_sample_arm.to_csv(os.path.join(args.outdir, 'TVR_TL_association_sample_arm_level.tsv'), sep='\t', index=False)

    per_arm = per_arm_assoc(sample_arm, feature='median_TVR_bp', y='median_TL_bp', min_samples=args.min_group_n)
    per_arm.to_csv(os.path.join(args.outdir, 'per_chr_arm_TVRlen_vs_TL_correlation.tsv'), sep='\t', index=False)

    generate_inference_text(assoc_allele, assoc_sample_arm, os.path.join(args.outdir, 'TVR_TL_inference_summary.txt'))

    draw_variant_map(
        alleles,
        os.path.join(args.outdir, 'merged_nonambiguous_p_variant_map.png'),
        os.path.join(args.outdir, 'merged_nonambiguous_p_variant_map.pdf'),
        'Merged non-ambiguous p-arm telomere variant map',
        arm='p', max_bp=args.max_bp, bin_bp=args.bin_bp
    )
    draw_variant_map(
        alleles,
        os.path.join(args.outdir, 'merged_nonambiguous_q_variant_map.png'),
        os.path.join(args.outdir, 'merged_nonambiguous_q_variant_map.pdf'),
        'Merged non-ambiguous q-arm telomere variant map',
        arm='q', max_bp=args.max_bp, bin_bp=args.bin_bp
    )

    plot_dir = os.path.join(args.outdir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    scatter_plot(good_alleles, 'tvr_len', 'TL_p75', os.path.join(plot_dir, 'allele_TL_vs_TVRlen.png'), 'Allele-level TL vs TVR length', 'TVR length (bp)', 'TL_p75 (bp)')
    scatter_plot(good_alleles, 'tvr_nonC_fraction', 'TL_p75', os.path.join(plot_dir, 'allele_TL_vs_TVR_nonC_fraction.png'), 'Allele-level TL vs TVR non-C fraction', 'TVR non-C fraction', 'TL_p75 (bp)')
    scatter_plot(good_alleles, 'tvr_symbol_entropy', 'TL_p75', os.path.join(plot_dir, 'allele_TL_vs_TVR_entropy.png'), 'Allele-level TL vs TVR entropy', 'TVR entropy', 'TL_p75 (bp)')
    scatter_plot(sample_arm, 'median_TVR_bp', 'median_TL_bp', os.path.join(plot_dir, 'sample_arm_TL_vs_TVRlen.png'), 'Sample-arm median TL vs TVR length', 'Median TVR length (bp)', 'Median TL (bp)')
    scatter_plot(sample_arm, 'median_tvr_nonC_fraction', 'median_TL_bp', os.path.join(plot_dir, 'sample_arm_TL_vs_TVR_nonC_fraction.png'), 'Sample-arm median TL vs TVR non-C fraction', 'Median TVR non-C fraction', 'Median TL (bp)')

    readme = os.path.join(args.outdir, 'README.txt')
    with open(readme, 'w', encoding='utf-8') as f:
        f.write('Outputs\n')
        f.write('=======\n')
        f.write('merged_allele_tvr_features.tsv: merged allele-level table across all successfully processed samples.\n')
        f.write('sample_qc_summary.tsv: per-sample allele counts and completeness after basic QC.\n')
        f.write('sample_arm_tvr_summary.tsv: per-sample per-chromosome-arm summary restricted to non-ambiguous p/q alleles passing QC.\n')
        f.write('cohort_arm_tvr_summary.tsv: cohort-level summary aggregated from sample-arm medians.\n')
        f.write('merged_nonambiguous_p_variant_map.* and merged_nonambiguous_q_variant_map.*: separate p-arm and q-arm TVR maps.\n')
        f.write('TVR_TL_association_allele_level.tsv: association of TVR features with TL_p75 across filtered alleles.\n')
        f.write('TVR_TL_association_sample_arm_level.tsv: association of TVR features with median TL at the sample-arm level.\n')
        f.write('per_chr_arm_TVRlen_vs_TL_correlation.tsv: per-arm correlation across samples for median TVR length vs median TL.\n')
        f.write('TVR_TL_inference_summary.txt: concise textual interpretation of the strongest within-arm associations.\n')
        f.write('plots/: selected scatter plots summarizing TVR-TL relationships.\n')
        f.write('\n')
        f.write(f'Basic QC thresholds used: min_support={args.min_support}, min_mapq={args.min_mapq}.\n')
        f.write('Ambiguous chromosome-arm assignments are excluded from the maps and the association analyses.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Multi-sample Telogator TVR map + TVR-TL association analysis')
    parser.add_argument('--manifest', help='TSV with columns: sample_id, tlens, optional final_zip')
    parser.add_argument('--final-manifest', help='Optional TSV with columns: sample_id, final_zip')
    parser.add_argument('-i', '--inputs', nargs='+', help='List of tlens_by_allele.tsv files')
    parser.add_argument('--scan-root', help='Recursively scan for *tlens_by_allele.tsv files')
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('--min-support', type=int, default=3)
    parser.add_argument('--min-mapq', type=float, default=0.0)
    parser.add_argument('--min-group-n', type=int, default=3, help='Minimum rows or samples for within-arm association')
    parser.add_argument('--max-bp', type=int, default=15000)
    parser.add_argument('--bin-bp', type=int, default=10)
    args = parser.parse_args()
    run(args)
