#!/usr/bin/env python3
"""05_gapfill_polish/gapfill_local_assembly.py — close assembly gaps with ONT local assembly.

Logic (manuscript STAR Methods, 'Gap filling and genome polishing'):
  1. ONT reads are aligned to the chromosome assembly (minimap2 -ax map-ont --secondary=no).
  2. A read is kept as gap-filling evidence if it is unmapped, poorly aligned
     (alignment identity < 80% or aligned read coverage < 80%), or uniquely spans both
     flanking sequences of a gap.
  3. For each gap, candidate reads from the flanking regions are collected and subjected
     to local iterative assembly (hifiasm).
  4. Locally assembled contigs are aligned back (minimap2 -ax asm5 --eqx) and a contig is
     accepted for gap replacement only when it bridges the gap and shows reliable overlap
     with BOTH flanking sequences.

Usage:
  python3 gapfill_local_assembly.py --gaps gaps.bed --paf ont2gap.paf \
      --reads ont.fastq --genome A-1.haphic.fa --out A-1.gapfilled.fa \
      --flank 50000 --threads 64 [--dry-run]
"""
import argparse, os, re, subprocess, sys
from collections import defaultdict

def read_fai(fa):
    fai = fa + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", fa], check=True)
    out = {}
    with open(fai) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            out[p[0]] = int(p[1])
    return out

def read_fasta(fa):
    seqs, name, buf = {}, None, []
    with open(fa) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs

def parse_gaps(bed):
    gaps = []
    with open(bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                gaps.append((p[0], int(p[1]), int(p[2])))
    return gaps

def parse_paf(paf, min_span=0.8):
    """yield (read, chrom, rstart, rend, aln_len, read_len) for usable alignments"""
    with open(paf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12: continue
            rlen, rstart, rend = int(p[1]), int(p[2]), int(p[3])
            alen, blen = int(p[9]), int(p[10])
            if rlen == 0: continue
            cov = (rend - rstart) / rlen
            ident = blen / alen if alen else 0
            if cov >= 0.8 and ident >= 0.8:
                yield (p[0], p[5], int(p[7]), int(p[8]), alen, rlen)

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaps", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    lens = read_fai(a.genome)
    seqs = read_fasta(a.genome)
    gaps = parse_gaps(a.gaps)
    if not gaps:
        print("no gaps found; copying genome", file=sys.stderr)
        with open(a.out, "w") as o:
            for n, s in seqs.items():
                o.write(f">{n}\n{s}\n")
        return

    # index alignments by chromosome
    aln = defaultdict(list)
    for r in parse_paf(a.paf):
        aln[r[1]].append(r)

    work = "_gapfill_work"; os.makedirs(work, exist_ok=True)
    replacements = {}

    for (chrom, gs, ge) in gaps:
        clen = lens[chrom]
        lf = (max(0, gs - a.flank), gs)          # left flank window
        rf = (ge, min(clen, ge + a.flank))       # right flank window
        # candidate reads spanning into both flanks (or unmapped/poor -> handled upstream)
        cand = set()
        for (read, c, s, e, alen, rlen) in aln.get(chrom, []):
            if (s <= lf[1] and e >= lf[0]) or (s <= rf[1] and e >= rf[0]):
                cand.add(read)
        tag = f"{chrom}_{gs}_{ge}"
        print(f"[{tag}] {len(cand)} candidate reads", file=sys.stderr)
        if a.dry_run or not cand:
            continue
        # dump candidate reads
        subfq = os.path.join(work, tag + ".fq")
        subprocess.run(f"seqkit grep -n -f <(printf '%s\\n' {' '.join(sorted(cand))}) {a.reads} > {subfq}",
                       shell=True, check=True, executable="/bin/bash")
        # local assembly
        outg = os.path.join(work, tag + ".hifiasm")
        subprocess.run(["hifiasm", "-o", outg, "-t", str(min(a.threads, 16)), subfq],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        gfa = outg + ".bp.p_ctg.gfa"
        if not os.path.exists(gfa):
            print(f"[{tag}] local assembly empty; gap left as-is", file=sys.stderr)
            continue
        subprocess.run(f"awk '/^S/{{print \">\"$2\"\\n\"$3}}' {gfa} > {outg}.fa",
                       shell=True, check=True, executable="/bin/bash")
        # align local contigs back, require bridging with overlap on both flanks
        paf2 = os.path.join(work, tag + ".back.paf")
        subprocess.run(f"minimap2 -ax asm5 --eqx {a.genome} {outg}.fa > {paf2}",
                       shell=True, check=True, executable="/bin/bash")
        best = None
        with open(paf2) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12 or p[5] != chrom: continue
                qs, qe = int(p[7]), int(p[8])
                if qs < gs - 100 and qe > ge + 100:   # bridges the gap with overlap
                    qname, qseq_len = p[0], int(p[1])
                    best = (qname, qs, qe, p[4])
                    break
        if best:
            replacements[(chrom, gs, ge)] = best
            print(f"[{tag}] bridged by {best[0]} ({best[2]-best[1]} bp)", file=sys.stderr)
        else:
            print(f"[{tag}] no bridging contig; gap left as-is", file=sys.stderr)

    # apply replacements
    loc = {}
    for (chrom, gs, ge), (qname, qs, qe, strand) in replacements.items():
        gfa_fa = os.path.join(work, f"{chrom}_{gs}_{ge}.hifiasm.fa")
        s = read_fasta(gfa_fa)[qname]
        insert_start, insert_end = (gs - qs), (ge - qs)
        piece = s[insert_start:insert_end]
        if strand == "-":
            piece = revcomp(piece)
        loc[chrom] = seqs[chrom][:gs] + piece + seqs[chrom][ge:]

    n_closed = 0
    with open(a.out, "w") as o:
        for name, s in seqs.items():
            s2 = loc.get(name, s)
            for (chrom, gs, ge), _ in replacements.items():
                if chrom == name:
                    n_closed += 1
            for i in range(0, len(s2), 60):
                pass
            o.write(f">{name}\n")
            for i in range(0, len(s2), 80):
                o.write(s2[i:i+80] + "\n")
    print(f"gaps closed: {n_closed}/{len(gaps)}", file=sys.stderr)

if __name__ == "__main__":
    main()
