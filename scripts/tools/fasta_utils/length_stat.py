#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Calculate the length of each chromosome from a FASTA file")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file to write chromosome lengths")
    
    args = parser.parse_args()
    fasta_file = args.fasta
    output_file = args.output

    with open(output_file, 'w') as out_f:
        out_f.write("Chromosome\tLength\n")
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom = record.id
            length = len(record.seq)
            out_f.write(f"{chrom}\t{length}\n")

if __name__ == "__main__":
    main()
