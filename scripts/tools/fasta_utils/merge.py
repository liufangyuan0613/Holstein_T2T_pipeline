#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def merge_fasta(file1, file2, output_file):
    with open(output_file, 'w') as out_f:
        for fasta_file in [file1, file2]:
            with open(fasta_file, 'r') as in_f:
                for line in in_f:
                    out_f.write(line)

def main():
    parser = argparse.ArgumentParser(description="Merge two FASTA files into one")
    parser.add_argument("-f1", "--file1", required=True, help="First input FASTA file")
    parser.add_argument("-f2", "--file2", required=True, help="Second input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    file1 = args.file1
    file2 = args.file2
    output_file = args.output

    merge_fasta(file1, file2, output_file)

if __name__ == "__main__":
    main()

