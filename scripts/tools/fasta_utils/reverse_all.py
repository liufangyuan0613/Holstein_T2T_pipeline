from Bio import SeqIO
import sys

def reverse_complement_sequences(input_file, output_file):
    """
    Generates reverse complement sequences for each sequence in the input fasta file
    and writes them to the output fasta file.

    Args:
    input_file (str): Path to the input fasta file.
    output_file (str): Path to the output fasta file.
    """
    # Load sequences from the input fasta file
    sequences = SeqIO.parse(input_file, 'fasta')
    
    # Generate reverse complement for each sequence
    rev_comp_sequences = []
    for sequence in sequences:
        rev_comp_seq = sequence.seq.reverse_complement()
        sequence.seq = rev_comp_seq
        sequence.id += "_revcomp"
        sequence.description = "Reverse complement of " + sequence.id
        rev_comp_sequences.append(sequence)
    
    # Write reverse complement sequences to the output file
    SeqIO.write(rev_comp_sequences, output_file, 'fasta')
    print(f"Reverse complement sequences have been saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta> <output_fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    reverse_complement_sequences(input_fasta, output_fasta)
