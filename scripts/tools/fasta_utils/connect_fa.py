from Bio import SeqIO
import sys

def concatenate_sequences(file1, file2, output_file):
    """
    Concatenates sequences from file1 with corresponding sequences from file2
    in order and saves them to an output fasta file.

    Args:
    file1 (str): Path to the first input fasta file, containing the first part of the sequences.
    file2 (str): Path to the second input fasta file, containing the second part of the sequences.
    output_file (str): Path to the output fasta file to save concatenated sequences.
    """
    # Load sequences from both fasta files
    seqs1 = list(SeqIO.parse(file1, 'fasta'))
    seqs2 = list(SeqIO.parse(file2, 'fasta'))

    # Ensure the number of sequences is the same in both files
    if len(seqs1) != len(seqs2):
        raise ValueError("The number of sequences in each file must be the same.")

    concatenated_seqs = []
    # Iterate over pairs of sequences and concatenate them
    for seq1, seq2 in zip(seqs1, seqs2):
        new_seq = seq1.seq + seq2.seq  # Concatenate sequences
        new_id = f"{seq1.id}_concat_{seq2.id}"
        new_description = "Concatenation of two sequences"
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        concatenated_seq = SeqRecord(Seq(new_seq), id=new_id, description=new_description)
        concatenated_seqs.append(concatenated_seq)

    # Write the concatenated sequences to the output fasta file
    SeqIO.write(concatenated_seqs, output_file, 'fasta')
    print("Concatenation complete. Results saved to", output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_fasta1> <input_fasta2> <output_fasta>")
        sys.exit(1)
    
    input_fasta1 = sys.argv[1]
    input_fasta2 = sys.argv[2]
    output_fasta = sys.argv[3]

    concatenate_sequences(input_fasta1, input_fasta2, output_fasta)
