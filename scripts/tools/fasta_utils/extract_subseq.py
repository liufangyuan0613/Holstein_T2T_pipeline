from Bio import SeqIO
import sys

def extract_sequence(input_file, output_file, seq_id, start, end):
    """
    Extracts a subsequence from a fasta file and saves it to a new file.

    Args:
    input_file (str): Path to the input fasta file.
    output_file (str): Path to the output fasta file.
    seq_id (str): The ID of the sequence to extract.
    start (int): The starting position of the subsequence (1-based).
    end (int): The ending position of the subsequence (inclusive).
    """
    # Load the sequences from the fasta file
    records = SeqIO.parse(input_file, 'fasta')
    
    # Look for the specified sequence and extract the subsequence
    for record in records:
        if record.id == seq_id:
            # Adjust for 0-based indexing used in Python
            subseq = record.seq[start-1:end]
            # Create a new record for the subsequence
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            new_record = SeqRecord(Seq(subseq), id=seq_id + "_extract", description="Extracted from {} to {}".format(start, end))
            # Write the subsequence to the output file
            SeqIO.write(new_record, output_file, 'fasta')
            print("Subsequence extracted and saved to", output_file)
            break
    else:
        print("Sequence ID not found in the provided fasta file.")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <input_fasta> <output_fasta> <sequence_id> <start_position> <end_position>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    sequence_id = sys.argv[3]
    start_pos = int(sys.argv[4])
    end_pos = int(sys.argv[5])

    extract_sequence(input_fasta, output_fasta, sequence_id, start_pos, end_pos)
