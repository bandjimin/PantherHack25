# main.py

from Bio import SeqIO

# Function to search for patterns in a DNA sequence
def find_patterns(sequence, pattern):
    if pattern in sequence:
        return f"Pattern '{pattern}' found in sequence!"
    else:
        return f"Pattern '{pattern}' not found in sequence."

# Function to load DNA sequences from a file
def load_dna_sequence(file_path):
    # Load DNA sequence from a FASTA file
    seq_record = SeqIO.read(file_path, "fasta")
    return str(seq_record.seq)

# Example usage
if __name__ == "__main__":
    # Load a sample DNA sequence (this can be a path to a .fasta file)
    dna_sequence = load_dna_sequence("example_sequence.fasta")  # Update with your file path
    print(f"Loaded DNA sequence: {dna_sequence}")

    # Define the pattern to search for
    pattern = "AGCT"

    # Search for the pattern in the sequence
    result = find_patterns(dna_sequence, pattern)
    print(result)
