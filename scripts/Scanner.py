import csv
import re

# def load_fasta_sequence(fasta_path):
#     """Load DNA sequence from a FASTA file, skipping header."""
#     sequence = ''
#     with open(fasta_path, 'r') as file:
#         for line in file:
#             if not line.startswith('>'):
#                 sequence += line.strip()
#     return sequence.upper()

def parse_mutation_name(mutation_name):
    # Regular expression to capture the coordinates
    match = re.search(r'c\.(\d+)_?(\d+)', mutation_name)
    
    if match:
        start = int(match.group(1))  # First captured group is the start position
        end = int(match.group(2))    # Second captured group is the end position
        return start, end
    else:
        raise ValueError("Mutation name does not contain valid coordinates.")

def extract_sequence(fasta_path, start, end):
    """
    Extracts a specific sequence from the FASTA file between the start and end positions.
    """
    with open(fasta_path, 'r') as file:
        sequence = ""
        in_sequence = False
        current_position = 0

        for line in file:
            # Ignore the header line starting with '>'
            if line.startswith('>'):
                continue
            
            # Process sequence lines
            line = line.strip()
            
            if current_position + len(line) >= start:  # Check if we've reached the start position
                if current_position < end:  # Only include sequence up to the end position
                    sequence += line[start - current_position: end - current_position]
                else:
                    break  # No need to process further lines once we reach the end position
            current_position += len(line)

        return sequence


# Use this function to extract the sequence from chr13 based on gene position (start and end)


def load_variants(csv_path):
    variants = []
    with open(csv_path, 'r') as file:
        reader = csv.DictReader(file)
        # Clean up column names by stripping leading/trailing spaces
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        
        for row in reader:
            variants.append(row)
    return variants


def analyze_variants(fasta_path, variants, chromosome_input):
    findings = []

    # Loop through the variants
    for variant in variants:
        if variant['Chromosome'] == chromosome_input:
            # Check if the variant matches with the patient's DNA sequence (simplified matching)
            mutation_type = variant['Type'].strip()  # Remove extra spaces if any
            if mutation_type == "Indel":
                    
            elif mutation_type == "Deletion":
                        
            elif mutation_type == "Deletion":
                        
            else: continue
            mutation_name = variant['Name'].strip()
            start_position, end_position = parse_mutation_name(mutation_name)
            print("Loading patient DNA...")
            print(variant["Chromosome"])
            extracted_sequence = extract_sequence(fasta_path, start_position, end_position)
            print(extracted_sequence)
            gene_symbol = variant['GeneSymbol'].strip()

            # Example logic to simulate DNA testing (this could be more complex)
            if gene_symbol in extracted_sequence:
                findings.append({
                    'Gene': gene_symbol,
                    'MutationType': mutation_type,
                    'MutationName': mutation_name,
                    'Position': variant['Start'],
                    'Condition': variant['Condition'],
                    'ClinicalSignificance': variant['ClinicalSignificance'],
                    # "Extracted_Sequence": extracted_sequence
                })

    return findings



def main():
    chromosome_input = input("Enter the gene being presented: ")
    fasta_path = input("Enter Patient Sequence File: ")
    variants_csv_path = "C:\\Users\\rogue\\PantherHack2025\\PantherHack25\\outputs\\pathogenic_variants.csv"

    print("Loading patient DNA...")
    # dna_sequence = load_fasta_sequence(fasta_path)

    print("Loading known pathogenic variants...")
    variants = load_variants(variants_csv_path)
    

    print("Analyzing for mutations...")
    findings = analyze_variants(fasta_path, variants, chromosome_input)
    # if findings:
    #     print(f"\nMutations found for Chromosome {chromosome_input}:")
    #     for finding in findings:
    #         print(finding)
    # else:
    #     print(f"\nNo mutations found for Chromosome {chromosome_input}.")


if __name__ == "__main__":
    main()
