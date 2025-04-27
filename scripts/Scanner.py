import csv
import re

def load_fasta_sequence(fasta_path):
    """Loads the whole FASTA sequence into one string."""
    sequence_lines = []
    with open(fasta_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence_lines.append(line.strip())
    return ''.join(sequence_lines)

def parse_mutation_name(mutation_name, mutation_type):
    """
    Parses the mutation name and returns the start and end coordinates.
    Additionally, for Indels, it will return the inserted sequence.
    """

    # Case 1: Single Nucleotide Variant (SNV)
    if mutation_type == 'single nucleotide variant':
        # Case 1.1: Standard SNV (e.g., c.694C>T, c.1234A>G)
        match = re.search(r'c\.(\d+)([A-Za-z])>([A-Za-z])', mutation_name)
        if match:
            start = int(match.group(1))  # Start position (e.g., 694 or 1234)
            ref_allele = match.group(2)  # Reference allele (e.g., C or A)
            alt_allele = match.group(3)  # Alternate allele (e.g., T or G)
            return start, start, ref_allele, alt_allele

        # Case 1.2: Splice site mutation (e.g., c.700+4G>T, c.581-176A>T)
        match = re.search(r'c\.(\d+)([+-]\d+)([A-Za-z])>([A-Za-z])', mutation_name)
        if match:
            exon_position = int(match.group(1))  # Exon base
            splice_offset = int(match.group(2))  # +4 or -176 (converted to int)
            final_position = exon_position + splice_offset  # Compute final position
            ref_allele = match.group(3)
            alt_allele = match.group(4)
            return final_position, final_position, ref_allele, alt_allele
        else:
            return -1, -1, -1 , -1

    # Case 2: Insertion (Indel)
    elif mutation_type == 'Indel':
        # Regular expression for Indels (e.g., c.80_83delinsTGCTGTAAACTGTAACTGTAAA)
        match = re.search(r'c\.(\d+)_?(\d+)(delins[A-Za-z]+)', mutation_name)
        if match:
            start = int(match.group(1))  # Start position (e.g., 80)
            end = int(match.group(2))    # End position (e.g., 83)
            inserted_sequence = match.group(3)[5:]  # Get the inserted bases after 'delins'
            return start, end, inserted_sequence
        else:
            return -1, -1, -1

    # Case 3: Deletion (e.g., c.1234_1235del)
    elif mutation_type == 'Deletion':
            match = re.search(r'c\.(\d+[+-]?\d*)_(\d+[+-]?\d*)del', mutation_name)
            if match:
                start = match.group(1)  # '361-5'
                end = match.group(2)    # '361-1'

                # Now, remove '+' or '-' for easier processing
                start_clean = int(re.sub(r'[+-]', '', start))
                end_clean = int(re.sub(r'[+-]', '', end))

                return start_clean, end_clean  # No sequence for deletion
            else:
                return -1, -1

    else:
        raise ValueError(f"Mutation type {mutation_type} is not recognized.")


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

    print("Loading patient DNA sequence into memory...")
    patient_sequence = load_fasta_sequence(fasta_path)
    print(f"Sequence loaded: {len(patient_sequence)} bases")

    for variant in variants:
        if variant['Chromosome'] == chromosome_input:
            # Check if the variant matches with the patient's DNA sequence (simplified matching)
            mutation_name = variant['Name'].strip()
            mutation_type = variant['Type'].strip()  # Remove extra spaces if any
            if mutation_type == "Indel":
                # Parse the name, not using the s & e, but if alternate allele not available, then insert_sequence is valuable
                s, e, ins_sequence  = parse_mutation_name(mutation_name, mutation_type)
                start_position = variant['Start']
                end_position = variant['Stop']
                new_sequence = variant['AlternateAllele'].strip()
                
                # Only checks for mutation if the insert section, start, and stop variables are known
                if start_position == -1 or end_position == -1 or (new_sequence == 'na' and ins_sequence == -1):
                    continue
                    
                # Assurance for no out of bounds error for the string
                if len(patient_sequence) >= end_position:
                    # if new_sequence == 'na' use ins_sequence
                    if new_sequence == 'na':
                        new_or_ins_sequence = ins_sequence
                    else:
                        new_or_ins_sequence = new_sequence
                    # Takes the sequence from the chromosome looking specifically at the start and stop sections
                    if new_or_ins_sequence in patient_sequence[start_position-1: end_position]:
                        '''
                        Put stuff into findings
                        '''
                # else there is no mutation present and the loop can continue
            elif mutation_type == "Deletion":
                start_position, end_position  = parse_mutation_name(mutation_name, mutation_type)
                if start_position == -1 or end_position == -1:
                    continue
            elif mutation_type == "single nucleotide variant":
                start_position, end_position, ref_allele, alt_allele  = parse_mutation_name(mutation_name, mutation_type)
                if start_position == -1 or end_position == -1 or ref_allele == -1 or alt_allele == -1:
                    continue
            else: continue
            
            
            # print("Loading patient DNA...")
            # print(variant["Chromosome"])
            # extracted_sequence = extract_sequence(fasta_path, start_position, end_position)
            extracted_sequence = patient_sequence[start_position-1:end_position]
            # print(extracted_sequence)
            gene_symbol = variant['GeneSymbol'].strip()

            # Example logic to simulate DNA testing (this could be more complex)
            if gene_symbol in extracted_sequence:
                findings.append({
                    'Gene': gene_symbol,
                    'MutationType': mutation_type,
                    'MutationName': mutation_name,
                    'Position': variant['Start'],
                    'Condition': variant['Condition'],
                    'ClinicalSignificance': variant['ClinicalSignificance']
                })

    return findings



def main():
    chromosome_input = input("Enter the gene being presented: ")
    fasta_path = input("Enter Patient Sequence File: ")
    variants_csv_path = "C:\\Users\\rogue\\PantherHack2025\\PantherHack25\\outputs\\pathogenic_variants_sorted.csv"

    print("Loading patient DNA...")
    # dna_sequence = load_fasta_sequence(fasta_path)

    print("Loading known pathogenic variants...")
    variants = load_variants(variants_csv_path)
    

    print("Analyzing for mutations...")
    findings = analyze_variants(fasta_path, variants, chromosome_input)
    if findings:
        print(f"\nMutations found for Chromosome {chromosome_input}:")
        for finding in findings:
            print(finding)
    else:
        print(f"\nNo mutations found for Chromosome {chromosome_input}.")


if __name__ == "__main__":
    main()
