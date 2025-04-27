# main.py

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML



def kmp_search(sequence, pattern):
    positions = []
    lps = compute_lps(pattern)
    i = j = 0
    while i < len(sequence):
        if sequence[i] == pattern[j]:
            i += 1
            j += 1
            if j == len(pattern):
                positions.append(i - j)
                j = lps[j - 1]
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return positions

# Function to load DNA sequences from a file
def load_dna_sequence(file_path):
    # Load DNA sequence from a FASTA file
    with open(file_path, "r") as f:
        lines = f.readlines()
        dna_sequence =''.join(line.strip() for line in lines if not line.startswith('>'))
    # seq_record = SeqIO.read(file_path, "fasta")
    return dna_sequence[:1000000]

# Function to perform BLAST search against a reference genome (example: Human genome)
def blast_search(input_sequence):
    # Send the input sequence to NCBI's BLAST server for alignment
    result_handle = NCBIWWW.qblast("blastn", "nt", input_sequence)
    return result_handle

# Function to parse BLAST results and detect mutations, insertions, deletions
def parse_blast_results(blast_result_handle, reference_sequence):
    # Parse the BLAST result
    blast_records = NCBIPGI.parse(blast_result_handle)
    
    # For simplicity, let's assume we're only dealing with the first match
    for blast_record in blast_records:
        alignment = blast_record.alignments[0]
        aligned_seq = alignment.hsps[0].sbjct  # aligned portion of the input sequence
        
        # Compare aligned sequences with the reference sequence
        mutations, insertions, deletions = compare_sequences(aligned_seq, reference_sequence)
        
        return mutations, insertions, deletions

def compare_sequences(input_sequence, reference_sequence):
    """
    Compare two sequences and identify mutations, insertions, and deletions.
    """
    mutations = []
    insertions = []
    deletions = []
    
    # Ensure both sequences are the same length for comparison
    min_length = min(len(input_sequence), len(reference_sequence))
    
    # Compare sequences base by base
    for i in range(min_length):
        if input_sequence[i] != reference_sequence[i]:
            mutations.append((i, input_sequence[i], reference_sequence[i]))  # (index, input base, reference base)
    
    # Handle insertions or deletions if sequences are of different lengths
    if len(input_sequence) > len(reference_sequence):
        insertions = [(i, input_sequence[i]) for i in range(len(reference_sequence), len(input_sequence))]
    elif len(input_sequence) < len(reference_sequence):
        deletions = [(i, reference_sequence[i]) for i in range(len(input_sequence), len(reference_sequence))]
    
    return mutations, insertions, deletions

def check_disorders(dna_sequence):
    dna_sequence = dna_sequence.upper()

    # Define the known disorder-related patterns and thresholds
    disorders = [
        {
            "name": "Huntington's Disease",
            "pattern": "CAG",
            "min_repeats": 36,  # pathological threshold
            "repeats": True,
            "MID":False,
            "Location": null,
            "gene": null
        },
        {
            "name": "Fragile X Syndrome",
            "pattern": "CGG",
            "min_repeats": 200,
            "repeats": True,
            "MID":False,
            "Location": null,
            "gene": null
        },
        {
            "name": "Myotonic Dystrophy (Type 1)",
            "pattern": "CTG",
            "min_repeats": 50,
            "repeats": True,
            "MID":False,
            "Location": null,
            "gene": null

        },
        {
            "name": "Sickle Cell Anemia",
            "pattern": "GTG",  # GTG is the mutation (instead of GAG)
            "min_repeats": 1,
            "repeats": True,
            "MID":False,
            "Location": null,
            "gene": null
        },
        {
            "name": "Tay-Sachs Disease",
            "pattern": "TATC",
            "min_repeats": 1,
            "repeats": True,
            "MID":False,
            "Location": null,
            "gene": null
        },
        {
            "name": "Cystic Fibrosis (ΔF508)",
            "pattern": "TTC",  # Deletion of TTT (phenylalanine)
            "min_repeats": 0,  # absence indicates mutation, special case
            "repeats": True,
            "MID":False,
            "Location": 508,
            "gene": "CFTR"
        },
            "name": "Duchenne Muscular Dystrophy (DMD)",
            "pattern": "",
            "min_repeats": ,
            "repeats": False,
            "MID": True,
            "Location": 45-55,
            "gene": "Xp21.2"
    ]

    findings = []

    for disorder in disorders:
        name = disorder["name"]
        pattern = disorder["pattern"]
        min_repeats = disorder["min_repeats"]
        repeat = disorder["repeats"]
        MID = disorder["MID"]
        if repeat:
            if min_repeats > 1:
                max_repeat = find_max_consecutive_repeats_kmp(dna_sequence, pattern)
                findings.append(f"{name}: max {max_repeat} consecutive '{pattern}' repeats.")

                if max_repeat >= min_repeats:
                    findings.append(f"⚠️ Pattern for {name} detected ({max_repeat} consecutive repeats)")
            else:
                if name.startswith("Cystic Fibrosis"):
                    matches = kmp_search(dna_sequence, pattern)
                    if not matches:
                        findings.append(f"⚠️ Possible {name} (missing '{pattern}' in sequence)")
                    else:
                        findings.append(f"{name}: {len(matches)} match(es) of '{pattern}'")
        elif MID:
            input_sequence = load_sequence(dna_sequence)
            reference_sequence = load_sequence(pattern)
            blast_result = blast_search(input_sequence)
            mutations, insertions, deletions = parse_blast_results(blast_results, reference_sequence)
            print("Mutations found:")
            for mutation in mutations:
                print(f"Position: {mutation[0]}, Input: {mutation[1]}, Reference: {mutation[2]}")
            
            print("Insertions found:")
            for insertion in insertions:
                print(f"Position: {insertion[0]}, Base: {insertion[1]}")
            
            print("Deletions found:")
            for deletion in deletions:
                print(f"Position: {deletion[0]}, Base: {deletion[1]}")
            
    if not findings:
        return "✅ No known disorder-related patterns detected."
    else:
        return "\n".join(findings)
    # for disorder in disorders:
    #     name = disorder["name"]
    #     pattern = disorder["pattern"]
    #     min_repeats = disorder["min_repeats"]

    #     # Count pattern repeats
    #     match_positions = find_pattern_positions(dna_sequence, pattern)
    #     count = len(match_positions)

    #     # Always show how many repeats were found
    #     findings.append(f"{name}: {count}x '{pattern}' found.")

    #     # Special case for CF: mutation = deletion of "TTT"
    #     if name.startswith("Cystic Fibrosis"):
    #         if count == 0:
    #             findings.append(f"⚠️ Possible {name} (missing '{pattern}' in sequence)")
        
    #     if name.startswith("Huntington's Disease"):
    #         if count >= min_repeats:
    #             pos_string = ', '.join(str(pos) for pos in match_positions[:10])
    #             more = "..." if count > 10 else ""
    #             findings.append(f"⚠️ Pattern for {name} detected ({count}x '{pattern}') at positions: {pos_string} {more}")

    #     if name.startswith("Fragile X Syndrome"):
    #         if count >= min_repeats:
    #             pos_string = ', '.join(str(pos) for pos in match_positions[:10])
    #             more = "..." if count > 10 else ""
    #             findings.append(f"⚠️ Pattern for {name} detected ({count}x '{pattern}') at positions: {pos_string} {more}")

    #     if name.startswith("Myotonic Dystrophy (Type 1)"):
    #         if count >= min_repeats:
    #             pos_string = ', '.join(str(pos) for pos in match_positions[:10])
    #             more = "..." if count > 10 else ""
    #             findings.append(f"⚠️ Pattern for {name} detected ({count}x '{pattern}') at positions: {pos_string} {more}")
        
    #     if name.startswith("Sickle Cell Anemia"):
    #         if count >= min_repeats:
    #             pos_string = ', '.join(str(pos) for pos in match_positions[:10])
    #             more = "..." if count > 10 else ""
    #             findings.append(f"⚠️ Pattern for {name} detected ({count}x '{pattern}') at positions: {pos_string} {more}")

    #     if name.startswith("Tay-Sachs Disease"):
    #         if count >= min_repeats:
    #             pos_string = ', '.join(str(pos) for pos in match_positions[:10])
    #             more = "..." if count > 10 else ""
    #             findings.append(f"⚠️ Pattern for {name} detected ({count}x '{pattern}') at positions: {pos_string} {more}")

    

# def find_pattern_positions(sequence, pattern):
#     positions = []
#     lps = compute_lps(pattern)
#     i = j = 0

#     while i < len(sequence):
#         if sequence[i] == pattern[j]:
#             i += 1
#             j += 1
#             if j == len(pattern):
#                 positions.append(i - j)
#                 j = lps[j - 1]  # continue search for next match
#         else:
#             if j != 0:
#                 j = lps[j - 1]
#             else:
#                 i += 1

#     return positions
def find_max_consecutive_repeats_kmp(sequence, pattern):
    positions = kmp_search(sequence, pattern)
    if not positions:
        return 0

    pattern_len = len(pattern)
    max_streak = 1
    current_streak = 1

    for i in range(1, len(positions)):
        if positions[i] == positions[i - 1] + pattern_len:
            current_streak += 1
            max_streak = max(max_streak, current_streak)
        else:
            current_streak = 1

    return max_streak



def compute_lps(pattern):
    lps = [0] * len(pattern)
    length = 0  # length of the previous longest prefix suffix
    i = 1

    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1

    return lps


# Example usage
if __name__ == "__main__":
    # Input DNA sequence (you can replace this with file reading logic)
    fasta = input("Input File Name -> ")
    dna_sequence = load_dna_sequence(fasta)
    print(dna_sequence)
    result = check_disorders(dna_sequence)
    print("\n--- SCREENING RESULT ---")
    print(result)


# Example usage
# if __name__ == "__main__":
#     # Load a sample DNA sequence (this can be a path to a .fasta file)
#     fasta = input("Input File Name -> ")
#     dna_sequence = load_dna_sequence(fasta)  # Update with your file path
#     # print(f"Loaded DNA sequence: {dna_sequence}")

#     # Define the pattern to search for
#     pattern = "AGCT"
#     HDpattern = "CAG" # Huntington's Disease if more then 36-40CAG repeats
#     FXSpattern = "CGG" # Fragile X Syndrome if more then 200 CGG repeats
#     MDpattern = "CTG" # Myotonic Dystrophy if more then 50-1500 repeats
#     SCApattern = "GTG" #sickle cell anemia mutation (valine)
#     TSDpattern = "TATC" # Tay-Sachs Disease, 4 base insertion in exon 11 of the HEXA gene
#     # Cystic Fibrosis - delection of 3 bases -> loss of phenylalanine (TTT or TTC)
#     # at position 508 in CFTR gene

#     # Search for the pattern in the sequence
#     matches = kmp_search(dna_sequence, pattern)
#     print(f"Pattern '{pattern}' found at indices: {matches}")
