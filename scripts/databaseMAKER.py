import csv
import gzip
# import requests

# def download_clinvar_data():
#     url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
#     filename = "variant_summary.txt.gz"
    
#     print(f"Downloading ClinVar data from {url}...")
#     response = requests.get(url)
#     with open(filename, 'wb') as f:
#         f.write(response.content)
#     print("Download complete.")

def extract_pathogenic_variants():
    print("Extracting pathogenic variants...")
    
    pathogenic_variants = []
    with gzip.open('variant_summary.txt.gz', mode='rt') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if "pathogenic" in row["ClinicalSignificance"].lower():
                variant = {
                    "GeneSymbol": row["GeneSymbol"],
                    "Type": row["Type"],
                    "Name": row["Name"],
                    "Chromosome": row["Chromosome"],
                    "Start": row["Start"],
                    "ReferenceAllele": row["ReferenceAllele"],
                    "AlternateAllele": row["AlternateAllele"],
                    "ClinicalSignificance": row["ClinicalSignificance"],
                    "Condition": row.get("Conditions", row.get("Condition(s)", "Unknown"))
                }
                pathogenic_variants.append(variant)
    
    print(f"Found {len(pathogenic_variants)} pathogenic variants.")
    return pathogenic_variants

def save_variants(variants):
    output_file = "pathogenic_variants.csv"
    print(f"Saving to {output_file}...")
    keys = variants[0].keys()
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(variants)
    print("Saved!")

def main():
    # download_clinvar_data()
    variants = extract_pathogenic_variants()
    save_variants(variants)

if __name__ == "__main__":
    main()
