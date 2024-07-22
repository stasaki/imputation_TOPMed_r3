import pandas as pd
import os
import shutil
import argparse
import glob

def backup_file(file_path):
    backup_path = f"{file_path}.bak"
    shutil.copyfile(file_path, backup_path)
    print(f"Backup created: {backup_path}")

def rename_variants(bim_file, problematic_variants, variant_counts):
    # Create a backup of the original .bim file
    backup_file(bim_file)

    # Read the .bim file
    bim_df = pd.read_csv(bim_file, sep='\t', header=None)

    # Ensure the .bim file has 6 columns
    if bim_df.shape[1] != 6:
        raise ValueError(f"File {bim_file} does not have 6 columns. Check the format.")

    # Rename variants that are in the problematic variants list
    for idx, row in bim_df.iterrows():
        variant_id = row[1]
        if variant_id in variant_counts:
            variant_counts[variant_id] += 1
            new_variant_id = f"{variant_id}-{variant_counts[variant_id]}"
            bim_df.at[idx, 1] = new_variant_id

    # Save the modified .bim file
    bim_df.to_csv(bim_file, sep='\t', header=False, index=False)
    print(f"Renamed variants in: {bim_file}")

def process_files(file_paths, problematic_variants):
    # Dictionary to keep track of counts for each problematic variant
    variant_counts = {variant: 0 for variant in problematic_variants}

    for bim_file in file_paths:
        if os.path.exists(bim_file):
            rename_variants(bim_file, problematic_variants, variant_counts)
        else:
            print(f"File {bim_file} does not exist. Skipping.")

def main():
    parser = argparse.ArgumentParser(description="Rename variants in PLINK .bim files to add index for problematic variants and create backups.")
    parser.add_argument('bim_files', nargs='+', help='List of .bim files to process')
    parser.add_argument('--problematic', required=True, help='File containing problematic variants')
    args = parser.parse_args()

    # Read problematic variants from file
    with open(args.problematic, 'r') as f:
        problematic_variants = [line.strip() for line in f.readlines()]

    # Use glob to expand wildcards
    file_paths = []
    for pattern in args.bim_files:
        file_paths.extend(glob.glob(pattern))

    process_files(file_paths, problematic_variants)

if __name__ == "__main__":
    main()
