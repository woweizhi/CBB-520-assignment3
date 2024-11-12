import os
from collections import defaultdict

def count_hydrogen_bonds(file_path):
    # Initialize a dictionary to count hydrogen bonds
    hbonds_count = defaultdict(int)

    # Open the PDB file and process line by line
    with open(file_path, "r", encoding="latin-1") as file:
        for line in file:
            # Count number of hydrogen bonds between each pair
            if line.startswith(("ACC", "DNR")):
                parts = line.split()
                amino_acid_1 = parts[1]  # e.g., 'GLU'
                amino_acid_2 = parts[6]  # e.g., 'ASN'

                # Create a sorted tuple (order doesn't matter)
                pair = tuple(sorted([amino_acid_1, amino_acid_2]))

                # Increment the count for this pair
                hbonds_count[pair] += 1

    return hbonds_count

def tally_hydrogen_bonds(folder_path):
    # Initialize a global dictionary to tally all counts
    total_hbonds_count = defaultdict(int)

    # Loop through all files in the specified folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".out"):  # Adjust the extension as necessary
            file_path = os.path.join(folder_path, filename)
            # Count hydrogen bonds in the current file
            hbonds_count = count_hydrogen_bonds(file_path)

            # Tally up the counts
            for pair, count in hbonds_count.items():
                total_hbonds_count[pair] += count

    return total_hbonds_count

def count_amino_acids(file_path):
    # Initialize a dictionary to count amino acids
    amino_acid_count = defaultdict(int)

    # Define a mapping from single-letter to three-letter amino acid codes
    amino_acid_mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    # Open the file and process line by line
    with open(file_path, "r", encoding="latin-1") as file:
        for line in file:
            if line.startswith("SEQ"):
                parts = line.split()
                # The sequence is usually in the third part of the line
                sequence = parts[2]

                # Count each amino acid in the sequence
                for amino_acid in sequence:
                    if amino_acid in amino_acid_mapping:
                        # Get the three-letter code and increment the count
                        three_letter_code = amino_acid_mapping[amino_acid]
                        amino_acid_count[three_letter_code] += 1

    return amino_acid_count

def tally_amino_acids(folder_path):
    # Initialize a global dictionary to tally all counts
    total_amino_acids = defaultdict(int)

    # Loop through all files in the specified folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".out"):  # Adjust the extension as necessary
            file_path = os.path.join(folder_path, filename)
            # Count amino acids in the current file
            amino_acid_count = count_amino_acids(file_path)

            # Tally up the counts
            for amino_acid, count in amino_acid_count.items():
                total_amino_acids[amino_acid] += count

    return total_amino_acids