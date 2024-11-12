"""
Code to replace amino acid with either B, J, O or U sequences.

Valine, Leucine, Isoleucine, Alanine, Methionine, with B
Phenylalanine, tryptophan, and Tyrosine with J
Aspartic Acid and Glutamic acid with O
Arginine, Histidine, and Lysine with U
Asparagine, Cysteine, Glutamine, Serine, and Threonine with Z

Sarah hui
"""

def replace_amino(protein_sequence: str) -> tuple[str, dict[str, int]]:
    """
    Replaces specific amino acids in a protein sequence with new labels
    and returns the modified sequence along with a dictionary of replacement counts.

    Replacements:
        - Valine, Leucine, Isoleucine, Alanine, Methionine -> B
        - Phenylalanine, Tryptophan, Tyrosine -> J
        - Aspartic Acid, Glutamic Acid -> O
        - Arginine, Histidine, Lysine -> U
        - Asparagine, Cysteine, Glutamine, Serine, Threonine -> Z

    Args:
        protein_sequence: string representing a sequence of amino acids

    Returns:
        tuple: (modified sequence, dict of replacement counts)
    """

    # Define the mappings for amino acids
    replacements = {
        'B': ['V', 'L', 'I', 'A', 'M'],  # Valine, Leucine, Isoleucine, Alanine, Methionine
        'J': ['F', 'W', 'Y'],            # Phenylalanine, Tryptophan, Tyrosine
        'O': ['D', 'E'],                 # Aspartic Acid, Glutamic Acid
        'U': ['R', 'H', 'K'],            # Arginine, Histidine, Lysine
        'Z': ['N', 'C', 'Q', 'S', 'T'],  # Asparagine, Cysteine, Glutamine, Serine, Threonine
    }

    # Create a dictionary to track the count of each replacement
    replacement_counts = {'B': 0, 'J': 0, 'O': 0, 'U': 0, 'Z': 0}

    # Convert the protein sequence to a list for easier replacement
    sequence_list = list(protein_sequence)

    # Iterate over the sequence and apply replacements
    for i, amino_acid in enumerate(sequence_list):
        for replacement, amino_acids in replacements.items():
            if amino_acid in amino_acids:
                sequence_list[i] = replacement
                replacement_counts[replacement] += 1

    # Join the list back into a string
    modified_sequence = ''.join(sequence_list)

    # Return both the modified sequence and the dictionary with the count of each replacement
    return modified_sequence, replacement_counts
