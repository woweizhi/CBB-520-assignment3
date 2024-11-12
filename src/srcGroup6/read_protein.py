"""
Functions to read in and parse yeast protein sequence data and return a dataframe

Amey Chaware
"""

import pandas as pd


def read_file(filepath: str) -> list[str]:
    """Reads in a protein sequence file and returns a list with each protein as an
    element.

    It is assumed that each protein line starts with a '>' character and ends with a '*'
    character.

    Args:
        filepath: path to the protein sequence file

    Returns:
        list containing the metadata and amino acid sequence of each protein
    """
    with open(filepath) as f:
        contents = f.read()

    proteins = contents.split("\n>")
    # remove the '>' character from the first element
    proteins[0] = proteins[0].replace(">", "")
    # remove the '*' character from the end of each element
    for i, protein in enumerate(proteins):
        if protein[-1] == "*":
            proteins[i] = protein[:-1]
    # proteins = [protein.split("*")[0] for protein in proteins]
    # in each element replace the new line character with an empty string
    proteins = [protein.replace("\n", "") for protein in proteins]

    return proteins


def parse_proteins(proteins: list[str]) -> pd.DataFrame:
    """Parses a list of protein sequences and returns a DataFrame with the metadata and
    amino acid sequence of each protein.

    Args:
        proteins: list containing the metadata and amino acid sequence of each protein

    Returns:
        DataFrame containing the metadata and amino acid sequence of each protein
    """
    parsed_proteins = []
    for i, protein in enumerate(proteins):
        protein_dict = {}
        # first split using the quotation mark to get the file comment and protein
        # sequence
        metadata, comment, protein_seq = protein.split('"')

        # if metadata contains "Dubious", ignore and move on
        if "Dubious" in metadata:
            continue

        # split the metadata using comma and space
        metadata = metadata.split(", ")[:-1]
        if len(metadata) > 5:
            raise ValueError("Unknown metadata format for protein {}".format(i))

        # get the first element of metadata and split using the space character
        # to get two codes and one ID
        code1, code2, sgd_id = metadata[0].split(" ")
        # add these to the dictionary
        protein_dict["code1"] = code1
        protein_dict["code2"] = code2
        protein_dict["sgd_id"] = sgd_id.split(":")[1]

        # second element of metadata is the location
        protein_dict["location"] = metadata[1]

        # third element of metadata is the release
        protein_dict["release"] = metadata[2]

        # fourth element of metadata is reverse complement, if present
        if metadata[3] == "reverse complement":
            protein_dict["reverse_complement"] = True
            curr_element = 4
        else:
            protein_dict["reverse_complement"] = False
            curr_element = 3

        # if there is a fifth element, it is the type
        protein_dict["type"] = metadata[curr_element]

        if len(metadata) > curr_element + 1:
            raise ValueError("Unknown metadata format for protein {}".format(i))

        # add comment and protein sequence to the dictionary
        protein_dict["comment"] = comment
        protein_dict["protein_seq"] = protein_seq

        parsed_proteins.append(protein_dict)

    parsed_proteins = pd.DataFrame(parsed_proteins)

    return parsed_proteins


if __name__ == "__main__":
    proteins = read_file("S288c_proteins")
    parsed_proteins = parse_proteins(proteins)
    print(parsed_proteins.head())
