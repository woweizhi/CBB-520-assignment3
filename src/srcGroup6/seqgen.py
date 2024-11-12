"""
Code to generate random amino acid sequences.

Amey Chaware
"""

import random


def generate_random_sequence(
    length: int, amino_acids: list[str], probs: list[float]
) -> str:
    """Generates a random amino acid sequence of a given length using a given list of
    amino acids and their probabilities using monte carlo simulation.

    NOT TESTED because of the random nature of the function.

    Args:
        length: length of the amino acid sequence
        amino_acids: list of amino acids
        probs: list of probabilities of each amino acid

    Returns:
        random amino acid sequence
    """
    sequence = random.choices(amino_acids, probs, k=length)

    return "".join(sequence)


def find_amino_acid_distribution(sequence: str) -> tuple[list[str], list[float]]:
    """Finds the probability distribution of amino acids in a given sequence.

    Args:
        sequence: amino acid sequence

    Returns:
        list of amino acids and a corresponding list of their probabilities
    """
    amino_acids = list(set(sequence))
    probs = [sequence.count(aa) / len(sequence) for aa in amino_acids]

    return amino_acids, probs
