"""
Functions to find sequences in a list of amino acids.
"""

from collections import Counter


def find_two_letter_combos(sequence: str) -> dict[str, int]:
    """Finds all present two amino acid combinations and their frequency in a given
    amino acid sequence.

    Args:
        sequence: amino acid sequence
        repeat_length: length of the repeat

    Returns:
        dictionary with two letter combos as keys and their counts as values
    """

    # zip the offsetted sequence with the original sequence
    combos = ["".join(i) for i in zip(sequence, sequence[1:])]
    repeats = Counter(combos)

    return repeats


def find_three_letter_combos(sequence: str) -> dict[str, int]:
    """Finds all present three amino acid combinations and their frequency in a given
    amino acid sequence.

    Args:
        sequence: amino acid sequence

    Returns:
        dictionary with three letter combos as keys and their counts as values
    """

    # zip the offsetted sequence with the original sequence
    combos = ["".join(i) for i in zip(sequence, sequence[1:], sequence[2:])]
    repeats = Counter(combos)

    return repeats


def find_four_letter_combos(sequence: str) -> dict[str, int]:
    """Finds all present four amino acid combinations and their frequency in a given
    amino acid sequence.

    Args:
        sequence: amino acid sequence

    Returns:
        dictionary with four letter combos as keys and their counts as values
    """

    # zip the offsetted sequence with the original sequence
    combos = [
        "".join(i) for i in zip(sequence, sequence[1:], sequence[2:], sequence[3:])
    ]
    repeats = Counter(combos)

    return repeats


def find_sequence_occurance_in_amino_acids(
    templates: list[str], amino_acids: str
) -> dict[str, int]:
    """Finds the number of times given sequences occur in an amino acid string.

    Args:
        templates: list of sequences to find
        amino_acids: amino acid sequence

    Returns:
        dictionary with sequence as key and count as value
    """

    # initialize dictionary with 0 counts
    counts = {template: 0 for template in templates}

    for template in templates:
        counts[template] = amino_acids.count(template)

    return counts
