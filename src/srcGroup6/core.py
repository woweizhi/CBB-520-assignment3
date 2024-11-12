from functools import partial
from multiprocessing import Pool
from tqdm import tqdm

import pandas as pd

from srcGroup6.replaceaa import replace_amino
from srcGroup6.seqfind import (
    find_four_letter_combos,
    find_sequence_occurance_in_amino_acids,
    find_three_letter_combos,
    find_two_letter_combos,
)
from srcGroup6.seqgen import find_amino_acid_distribution, generate_random_sequence


def process_one_sequence(sequence: str, replace: bool = False) -> list[str]:
    """Processes a single amino acid sequence to find overrepresented two and three
    letter combos which are returned as a list of strings.

    Args:
        sequence: amino acid sequence
        replace: whether to replace specific amino acids with new labels

    Returns:
        list of strings containing the overrepresented two and three letter combos
    """
    if replace:
        sequence, _ = replace_amino(sequence)

    # find two and three letter combos in the sequence
    two_letter_combos = find_two_letter_combos(sequence)
    three_letter_combos = find_three_letter_combos(sequence)

    # randomize the sequence 100 times and count how many times each of the two and
    # three letter combos appear in the randomized sequences
    distribution = find_amino_acid_distribution(sequence)

    # generate 100 random sequences
    random_sequences = [
        generate_random_sequence(len(sequence), *distribution) for _ in range(100)
    ]
    # initialize lists to store the counts of two and three letter combos
    two_letter_counts = []
    three_letter_counts = []
    if replace:
        # in problem 2, we also find four letter combos
        four_letter_combos = find_four_letter_combos(sequence)
        four_letter_counts = []

    # count the number of times each two and three letter combo appears in the random
    # sequences
    for random_sequence in random_sequences:
        two_letter_counts.append(
            find_sequence_occurance_in_amino_acids(two_letter_combos, random_sequence)
        )

        three_letter_counts.append(
            find_sequence_occurance_in_amino_acids(three_letter_combos, random_sequence)
        )
        if replace:
            four_letter_counts.append(
                find_sequence_occurance_in_amino_acids(
                    four_letter_combos, random_sequence
                )
            )

    # convert the counts dictionaries to a dataframe where each row is a combo and each
    # column is a count in a random sequence
    two_letter_df = pd.DataFrame(two_letter_counts)
    three_letter_df = pd.DataFrame(three_letter_counts)

    # find mean and std for each two and three letter combo
    two_letter_mean = two_letter_df.mean()
    two_letter_std = two_letter_df.std()

    three_letter_mean = three_letter_df.mean()
    three_letter_std = three_letter_df.std()

    # compare the counts of two and three letter combos in the original sequence to the
    # mean and std of the counts in the random sequences
    # original sequence counts are in the dictionaries two_letter_combos and three_letter_combos
    overrepresented_two_letter_combos = [
        combo
        for combo, count in two_letter_combos.items()
        if count > two_letter_mean[combo] + 3 * two_letter_std[combo]
    ]
    overrepresented_three_letter_combos = [
        combo
        for combo, count in three_letter_combos.items()
        if count > three_letter_mean[combo] + 3 * three_letter_std[combo]
    ]

    # for problem 2, we also return the overrepresented four letter combos
    if replace:
        four_letter_df = pd.DataFrame(four_letter_counts)

        four_letter_mean = four_letter_df.mean()
        four_letter_std = four_letter_df.std()

        overrepresented_four_letter_combos = [
            combo
            for combo, count in four_letter_combos.items()
            if count > four_letter_mean[combo] + 3 * four_letter_std[combo]
        ]

        return [
            overrepresented_two_letter_combos,
            overrepresented_three_letter_combos,
            overrepresented_four_letter_combos,
        ]

    return [overrepresented_two_letter_combos, overrepresented_three_letter_combos]


def process_sequences(
    sequence_df: pd.DataFrame, replace: bool = False, num_proc: int = 1
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int]]:
    """Processes a dataframe containing amino acid sequences to find overrepresented two
    and three letter combos in each sequence.

    Args:
        sequence_df: dataframe containing amino acid sequences
        replace: whether to replace specific amino acids with new labels
        num_proc: number of processes to use for parallel processing

    Returns:
        tuple containing the input dataframe with two additional columns containing the
        overrepresented two and three letter combos, a dictionary containing the counts
        of each overrepresented two letter combo, and a dictionary containing the counts
        of each overrepresented three letter
    """
    # process each sequence in parallel
    proc_fn = partial(process_one_sequence, replace=replace)
    with Pool(num_proc) as pool:
        #results = pool.map(proc_fn, sequence_df["protein_seq"])
        # TODO: add progress bar for multi-processing
        results = list(tqdm(pool.imap_unordered(proc_fn, sequence_df["protein_seq"]), total=len(sequence_df["protein_seq"])))

    # add column for overrepresented two and three letter combos
    sequence_df["overrepresented_two_letter_combos"] = [result[0] for result in results]
    sequence_df["overrepresented_three_letter_combos"] = [
        result[1] for result in results
    ]

    # for problem 2, we also find overrepresented four letter combos
    if replace:
        sequence_df["overrepresented_four_letter_combos"] = [
            result[2] for result in results
        ]

    # find how many times each combo has been overrepresented
    overrepresented_two_letter_combos = {}
    overrepresented_three_letter_combos = {}

    # for problem 2, we also find overrepresented four letter combos
    if replace:
        overrepresented_four_letter_combos = {}

    for result in results:
        for combo in result[0]:
            overrepresented_two_letter_combos[combo] = (
                overrepresented_two_letter_combos.get(combo, 0) + 1
            )
        for combo in result[1]:
            overrepresented_three_letter_combos[combo] = (
                overrepresented_three_letter_combos.get(combo, 0) + 1
            )
        # for problem 2, we also find overrepresented four letter combos
        if replace:
            for combo in result[2]:
                overrepresented_four_letter_combos[combo] = (
                    overrepresented_four_letter_combos.get(combo, 0) + 1
                )

    if replace:
        return (
            sequence_df,
            overrepresented_two_letter_combos,
            overrepresented_three_letter_combos,
            overrepresented_four_letter_combos,
        )
    
    return (
        sequence_df,
        overrepresented_two_letter_combos,
        overrepresented_three_letter_combos,
    )
