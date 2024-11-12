import pandas as pd
import ast
import os
import logging

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

from srcGroup6.core import process_sequences
from srcGroup6.read_protein import read_file, parse_proteins

def get_max_proc_number():
    try:
        return os.sysconf('SC_NPROCESSORS_CONF')
    except (AttributeError, ValueError):
        return os.cpu_count()

def main():
    input_dir = "./data"
    output_dir = "./results"
    output_prefix = "group6"

    if output_prefix:
        out_file_path = os.path.join(output_dir, output_prefix)
        if not os.path.isdir(out_file_path):
            os.mkdir(out_file_path)

    # read in the protein sequences
    logging.info("Starting to read and process S288c_proteins file!")

    proteins = read_file(os.path.join(input_dir, "S288c_proteins"))
    # parse the protein sequences
    sequence_df = parse_proteins(proteins)

    # get the system maximum CPU process number
    max_procs = get_max_proc_number()
    print("Maximum system process number:", max_procs)

    # process the sequences
    logging.info("Find 2 and 3 amino acid with high abundance!")
    sequence_df, two_letter_counts, three_letter_counts = process_sequences(sequence_df, num_proc=max_procs)
    # save the processed sequences to a new file

    logging.info("Save problem 1 results to dir {}".format(out_file_path))
    sequence_df.to_csv(os.path.join(out_file_path, "processed_protein_sequences.csv"), index=False)
    # save the two letter counts to a new file
    pd.Series(two_letter_counts).sort_values(ascending=False).to_csv(os.path.join(out_file_path,"two_letter_counts_p1.csv"))
    # save the three letter counts to a new file
    pd.Series(three_letter_counts).sort_values(ascending=False).to_csv(os.path.join(out_file_path,"three_letter_counts_p2.csv"))

    # process the sequences with replacement
    logging.info("Find 2 and 3 amino acid with high abundance after Replacement!")
    sequence_df, two_letter_counts, three_letter_counts, four_letter_counts = (
        process_sequences(sequence_df, replace=True, num_proc=max_procs)
    )
    # save the processed sequences to a new file
    logging.info("Save problem 2 results to dir {}".format(out_file_path))
    sequence_df.to_csv(os.path.join(out_file_path, "processed_protein_sequences_p2.csv"), index=False)
    # save the two letter counts to a new file
    pd.Series(two_letter_counts).sort_values(ascending=False).to_csv(
        os.path.join(out_file_path, "two_letter_counts_p2.csv")
    )
    # save the three letter counts to a new file
    pd.Series(three_letter_counts).sort_values(ascending=False).to_csv(
        os.path.join(out_file_path, "three_letter_counts_p2.csv")
    )
    # save the four letter counts to a new file
    pd.Series(four_letter_counts).sort_values(ascending=False).to_csv(
        os.path.join(out_file_path, "four_letter_counts_p2.csv")
    )

    # TODO: run find_ortholog.py file to get the gossypii_protein_sequences.csv output

    logging.info("Find the same pattern in gossypii!")
    # Read CSV file
    df = pd.read_csv(os.path.join(input_dir, 'gossypii_protein_sequences.csv'))
    # Rename the column
    df = df.rename(columns={"Sequence": "protein_seq"})
    df = df.rename(columns={"Protein Name": "code1"})
    # Filter out rows where 'Sequence' is 'Not Found' and 'Timeout'
    proteins = df[~df['protein_seq'].isin(['Not Found', 'Timeout'])]
    # process the sequences

    sequence_df, two_letter_counts, three_letter_counts = process_sequences(proteins, num_proc=1)
    # save the processed sequences to a new file
    sequence_df.to_csv(os.path.join(out_file_path,"gossypii_processed_protein_sequences.csv"), index=False)
    # save the two letter counts to a new file

    logging.info("Save the 2,3 high abundant amino acid in gossypii to directory {}!".format(out_file_path))
    pd.Series(two_letter_counts).sort_values(ascending=False).to_csv(os.path.join(out_file_path, "gossypii_two_letter_counts.csv"))
    # save the three letter counts to a new file
    pd.Series(three_letter_counts).sort_values(ascending=False).to_csv(os.path.join(out_file_path, "gossypii_three_letter_counts.csv"))

    logging.info("Find the overlapped protein and amino acid!")
    df_processed = pd.read_csv(os.path.join(out_file_path, 'processed_protein_sequences.csv'))
    # Initialize an empty list to store the overlaps
    overlap_records_two = []
    overlap_records_three = []

    # Iterate through each row of the gossypii dataframe
    for index, row in sequence_df.iterrows():
        # Get the Protein Name and overrepresented_two_letter_combos from the current row
        protein_name = row['code1']
        gossypii_overrep_two = row['overrepresented_two_letter_combos']
        gossypii_overrep_three = row['overrepresented_three_letter_combos']

        # Locate the corresponding row in df_processed by Protein Name
        matching_row = df_processed[df_processed['code1'] == protein_name]

        # If a matching row is found, check for overlap
        if not matching_row.empty:
            yeast_overrep_two = ast.literal_eval(matching_row['overrepresented_two_letter_combos'].values[0])
            yeast_overrep_three = ast.literal_eval(matching_row['overrepresented_three_letter_combos'].values[0])

            # Find the overlap between the two lists
            overlap_for_twoletter = set(gossypii_overrep_two).intersection(set(yeast_overrep_two))
            overlap_for_threeletter = set(gossypii_overrep_three).intersection(set(yeast_overrep_three))

            # If there is overlap, record it
            if overlap_for_twoletter:
                overlap_records_two.append([protein_name, list(overlap_for_twoletter)])
            if overlap_for_threeletter:
                overlap_records_three.append([protein_name, list(overlap_for_threeletter)])

    logging.info("Save the results!")
    # Convert the overlap records into a dataframe and save it to a new CSV
    overlap_two_df = pd.DataFrame(overlap_records_two, columns=['Protein Name', 'Overlap'])
    overlap_two_df.to_csv(os.path.join(out_file_path,'overlap_two_letter.csv'), index=False)

    overlap_three_df = pd.DataFrame(overlap_records_three, columns=['Protein Name', 'Overlap'])
    overlap_three_df.to_csv(os.path.join(out_file_path,'overlap_three_letter.csv'), index=False)


if __name__ == "__main__":
    main()