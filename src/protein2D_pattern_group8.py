from srcGroup8.util import *
import csv
import logging
import pandas as pd
import numpy as np
import scipy.stats as stats

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

def main():
    input_dir = "./data"
    output_dir = "./results"
    output_prefix = "group8"

    if output_prefix:
        out_file_path = os.path.join(output_dir, output_prefix)
        if not os.path.isdir(out_file_path):
            os.mkdir(out_file_path)

    logging.info("Count hydrogen bonds")
    # Count hydrogen bonds
    folder_path = os.path.join(input_dir, "stride_output")
    total_hbonds_count = tally_hydrogen_bonds(folder_path)
    total_amino_acid_counts = tally_amino_acids(folder_path)

    # Print the total results
    print("Total counts of hydrogen bonds for each pair of amino acids:")
    for pair, count in total_hbonds_count.items():
        print(f'{pair}: {count} hydrogen bonds')

    print("Total counts of amino acids across all sequences:")
    for amino_acid, count in total_amino_acid_counts.items():
        print(f'{amino_acid}: {count}')

    # Define the output CSV file name
    output_file = os.path.join(out_file_path, "hydrogen_bonds_counts.csv")

    # Open the CSV file in write mode
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write the header to the CSV file
        writer.writerow(["Amino Acid Pair", "Count of Hydrogen Bonds"])

        # Write each pair and its count to the CSV file
        for pair, count in total_hbonds_count.items():
            writer.writerow([pair, count])

    print(f"Results saved to {output_file}")

    #%%

    # Define the output CSV file name
    output_file = os.path.join(out_file_path, "amino_acids_counts.csv")

    # Open the CSV file in write mode
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write the header for amino acid counts
        writer.writerow(["Amino Acid", "Count"])

        # Write total counts of amino acids across all sequences
        for amino_acid, count in total_amino_acid_counts.items():
            writer.writerow([amino_acid, count])

    print(f"Amino acid counts saved to {output_file}")

    logging.info("Statistical Analysis")
    # Statistical Analysis

    # Data input
    logging.info("Read amino_acids_counts and hydrogen_bonds_counts data input")
    df_AA = pd.read_csv(os.path.join(out_file_path, "amino_acids_counts.csv"))
    df_pairs= pd.read_csv(os.path.join(out_file_path, "hydrogen_bonds_counts.csv"))

    # Compute expectation
    n_AA = np.sum(df_AA['Count'])
    n_bonds = np.sum(df_pairs['Count of Hydrogen Bonds'])
    def expected_n_pairs(AA1, AA2):
        n1 = df_AA.loc[df_AA['Amino Acid']==AA1, 'Count'].squeeze()
        n2 = df_AA.loc[df_AA['Amino Acid']==AA2, 'Count'].squeeze()
        return (2 * (n1/n_AA) * (n2/n_AA) * n_bonds)

    df_pairs['Expected Count'] = 0
    for i in range(df_pairs.shape[0]):
        AA1, AA2 = eval(df_pairs.loc[i, 'Amino Acid Pair'])
        df_pairs.loc[i, 'Expected Count'] = expected_n_pairs(AA1, AA2)

    logging.info("Perform Chi square test!")
    # Chi square test
    def chi_test(observation, expectation):
        chi2 = (observation - expectation) ** 2 / expectation
        p = stats.chi2.sf(chi2, df=1)
        return chi2, p

    chi2, p = chi_test(df_pairs['Count of Hydrogen Bonds'], df_pairs['Expected Count'])
    df_pairs['chi2_stat'] = chi2
    df_pairs['p_value'] = p

    alpha = 0.05
    df_pairs['significance'] = (df_pairs['p_value'] < alpha)

    df_pairs['Comparison'] = 'Not significant'
    df_pairs.loc[(df_pairs['Count of Hydrogen Bonds'] <= df_pairs['Expected Count']) & (df_pairs['significance']),
                 'Comparison'] = 'Less likely'
    df_pairs.loc[(df_pairs['Count of Hydrogen Bonds'] > df_pairs['Expected Count']) & (df_pairs['significance']),
                 'Comparison'] = 'More likely'

    logging.info(df_pairs)

    logging.info("save test results!")
    df_pairs.to_csv(os.path.join(out_file_path, 'Test results.csv'))

if __name__ == "__main__":
    main()

#%% md

# Visualization

## TODO: R code for result visualization
#
# %load_ext rpy2.ipython
#
# #%%
#
# %%R
# knitr::opts_chunk$set(echo = TRUE)
# library(dplyr)
# library(tidyverse)
# library(ggplot2)
#
# #%%
#
# %%R
# df <- read.csv('./results/Test results.csv', header = TRUE)
#
# df1 <- df %>%
#   mutate(AA1 = gsub("\\('([A-Z]+)',.*", "\\1", Amino.Acid.Pair),
#          AA2 = gsub(".*, '([A-Z]+)'\\)", "\\1", Amino.Acid.Pair)) %>%
#   select(AA1, AA2, p_value, significance, Comparison)
#
# #%%
#
# %%R
# amino_acids <- c('ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
#                  'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR')
#
# all_pairs <- expand.grid(AA1 = amino_acids, AA2 = amino_acids, stringsAsFactors = FALSE) %>%
#   filter(AA1 <= AA2)
#
# #Find the missing pair
# missing_pairs <- anti_join(all_pairs, df1, by = c("AA1", "AA2"))
#
# print(missing_pairs)
# # PRO - PRO is missing
#
# new_row <- data.frame(AA1 = "PRO", AA2 = "PRO", p_value = NA, significance = NA, Comparison = NA)
# df2 <- rbind(df1, new_row)
#
# df2_1 <- df2 %>%
#   mutate(Comparison2 = Comparison) %>%
#   mutate(Comparison2 = ifelse(p_value > 0.01, "Not Significant", Comparison2))
#
#
# complete_data <- all_pairs %>%
#   left_join(df2_1, by = c("AA1", "AA2")) %>%
#   arrange(AA1) %>%
#   select(AA1, AA2, Comparison2)
#
# heatmap_data <- complete_data %>%
#   spread(key = AA2, value = Comparison2)
#
#
# p = ggplot(complete_data, aes(x = AA1, y = AA2, fill = Comparison2)) +
#   geom_tile(color = "black") +
#   scale_fill_manual(
#     values = c(
#       "Less likely" = "#377eb8",
#       "More likely" = "#e41a1c",
#       "Not Significant" = "#f7f7f7"
#     ),
#     na.value = "grey80"
#   ) +
#   theme_minimal(base_size = 15) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
#     axis.text.y = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     legend.position = "right",
#     legend.title = element_text(size = 14, face = "bold"),
#     legend.text = element_text(size = 12)
#   ) +
#   labs(
#     title = "Heatmap of Abnormal Pairs",
#     x = "Amino Acid 1",
#     y = "Amino Acid 2",
#     fill = "Comparison"
#   ) +
#   coord_fixed()
#
#
# ggsave("CBBHW.pdf", p, width = 10, height = 10, dpi = 800)
