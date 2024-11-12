import os
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')
import numpy as np
from collections import defaultdict
from scipy.stats import chisquare
import json

stride_output_dir = "./data/stride_output/"
result_dir = "./results/group9"
os.makedirs(result_dir, exist_ok=True)
stride_output_files = os.listdir(stride_output_dir)
structure_counts = defaultdict(lambda: defaultdict(int))
total_counts = defaultdict(int)
structure_distribution = defaultdict(lambda: defaultdict(float))
total_distribution = defaultdict(float)

for stride_file in stride_output_files:
    substructure = {}
    with open(os.path.join(stride_output_dir, stride_file), 'r', encoding="latin-1") as f:  # add encoding to latin-1
    #with open(os.path.join(stride_output_dir, stride_file), 'r') as f:
        for line in f:
            if line.startswith('LOC'):
                tokens = line.split()
                st = tokens[1]
                start = int(tokens[3])
                end = int(tokens[6]) + 1
                for k in range(start, end):
                    substructure[k] = st
            if line.startswith("ASG"):
                tokens = line.split()
                amino_acid = tokens[1]  # Amino acid type
                loc = tokens[3]
                structure = tokens[6]  # Secondary structure type (H, C, etc.)
                if structure == 'Turn':
                    structure = substructure[int(loc)]
                structure_counts[structure][amino_acid] += 1
                total_counts[amino_acid] += 1

for structure, structure_dict in structure_counts.items():
    structure_residues = sum(structure_dict.values())
    for amino_acid, count in structure_dict.items():
        structure_distribution[structure][amino_acid] = count / structure_residues
total_residues = sum(total_counts.values())
for amino_acid, count in total_counts.items():
    total_distribution[amino_acid] = count / total_residues

amino_acid_names = list(total_counts.keys())
amino_acid_counts = list(total_counts.values())
plt.bar(amino_acid_names, amino_acid_counts, alpha=0.6)
plt.xlabel('Amino Acid')
plt.ylabel('Count')
plt.title('Total Amino Acid Counts')
plt.gcf().autofmt_xdate()
plt.savefig(os.path.join(result_dir, 'total_amino_acids_counts.png'))
plt.close()

for structure, structure_dict in structure_counts.items():
    amino_acid_names = list(structure_dict.keys())
    amino_acid_counts = list(structure_dict.values())
    plt.bar(amino_acid_names, amino_acid_counts, alpha=0.6)
    plt.xlabel('Amino Acid')
    plt.ylabel('Count')
    plt.title(f'{structure} Amino Acid Counts')
    plt.gcf().autofmt_xdate()
    plt.savefig(os.path.join(result_dir, f'{structure}_amino_acids_counts.png'))
    plt.close()

amino_acid_names = list(total_distribution.keys())
amino_acid_distribution = list(total_distribution.values())
plt.bar(amino_acid_names, amino_acid_distribution, alpha=0.6)
plt.xlabel('Amino Acid')
plt.ylabel('Distribution')
plt.title('Total Amino Acid Distribution')
plt.gcf().autofmt_xdate()
plt.savefig(os.path.join(result_dir, 'total_amino_acids_distribution.png'))
plt.close()

for structure, structure_dict in structure_distribution.items():
    amino_acid_names = list(structure_dict.keys())
    amino_acid_distribution = list(structure_dict.values())
    plt.bar(amino_acid_names, amino_acid_distribution, alpha=0.6)
    plt.xlabel('Amino Acid')
    plt.ylabel('Distribution')
    plt.title(f'{structure} Amino Acid Distribution')
    plt.gcf().autofmt_xdate()
    plt.savefig(os.path.join(result_dir, f'{structure}_amino_acids_distribution.png'))
    plt.close()

total_alpha_helix_residues = sum(structure_counts['AlphaHelix'].values())
observed_alphahelix = []
expected_alphahelix = []
txt_file = os.path.join(result_dir, 'chi_test_results.txt')
with open(txt_file, 'w') as f:
    for amino_acid in amino_acid_names:
        observed_value = structure_counts['AlphaHelix'][amino_acid]
        expected_value = total_counts[amino_acid] * (total_alpha_helix_residues / total_residues)
        observed_alphahelix.append(observed_value)
        expected_alphahelix.append(expected_value)
        standardized_residuals = (observed_value - expected_value) / np.sqrt(expected_value)

        outstring = f"Amino acid: {amino_acid}, Observed value:{observed_value}, Expected value: {expected_value}, Standardized residual: {standardized_residuals}.\n"
        f.write(outstring)
        print(outstring)
        if standardized_residuals > 1.96:
            outstring = f"{amino_acid} is significantly more likely to be in Alphahelix.\n"
            f.write(outstring)
            print(outstring)
        elif standardized_residuals < -1.96:
            outstring = f"{amino_acid} is significantly less likely to be in Alphahelix.\n"
            f.write(outstring)
            print(outstring)
        else:
            outstring = f"The distribution of {amino_acid} has no significant difference.\n"
            f.write(outstring)
            print(outstring)
    chi2, p_value = chisquare(observed_alphahelix, expected_alphahelix)
    if p_value < 0.05:
        outstring = "The distribution of amino acids in Alphahelix is significantly different.\n"
        f.write(outstring)
        print(outstring)
    else:
        outstring = "The distribution of amino acids in Alphahelix is not significantly different.\n"
        f.write(outstring)
        print(outstring)

width = 0.4
x_pos = np.array(range(len(amino_acid_names)))
plt.bar(x_pos - width / 2, observed_alphahelix, width=width, alpha=0.6, label='Observed')
plt.bar(x_pos + width / 2, expected_alphahelix, width=width, alpha=0.6, label='Expected')
plt.xlabel('Amino Acid')
plt.ylabel('Count')
plt.title('Alpha Helix Amino Acid Distribution')
plt.xticks(x_pos, amino_acid_names)
plt.legend()
plt.gcf().autofmt_xdate()
plt.savefig(os.path.join(result_dir, 'Chi_test_alpha_helix_amino_acids_distribution.png'))
plt.close()

dict_list = [structure_counts, total_counts, structure_distribution, total_distribution]
with open(os.path.join(result_dir, 'data.json'), 'w') as f:
    json.dump(dict_list, f, indent=4)
