import pandas as pd
import ast
import os
import logging
import collections
import numpy as np
from collections import Counter
from collections import defaultdict
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq

# Mapping for amino acids
aa_map = {
    'V': 'B', 'L': 'B', 'I': 'B', 'A': 'B', 'M': 'B',
    'F': 'J', 'W': 'J', 'Y': 'J',
    'D': 'O', 'E': 'O',
    'R': 'U', 'H': 'U', 'K': 'U',
    'N': 'Z', 'C': 'Z', 'Q': 'Z', 'S': 'Z', 'T': 'Z'
}

def get_fasta_dictionary(fasta_path):

    yeast_aa = {i.name:str(i.seq) for i in SeqIO.parse(fasta_path, "fasta")}
    print("calculating aa frequencies for %s proteins"%len(yeast_aa))

    return yeast_aa

def get_aminoacid_frequencies(aminoacid_dictionary):

  total_aa_count = Counter()
  for seq in aminoacid_dictionary.values():
      total_aa_count.update(seq) # Add counts of all amino acids found in each sequence
  print("total_aa_count:",total_aa_count)

  # Remove the stop codon from the amino acid frequency distribution
  del total_aa_count['*']

  # Sum the total counts of all amino acids to get the entire length of all sequences combined
  total_length = sum(total_aa_count.values())

  # Calculate the frequency of each amino acid by dividing its count by the total length of amino acids
  aa_frequencies = {aa: count / total_length for aa, count in total_aa_count.items()}

  print("Amino Acid Frequencies:", aa_frequencies)

  return  collections.OrderedDict(aa_frequencies), aminoacid_dictionary

### Generate random proteins
# Assumption that expected protein length ~ exp(observed protein length)

def get_random_protein(length_aa,aa_frequencies):

  result = np.random.choice( len(aa_frequencies) , length_aa  , p= list(aa_frequencies.values()) )
  result = ''.join([ list(aa_frequencies.keys())[a] for a in result])

  return result

# Tabulate the repeats founds in the sampled proteins
# Dictionary that is the count of the repeated amino acids motifs

def find_triplet_motifs(protein_sequence,range_list=[1,2,3]):

  repeat_dictionary = defaultdict( lambda : 0 )

  for len_motif in range_list:

    for i in range(len(protein_sequence) - len_motif*3 ):
        if (protein_sequence[i+len_motif*0:i+len_motif*1] == protein_sequence[i+len_motif*1:i+len_motif*2]) and \
           (protein_sequence[i+len_motif*1:i+len_motif*2] == protein_sequence[i+len_motif*2:i+len_motif*3]) :
           repeat_dictionary[protein_sequence[i:i+len_motif]] += 1

  return dict(repeat_dictionary)


def replace_aas(sequence):
    replaced_sequence = []

    # Iterate over each amino acid in the sequence
    for aa in sequence:
        replaced_sequence.append(aa_map.get(aa, aa))  # Append to the list

    return ''.join(replaced_sequence)

def find_motifs(aminoacid,length=4):
  substring_dict = dict(Counter([aminoacid[i:i+length] for i in range(len(aminoacid))]))
  #repeat_dict = {i:substring_dict[i] for i in substring_dict if i[0]*length == i}
  return substring_dict