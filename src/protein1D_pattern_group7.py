import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from srcGroup7.util import *
from Bio import SeqIO
import os
import logging
from Bio.Seq import Seq

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

def main():
    input_dir = "./data"
    output_dir = "./results"
    output_prefix = "group7"

    if output_prefix:
        out_file_path = os.path.join(output_dir, output_prefix)
        if not os.path.isdir(out_file_path):
            os.mkdir(out_file_path)
    logging.info("Task A: Starting with the S288c_proteins, find [A]3,[AB]3,[ABC]*3 motifs that are 3 orders of magnitude differentially abundant in the protein sequence than in 100 randomized protein sequences. The randomization should maintain the frequency of each amino acid in each protein, but we will not maintain diamino acid frequencies.")

    logging.info("Calculate the observed amino acid frequencies")
    yeast_aa = get_fasta_dictionary(os.path.join(input_dir, "S288c_proteins"))

    aa_frequencies, yeast_aa = get_aminoacid_frequencies(yeast_aa)


    logging.info("Determine the protein length distribution")
    aa_lens = pd.Series([len(aa) for aa in yeast_aa.values()])
    observed_mean = np.mean(aa_lens)

    ax = aa_lens.plot(kind='kde')
    ax.set_title("Observed Protein Length Distribution (n=%s)"%(len(yeast_aa)))
    ax.set_xlabel("Length (amino acids)")
    ax.set_ylabel("Frequency")

    print('observed_mean:',observed_mean)

    ax.axvline(observed_mean,color='r')

    # Save the ax object as a file (use plt.gcf() to get the current figure)
    fig = ax.get_figure()
    fig.savefig(os.path.join(out_file_path, 'protein_len_dist.png'))  # Save the figure
    #plt.show()

    # sanity check that the expected geometric distribution mimics observed exponential distribution at large n
    samples = np.random.default_rng().exponential(scale=observed_mean, size=10000)
    ax = pd.Series(samples).plot(kind='hist')

    # Save the ax object as a file (use plt.gcf() to get the current figure)
    fig = ax.get_figure()
    fig.savefig(os.path.join(out_file_path, 'sanity_check.png'))  # Save the figure
    #plt.show()

    logging.info("Generate random proteins")
    # Should the number of sample proteins be 100 or the number of yeast proteins?
    #length_sample = list(map(int,np.random.default_rng().exponential(scale= observed_mean , size= 100 )))
    length_sample = list(map(int,np.random.default_rng().exponential(scale= observed_mean , size= len(yeast_aa) )))
    random_proteins = [get_random_protein(random_length, aa_frequencies) for random_length in length_sample]
    len(random_proteins)

    logging.info("Tabulate triplet motifs")
    ### Tabulate triplet motifs
    # Merge the triplet motif dictionaries for each random protein
    merge_dictionaries = lambda dict_list : pd.DataFrame(dict_list).agg("sum").to_dict()

    yeast_random_triplet_motifs = merge_dictionaries([find_triplet_motifs(aa) for aa in random_proteins])
    yeast_observed_triplet_motifs = merge_dictionaries([find_triplet_motifs(aa) for aa in yeast_aa.values()])

    print(yeast_random_triplet_motifs)
    print(yeast_observed_triplet_motifs)

    logging.info("Determine overabundant triplet motifs")
    ### Determine overabundant triplet motifs
    #We chose an abundance threshold of 10x above random frequency rather than calculating 3x above the  sample motif abundance standard deviation. The sample motif abundance STD calculation requires normalization of the motif frequency, but we did not dertermine what this factor should be.

    abundance_threshold = 10

    ctr = 0

    yeast_differential_motifs = set()

    for motif in yeast_observed_triplet_motifs:

        random_expected_freq = yeast_random_triplet_motifs.get(motif, 1)

        if yeast_observed_triplet_motifs[motif]  >= random_expected_freq * abundance_threshold :
          print("[%s] %dx over-abundant"%(motif, yeast_observed_triplet_motifs[motif] / random_expected_freq ))
          yeast_differential_motifs.add(motif)

        if yeast_observed_triplet_motifs[motif]  <= random_expected_freq / abundance_threshold:
          print("[%s] %.2gx under-abundant"%(motif, yeast_observed_triplet_motifs[motif] / random_expected_freq  ))
          yeast_differential_motifs.add(motif)

    print('%s differentially expressed triplet aminoacid motifs at %sx abundance threshold '%(len(yeast_differential_motifs), abundance_threshold))

    # Find proteins in yeast proteome that contains the triplet motif

    yeast_genes_different_motifs = set()
    for sequence in yeast_aa:
      for motif in yeast_differential_motifs:
        search_result = yeast_aa[sequence].find(motif * 3 )
        if search_result > 0: yeast_genes_different_motifs.add(sequence)

    print("%s genes with a differential abundant amino acid motif in yeast"%len(yeast_genes_different_motifs))

    # Our result concord with published data:
    #
    #
    # > Sen A, Hsieh WC, Aguilar RC. The Information Content of Glutamine-Rich Sequences Define Protein Functional Characteristics. Proc IEEE Inst Electr Electron Eng. 2017 Feb;105(2):385-393. doi: [10.1109/JPROC.2016.2613076](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7505158/). Epub 2016 Dec 1. PMID: 32963411; PMCID: PMC7505158.
    #
    # > One interesting example of message coding is presented by proteins displaying sequence regions highly-enriched in the **AA** glutamine (or Q, in one-letter code). These proteins play a central role in several neurodegenerative disorders such as Huntington disease and spinocerebellar ataxias [3]. [...] These regions, collectively known as coiled coils (CC), introduce a 3-dimensional factor that is crucial to the interpretation of the information coded in the QR message.

    ## Task B
    #Replacing hydrophobic amino acids Valine, Leucine, Isoleucine, Alanine, Methionine, with B, and replacing the aromatic amino acids Phenylalanine, tryptophan, and Tyrosine with J, and replacing the acidic amino acids Aspartic Acid and Glutamic acid with O, and replacing the basic amino acids Arginine, Histidine, and Lysine with U, and replacing the amino acids with polar neutral side chains, Asparagine, Cysteine, Glutamine, Serine, and Threonine with Z, repeat the analysis in 1), though in this case look for 2 and 3 and 4 amino acids or amino acid groups that are 3 orders of magnitude more abundant in the protein sequence than in 100 randomized protein sequences.
    logging.info("Task B: Replacing hydrophobic amino acids Valine, Leucine, Isoleucine, Alanine, Methionine, with B, and replacing the aromatic amino acids Phenylalanine, tryptophan, and Tyrosine with J, and replacing the acidic amino acids Aspartic Acid and Glutamic acid with O, and replacing the basic amino acids Arginine, Histidine, and Lysine with U, and replacing the amino acids with polar neutral side chains, Asparagine, Cysteine, Glutamine, Serine, and Threonine with Z, repeat the analysis in 1), though in this case look for 2 and 3 and 4 amino acids or amino acid groups that are 3 orders of magnitude more abundant in the protein sequence than in 100 randomized protein sequences.")

    replaced_proteins = { protein:replace_aas(yeast_aa[protein]) for protein in yeast_aa}

    new_distribution, replaced_proteins = get_aminoacid_frequencies(replaced_proteins)

    length_sample = list(map(int,np.random.default_rng().exponential(scale= observed_mean , size= len(replaced_proteins) )))
    random_proteins = [get_random_protein(random_length, new_distribution) for random_length in length_sample]
    len(random_proteins)

    random_repeat_aa_groups = merge_dictionaries([find_triplet_motifs(aa) for aa in random_proteins])
    observed_repeat_aa_groups = merge_dictionaries([find_triplet_motifs(aa) for aa in replaced_proteins.values()])

    abundance_threshold = 3

    ctr = 0

    for motif in observed_repeat_aa_groups:

        random_expected_freq = random_repeat_aa_groups.get(motif, 1)

        if observed_repeat_aa_groups[motif]  >= random_expected_freq * abundance_threshold :
          ctr += 1
          print("[%s] %dx over-abundant"%(motif, observed_repeat_aa_groups[motif] / random_expected_freq ))
        if observed_repeat_aa_groups[motif]  <= random_expected_freq / abundance_threshold:
          ctr += 1
          print("[%s] %.2gx under-abundant"%(motif, observed_repeat_aa_groups[motif] / random_expected_freq  ))

    print('%s differentially expressed triplet aminoacid motifs'%(ctr))

    # We need to modify the motif finder, since we are looking for whether {2,3,4}-mers of aa groups are more abundant
    # Find all substrings of length

    random_repeat_aa_groups = merge_dictionaries([ merge_dictionaries([ find_motifs(aa,2),find_motifs(aa,3),find_motifs(aa,4)]) for aa in random_proteins])
    observed_repeat_aa_groups = merge_dictionaries([ merge_dictionaries([ find_motifs(aa,2),find_motifs(aa,3),find_motifs(aa,4)]) for aa in replaced_proteins.values()])

    abundance_threshold = 3

    ctr = 0

    for motif in observed_repeat_aa_groups:

        random_expected_freq = random_repeat_aa_groups.get(motif, 1)

        if observed_repeat_aa_groups[motif]  >= random_expected_freq * abundance_threshold :
          ctr += 1
          print("[%s] %dx over-abundant"%(motif, observed_repeat_aa_groups[motif] / random_expected_freq ))
        if observed_repeat_aa_groups[motif]  <= random_expected_freq / abundance_threshold:
          ctr += 1
          print("[%s] %.2gx under-abundant"%(motif, observed_repeat_aa_groups[motif] / random_expected_freq  ))

    print('%s differentially expressed {2,3,4}-mer aminoacid group motifs'%(ctr))

    # Find proteins in yeast proteome that contains the triplet motif

    motifs_to_find = ['SSS','QQQ','NNN','DDD']

    has_motif = set()
    for sequence in yeast_aa:
      for motif in motifs_to_find:
        search_result = yeast_aa[sequence].find(motif * 3 )
        if search_result > 0: has_motif.add(sequence)

    print(has_motif)

    logging.info("Task C: Protein sequences tend to be much more conserved than DNA sequences. Thus for proteins found with significantly overrepresented patterns, are those overrepresented, by the same measure, in the ortholog in gossypii, if there is one?")
    ## Task C

    # Protein sequences tend to be much more conserved than DNA sequences. Thus for proteins found with significantly overrepresented patterns, are those overrepresented, by the same measure, in the ortholog in gossypii, if there is one?
    #
    # Group 7
    # 	•	As for group 6, but instead of looking for short sequences, look for direct repeats of at least three copies of 1,2,3 amino acids, such as GHTGHTGHT and answer questions
    # 	1) (Amino acids)
    # 	2) (combining similar amino acids together)
    # 	3) (is the protein similar in A. gossypii?) for those repeats
    #
    # Ortholog mapping from Dietrich 2004: https://www.science.org/doi/10.1126/science.1095781
    #
    #
    # Yuxin: Below I am adding more to Task C. What I did: 1) upload the ashbya file that Jason got from Fred, 2) extract the sequence from the file like what we did with the yeast, 3) screen for the triplet amino acid motifs, using the find_triplet_motif, to detect the overrepresented motifs in ashbya, 4) and then compare the frequencies between yeast and ashbya.
    #
    # Instead of taking the overrepresented motifs of the yeast we obtained from Task B and searching those in the ashbya proteome, I determined the triplet motifs in ashbya and then compared the frequencies between the two species.

    # Load ashbya_Sc_orthologs

    ashbya_path = os.path.join(input_dir, 'Ashbya_gossypii_proteome.faa.fasta')
    #ashbya_path = './data/UP000000591.fasta'

    # Modify regular expression inorder to extract Ashbya protein names
    # Note not perfect as there is a common gene name as well gene id
    ashbya_aa = {i.description[i.description.find('GN'):].split(' ')[0].split('=')[1]:str(i.seq) for i in SeqIO.parse(ashbya_path, "fasta")}
    print("calculating aa frequencies for %s proteins"%len(ashbya_aa))

    # Generate random proteins based off ashyba aa distribution
    new_distribution, ashbya_aa = get_aminoacid_frequencies(ashbya_aa)

    merge_dictionaries = lambda dict_list : pd.DataFrame(dict_list).agg("sum").to_dict()
    #%%

    aa_lens = pd.Series([len(aa) for aa in ashbya_aa.values()])
    observed_mean = np.mean(aa_lens)

    length_sample = list(map(int,np.random.default_rng().exponential(scale= observed_mean , size= len(ashbya_aa) )))
    random_proteins = [get_random_protein(random_length, aa_frequencies) for random_length in length_sample]
    len(random_proteins)

    #%%

    ashbya_random_triplet_motifs = merge_dictionaries([find_triplet_motifs(aa) for aa in random_proteins])
    ashbya_observed_triplet_motifs = merge_dictionaries([find_triplet_motifs(aa) for aa in ashbya_aa.values()])

    #%%

    # Determine overrepresentation in Ashbya
    abundance_threshold = 10
    ashbya_differential_motifs = set()

    for motif in ashbya_observed_triplet_motifs:
        ashbya_expected_freq = ashbya_random_triplet_motifs.get(motif, 1)

        if ashbya_observed_triplet_motifs[motif] >= ashbya_expected_freq * abundance_threshold:
            print("[%s] %dx over-abundant" % (motif, ashbya_observed_triplet_motifs[motif] / ashbya_expected_freq))
            ashbya_differential_motifs.add(motif)
        if ashbya_observed_triplet_motifs[motif] <= ashbya_expected_freq / abundance_threshold:
            print("[%s] %.2gx under-abundant" % (motif, ashbya_observed_triplet_motifs[motif] / ashbya_expected_freq))
            ashbya_differential_motifs.add(motif)

    print('%s differentially expressed triplet aminoacid motifs in ashbya' % len(ashbya_differential_motifs))

    #%%

    # Find proteins in ashbya proteome that contains the triplet motif

    ashbya_genes_different_motifs = set()
    for sequence in ashbya_aa:
      for motif in ashbya_differential_motifs:
        search_result = ashbya_aa[sequence].find(motif * 3 )
        if search_result > 0: ashbya_genes_different_motifs.add(sequence)

    print("%s genes with a differential abundant amino acid motif in ashbya"%len(ashbya_genes_different_motifs))


    # Check if orthologs were abundant in each analysis
    # ortholog_mapping = ! cat ashbya_Sc_orthologs
    # ortholog_mapping =  {mapping.split(' ')[1]:mapping.split(' ')[7:] for mapping in ortholog_mapping}
    # len(ortholog_mapping)

    # Read the content of the file directly using Python file handling
    with open(os.path.join(input_dir, 'ashbya_Sc_orthologs'), 'r') as file:
        ortholog_mapping = file.readlines()

    # Create a dictionary where the key is the second word and the value is the seventh onward
    ortholog_mapping = {mapping.split()[1]: mapping.split()[7:] for mapping in ortholog_mapping if len(mapping.split()) > 7}

    # Print the length of the dictionary
    print(len(ortholog_mapping))

    #%%

    ctr = 0
    for gene in ashbya_genes_different_motifs:
      if gene in ortholog_mapping:
        for c in ortholog_mapping[gene]:
          if c in yeast_genes_different_motifs:
            print(gene,c)
            ctr += 1

    #%%

    print("%s ashbya proteins with differentially abundant aa repeat motifs had yeast orthologs with same property"%ctr)

    #%%

    print(len(ashbya_genes_different_motifs),len(yeast_genes_different_motifs))

if __name__ == "__main__":
    main()
