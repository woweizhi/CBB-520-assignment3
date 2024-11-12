# CBB520 Assignment 3 - Group 1 
This repository is a code wrapper for CBB520 assignment-2 project 6-10 in protein pattern discovery.

## Pre-requisite: python packages to run the code
* Before running the code, please have the following python packages installed on your machine: 
```biopython, matplotlib, numpy, pandas, scipy, tqdm```.
  
* For project 8 & 9, we also need stride installation, please follow the [instruction](https://webclu.bio.wzw.tum.de/stride/install.html) to install.

## Code & Folder structure
The **data/** folder contains a list of input files for the projects: `S288c_proteins` and `gossypii_protein_sequences.csv` are the input and metadata 
for project 6; `S288c_proteins`, `Ashbya_gossypii_proteome.faa.fasta` and `ashbya_Sc_orthologs` are the input and metadata for project 7;
`UP000002311_559292_YEAST_v4` folder is the [Alphafold protein structures from S. cerevisiae](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/).
The pdb files need to be unzipped with command: `gzip -d *.gz`.
`stride_output` folder is the output folder after stride processed. `test` folder is the testing input data for project 10. 
```
data/
├── Ashbya_gossypii_proteome.faa.fasta
├── S288c_proteins
├── UP000002311_559292_YEAST_v4
├── ashbya_Sc_orthologs
├── gossypii_protein_sequences.csv
├── stride_output
└── test
```

The **src/** folder wraps up the original source code and util functions for each group. 

```
src/
├── protein1D_pattern_group6.py
├── protein1D_pattern_group7.py
├── protein2D_pattern_group8.py
├── protein2D_pattern_group9.py
├── protein3D_pattern_group10.py
├── srcGroup10
│   └── util.py
├── srcGroup6
│   ├── core.py
│   ├── find_ortholog.py
│   ├── read_protein.py
│   ├── replaceaa.py
│   ├── seqfind.py
│   └── seqgen.py
├── srcGroup7
│   └── util.py
├── srcGroup8
│   └── util.py
└── srcGroup9
    └── stride_preprocess.py
```

The sub-directory in **results/** folder contains the results for each project.

```
results/
├── group10
├── group6
├── group7
├── group8
└── group9
```

## Running
To run the code for all the project: 
```shell
bash run.sh
```
which runs the source code for each project and also contains the [stride processing](stride.sh). 
Notice that you need to change the stride binary file directory to your own binary file address in your machine. 
Here we use `/Users/mac/Downloads/stride/stride` as an example in `stride.sh`.

