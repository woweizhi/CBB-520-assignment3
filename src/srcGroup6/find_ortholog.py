"""
Find the ortholog protein in gossypii by calling uniprot's API
Yuyan Zhang
"""
import pandas as pd
import requests
import time

df = pd.read_csv('processed_protein_sequences.csv')

# get the first column
protein_names = df['code1'].tolist()

# UniProt API URL
uniprot_url = "https://rest.uniprot.org/uniprotkb/search"

#Save a dictionary of the queried protein sequences
protein_sequences = {}

output_file = 'gossypii_protein_sequences.csv'
with open(output_file, 'w') as f:
    f.write('Protein Name,Sequence\n')

# Find whether each protein is in Eremothecium gossypii
for protein in protein_names:
    try:
        query = f'(organism_id:33169) AND ({protein})'
        params = {'query': query, 'format': 'json', 'limit': 1}

        response = requests.get(uniprot_url, params=params, timeout=60)
        print(response.text)
        if response.status_code == 200:
            results = response.json()
            if results['results']:
                sequence = results['results'][0]['sequence']['value']
                protein_sequences[protein] = sequence
                print(f'Found: {protein}')
            else:
                protein_sequences[protein] = 'Not Found'
                print(f'Not Found: {protein}')
        else:
            protein_sequences[protein] = 'Error'
            print(f'Error for {protein}')

        with open(output_file, 'a') as f:
            f.write(f'{protein},{protein_sequences[protein]}\n')

    except requests.exceptions.Timeout:
        print(f"Timeout for {protein}, skipping...")
        with open(output_file, 'a') as f:
            f.write(f'{protein},Timeout\n')

    except requests.exceptions.ConnectionError:
        print(f"Connection error for {protein}, retrying...")
        time.sleep(5)
        continue
print(f"Query completed, the result are saved in the '{output_file} file")
