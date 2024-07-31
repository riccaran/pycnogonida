print("Script configuration")
import os
import requests
import json
from tqdm import tqdm
import pandas as pd
import xml.etree.ElementTree as ET

# Ottieni AphiaID per Pycnogonida
url = "https://www.marinespecies.org/rest/AphiaIDByName/Pycnogonida"
response = requests.get(url)
aphia_id = response.json()

# Funzione per ottenere tutti i figli tassonomici ricorsivamente
def get_all_taxa(aphia_id):
    taxa = []
    offset = 1
    limit = 50
    while True:
        url = f"https://www.marinespecies.org/rest/AphiaChildrenByAphiaID/{aphia_id}?marine_only=true&offset={offset}&limit={limit}"
        response = requests.get(url)
        
        # Controlla se lo stato è 204, esci dal ciclo
        if response.status_code == 204:
            break
        
        try:
            data = response.json()
        except json.JSONDecodeError:
            break
        
        if not data:
            break
        
        taxa.extend(data)
        offset += limit

    # Ricorsione per ottenere tutti i figli di ciascun taxon        
    for taxon in taxa:
        if taxon.get('AphiaID'):
            print("\rRetrieving {}                                           ".format(taxon['scientificname']), end='')
            child_taxa = get_all_taxa(taxon['AphiaID'])
            taxon['children'] = child_taxa
    return taxa

# Funzione per ottenere tutte le occorrenze di una specie dato il suo nome scientifico con paginazione
def get_all_occurrences(scientific_name):
    occurrences = []
    offset = 0
    limit = 300  # Numero di occorrenze per pagina

    while True:
        url = f"https://api.gbif.org/v1/occurrence/search?scientificName={scientific_name}&limit={limit}&offset={offset}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'results' in data and len(data['results']) > 0:
                occurrences.extend(data['results'])
                offset += limit
            else:
                break
        else:
            break
    return occurrences

# Funzione per appiattire la struttura nidificata
def flatten_taxa(taxa, parent=None):
    flat_list = []
    for taxon in tqdm(taxa, desc='Taxa flattening for dataframe making'):
        flat_taxon = {
            'AphiaID': taxon.get('AphiaID'),
            'url': taxon.get('url'),
            'scientificname': taxon.get('scientificname'),
            'authority': taxon.get('authority'),
            'status': taxon.get('status'),
            'unacceptreason': taxon.get('unacceptreason'),
            'taxonRankID': taxon.get('taxonRankID'),
            'rank': taxon.get('rank'),
            'valid_AphiaID': taxon.get('valid_AphiaID'),
            'valid_name': taxon.get('valid_name'),
            'valid_authority': taxon.get('valid_authority'),
            'parentNameUsageID': taxon.get('parentNameUsageID'),
            'kingdom': taxon.get('kingdom'),
            'phylum': taxon.get('phylum'),
            'class': taxon.get('class'),
            'order': taxon.get('order'),
            'family': taxon.get('family'),
            'genus': taxon.get('genus'),
            'citation': taxon.get('citation'),
            'lsid': taxon.get('lsid'),
            'isMarine': taxon.get('isMarine'),
            'isBrackish': taxon.get('isBrackish'),
            'isFreshwater': taxon.get('isFreshwater'),
            'isTerrestrial': taxon.get('isTerrestrial'),
            'isExtinct': taxon.get('isExtinct'),
            'match_type': taxon.get('match_type'),
            'modified': taxon.get('modified'),
            'parent': parent
        }
        flat_list.append(flat_taxon)
        if 'children' in taxon and taxon['children']:
            flat_list.extend(flatten_taxa(taxon['children'], parent=taxon.get('AphiaID')))
    return flat_list

def retrieve_geographical_data(df):
    # Lista per memorizzare tutte le occorrenze
    all_occurrences = dict()

    # Itera attraverso il dataset e recupera tutte le occorrenze per ogni specie
    for scientific_name in tqdm(df['scientificname'], desc='Adding geographical data'):
        all_occurrences[scientific_name] = get_all_occurrences(scientific_name)

    return all_occurrences

def ensure_occurrencies_folder_exists(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Cartella '{folder_name}' creata.")
    else:
        print(f"Cartella '{folder_name}' già esistente.")

def fetch_genbank_data(species_name):
    # Struttura dati per salvare le informazioni
    species_data = []

    # Step 1: Search for the species in the NCBI nucleotide database
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        'db': 'nucleotide',
        'term': species_name,
        'retmode': 'xml',
        'retmax': '10'  # Adjust the number of results as needed
    }
    search_response = requests.get(search_url, params=search_params)
    search_tree = ET.fromstring(search_response.content)
    id_list = [id_elem.text for id_elem in search_tree.findall(".//Id")]

    for genbank_id in id_list:
        genbank_record = {'GenBank ID': genbank_id}
        
        # Step 2: Fetch summary for each GenBank ID
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            'db': 'nucleotide',
            'id': genbank_id,
            'retmode': 'xml'
        }
        summary_response = requests.get(summary_url, params=summary_params)
        summary_tree = ET.fromstring(summary_response.content)
        docsum = summary_tree.find(".//DocSum")
        if docsum is not None:
            for item in docsum.findall(".//Item"):
                name = item.get('Name')
                content = item.text
                genbank_record[name] = content

        # Step 3: Fetch detailed information (full record) for each GenBank ID
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            'db': 'nucleotide',
            'id': genbank_id,
            'retmode': 'xml',
            'rettype': 'gb'
        }
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_tree = ET.fromstring(fetch_response.content)
        
        for seq in fetch_tree.findall(".//GBSeq"):
            for child in seq:
                if child.tag == "GBSeq_feature-table":
                    features = []
                    for feature in child:
                        feature_data = {
                            'Feature Key': feature.find('GBFeature_key').text
                        }
                        for qual in feature.findall("GBFeature_quals/GBQualifier"):
                            qual_name = qual.find("GBQualifier_name").text
                            qual_value = qual.find("GBQualifier_value")
                            feature_data[qual_name] = qual_value.text if qual_value is not None else None
                        features.append(feature_data)
                    genbank_record['Features'] = features
                else:
                    genbank_record[child.tag] = child.text

        species_data.append(genbank_record)

    return species_data

def get_bold_data(species_name):
    """
    Ottiene dati dal BOLD Systems API per una data specie.
    
    :param species_name: Nome scientifico della specie
    :return: Risultati della query in formato JSON o None se la risposta è vuota o non valida.
    """
    url = "http://www.boldsystems.org/index.php/API_Public/combined"
    params = {
        'taxon': species_name,
        'format': 'json'
    }
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        try:
            return response.json()  # Prova a deserializzare il JSON
        except ValueError as e:
            return None
    else:
        return None

def fetch_obis_data(scientific_name):
    base_url = "https://api.obis.org/v3"
    endpoint = "/occurrence"
    
    params = {
        "scientificname": scientific_name
    }
    
    response = requests.get(base_url + endpoint, params=params)
    
    if response.status_code == 200:
        data = response.json()
        records = data['results']
        
        # Convert the results to a pandas DataFrame
        df = pd.DataFrame(records)
        return df
    else:
        return None

# Scarica tutti i taxa
pycnogonida_taxa = get_all_taxa(aphia_id)

# Appiattisci i dati
flat_data = flatten_taxa(pycnogonida_taxa)

print("Dataframe making")
# Crea un DataFrame da dati appiattiti
df = pd.DataFrame(flat_data)

print("Dataset processing")
pycnogonida = df[(df['rank'] == 'Species') & (df['status'] == 'accepted') & (df['isExtinct'].isna() | (df['isExtinct'] == 0))].reset_index(drop=True)

all_occurrences = retrieve_geographical_data(pycnogonida)

# GeneBank
ensure_occurrencies_folder_exists('genbank')
for species_name in tqdm(pycnogonida["scientificname"], desc='GenBank data retrieving'):
    species_data = fetch_genbank_data(species_name)
    species_df = pd.json_normalize(species_data)
    if not species_df.empty:
        species_df.to_csv('genbank/{}.tsv'.format(species_name.lower().replace(' ', '_')), sep='\t', index=False)

# BOLD
ensure_occurrencies_folder_exists('bold')
# Assumi che pycnogonida["scientificname"] contenga i nomi delle specie
for species_name in tqdm(pycnogonida["scientificname"], desc='BOLD data retrieving'):
    species_data = get_bold_data(species_name)
    if species_data:  # Controlla che i dati non siano None
        species_df = pd.json_normalize(species_data)
        if not species_df.empty:
            filename = f'bold/{species_name.lower().replace(" ", "_")}.tsv'
            species_df.to_csv(filename, sep='\t', index=False)
    else:
        pass

# OBIS
ensure_occurrencies_folder_exists('obis')
for species_name in tqdm(pycnogonida["scientificname"], desc='OBIS data retrieving'):
    species_df = fetch_obis_data(species_name)
    if not species_df.empty:
        species_df.to_csv('obis/{}.tsv'.format(species_name.lower().replace(' ', '_')), sep='\t', index=False)

print("Dataset saving")
pycnogonida.to_csv('pycnogonida.csv', index=False)

with open('all_occurrences.json', 'w') as json_file:
    json.dump(all_occurrences, json_file)

ensure_occurrencies_folder_exists('occurrencies')
for species, occurrences in all_occurrences.items():
    species_df = pd.DataFrame(occurrences)
    if species_df.shape != (0, 0):
        species_df.to_csv("gbif/{}.tsv".format(species.lower().replace(' ', '_')), sep = '\t', index = False)