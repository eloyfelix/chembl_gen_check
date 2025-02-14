from rdkit.Chem.Scaffolds import MurckoScaffold
from molbloom import CustomFilter
from rdkit import Chem, RDLogger
from io import BytesIO
import pandas as pd
import requests
import gzip
import math
import re

RDLogger.DisableLog('rdApp.warning')

def get_latest_chembl_file_url():
    base_url = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    response = requests.get(base_url)

    pattern = re.compile(r"chembl_(\d+)_chemreps\.txt\.gz")    
    if response.status_code == 200:
        for line in response.text.splitlines():
            match = pattern.search(line)
            if match:
                file_url = base_url + match.group(0)
                print(f"Found latest file: {file_url}")
                return file_url
    return None

def download_and_extract_smiles(file_url):
    response = requests.get(file_url, stream=True)
    response.raise_for_status()
    
    with gzip.open(BytesIO(response.content), 'rt') as f:
        data = pd.read_csv(f, sep="\t", usecols=[1], names=["canonical_smiles"], header=0)
    return data["canonical_smiles"].dropna()

def get_unique_scaffolds(smiles_list):
    unique_scaffolds = set()
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                unique_scaffolds.add(scaffold_smiles)
            except:
                continue
    return unique_scaffolds


def calc_m(epsilon, N):
    M = - (N * math.log(epsilon)) / (math.log(2) ** 2)
    return math.ceil(M)


def create_bloom_filter(scaffolds):
    m = calc_m(0.000025, len(scaffolds))
    n = len(scaffolds)
    print(m, n)
    # number of bits, number of molecules
    bf = CustomFilter(m, n, 'chembl_scaffolds')
    for scaffold in scaffolds:
        bf.add(scaffold)
    
    bf.save('chembl_scaffolds.bloom')
    print("Bloom filter saved to chembl_scaffolds.bloom")

file_url = get_latest_chembl_file_url()
if file_url:
    smiles_list = download_and_extract_smiles(file_url)
    unique_scaffolds = get_unique_scaffolds(smiles_list)
    create_bloom_filter(unique_scaffolds)
