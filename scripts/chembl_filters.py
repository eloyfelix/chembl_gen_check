from rdkit.Chem.Scaffolds import MurckoScaffold
from molbloom import CustomFilter
from rdkit import Chem, RDLogger
from io import BytesIO
import requests
import gzip
import math
import re
import logging
from genval.ring_systems import RingSystemFinder
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


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


def download_and_extract_smiles(file_url):
    response = requests.get(file_url, stream=True)
    response.raise_for_status()

    smiles_list = []
    with gzip.open(BytesIO(response.content), "rt") as f:
        # Skip the header line
        next(f)
        for line in f:
            try:
                smiles = line.split("\t")[1].strip()
                if smiles:
                    smiles_list.append(smiles)
            except IndexError:
                continue
    return smiles_list


def get_unique_scaffolds(mol_list):
    unique_scaffolds = set()
    for mol in tqdm(mol_list, desc="Processing scaffolds"):
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            unique_scaffolds.add(scaffold_smiles)
        except:
            continue
    return unique_scaffolds


def get_unique_ring_systems(mol_list):
    ring_finder = RingSystemFinder()
    unique_rings = set()
    for mol in tqdm(mol_list, desc="Processing ring systems"):
        try:
            rings = ring_finder.find_ring_systems(mol, as_mols=False)
            for r in rings:
                unique_rings.add(r)
        except Exception:
            continue
    return unique_rings


def calc_m(epsilon, N):
    M = -(N * math.log(epsilon)) / (math.log(2) ** 2)
    return math.ceil(M)


def create_bloom_filter(items, filter_id):
    m = calc_m(0.000025, len(items))
    n = len(items)
    logging.info(f"Creating bloom filter {filter_id} with {m} bits for {n} items")
    bf = CustomFilter(m, n, filter_id)
    for item in items:
        bf.add(item)
    bf.save(f"{filter_id}.bloom")
    logging.info(f"Bloom filter saved to {filter_id}.bloom")


if __name__ == "__main__":
    logging.info("Starting ChEMBL filter creation")
    file_url = get_latest_chembl_file_url()
    if file_url:
        logging.info(f"Downloading data from {file_url}")
        smiles_list = download_and_extract_smiles(file_url)

        mol_list = []
        for smiles in tqdm(smiles_list, desc="Parsing SMILES"):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_list.append(mol)

        unique_scaffolds = get_unique_scaffolds(mol_list)
        create_bloom_filter(unique_scaffolds, "chembl_scaffold")

        unique_ring_systems = get_unique_ring_systems(mol_list)
        create_bloom_filter(unique_ring_systems, "chembl_ring_system")
    else:
        logging.error("Could not find a valid file URL.")
