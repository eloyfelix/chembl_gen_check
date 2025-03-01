from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from molbloom import CustomFilter
from chembl_gen_check.lacan import mol_to_pairs
from chembl_gen_check.ring_systems import RingSystemFinder
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import chembl_downloader
from tqdm import tqdm
import logging
import pickle
import math

RDLogger.DisableLog("rdApp.*")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def get_lacan_profile_for_mol(mol):
    idx_counter = Counter()
    pair_counter = Counter()
    pairs = mol_to_pairs(mol)
    pair_counter.update(pairs)
    for a, b in pairs:
        idx_counter[a] += 1
        idx_counter[b] += 1

    return idx_counter, pair_counter, len(pairs)


def combine_lacan_profiles(profiles, size=1024):
    idx_counter = Counter()
    pair_counter = Counter()
    setsize = 0

    for i_count, p_count, s in profiles:
        idx_counter.update(i_count)
        pair_counter.update(p_count)
        setsize += s

    idx_occurrences = dict(idx_counter.most_common(size - 1))
    return {"idx": idx_occurrences, "pairs": dict(pair_counter), "setsize": setsize}


def get_lacan_profile_for_mols(mol_list, size=1024, n_workers=None, chunk_size=1000):
    profiles = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for i in tqdm(
            range(0, len(mol_list), chunk_size),
            desc="Processing molecules for LACAN profile",
        ):
            chunk = mol_list[i : i + chunk_size]
            futures = [executor.submit(get_lacan_profile_for_mol, mol) for mol in chunk]
            for future in futures:
                profiles.append(future.result())
    return combine_lacan_profiles(profiles, size)


def get_unique_scaffolds(mol_list, n_workers=None, chunk_size=1000):
    unique_scaffolds = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for i in range(0, len(mol_list), chunk_size):
            chunk = mol_list[i : i + chunk_size]
            futures.append(executor.submit(process_scaffold_chunk, chunk))

        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Processing scaffolds"
        ):
            try:
                scaffolds = future.result()
                unique_scaffolds.update(scaffolds)
            except Exception as e:
                logging.error(f"Error processing a chunk: {e}")
                continue

    return unique_scaffolds


def process_scaffold_chunk(mol_list):
    scaffolds = set()
    for mol in mol_list:
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            scaffolds.add(scaffold_smiles)
        except:
            continue
    return scaffolds


def get_unique_ring_systems(mol_list, n_workers=None, chunk_size=1000):
    ring_finder = RingSystemFinder()
    unique_rings = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for i in range(0, len(mol_list), chunk_size):
            chunk = mol_list[i : i + chunk_size]
            futures.append(executor.submit(process_ring_systems, chunk, ring_finder))

        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Processing ring systems"
        ):
            try:
                rings = future.result()
                unique_rings.update(rings)
            except Exception as e:
                logging.error(f"Error processing a chunk: {e}")
                continue

    return unique_rings


def process_ring_systems(mol_list, ring_finder):
    rings = set()
    for mol in mol_list:
        try:
            found_rings = ring_finder.find_ring_systems(mol, as_mols=False)
            for r in found_rings:
                rings.add(r)
        except Exception:
            continue
    return rings


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


def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol


if __name__ == "__main__":
    logging.info("Downloading/extracting data from ChEMBL")
    path = chembl_downloader.download_extract_sqlite()

    smiles_list = []
    with chembl_downloader.connect() as conn:
        cursor = conn.cursor()
        query = """
        SELECT canonical_smiles
        FROM compound_structures
        WHERE molregno IN (
            SELECT DISTINCT parent_molregno
            FROM molecule_hierarchy
        ) AND canonical_smiles IS NOT NULL
        """
        cursor.execute(query)
        for row in cursor.fetchall():
            smiles_list.append(row[0])

    if not smiles_list:
        logging.error("No SMILES data extracted. Exiting.")
        exit(1)

    with ProcessPoolExecutor() as executor:
        mol_list = list(
            tqdm(
                executor.map(smiles_to_mol, smiles_list, chunksize=10000),
                total=len(smiles_list),
                desc="Parsing SMILES",
                smoothing=0,
            )
        )
    mol_list = [mol for mol in mol_list if mol is not None]

    unique_scaffolds = get_unique_scaffolds(mol_list)
    create_bloom_filter(unique_scaffolds, "chembl_scaffold")

    unique_ring_systems = get_unique_ring_systems(mol_list)
    create_bloom_filter(unique_ring_systems, "chembl_ring_system")
    lacan_profile = get_lacan_profile_for_mols(mol_list)
    with open("chembl_lacan.pkl", "wb") as file:
        pickle.dump(lacan_profile, file)
