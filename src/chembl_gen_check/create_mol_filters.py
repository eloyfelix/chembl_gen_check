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
import argparse
import csv

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
        pbar = tqdm(total=len(mol_list), desc="Processing molecules for LACAN profile")
        for i in range(0, len(mol_list), chunk_size):
            chunk = mol_list[i : i + chunk_size]
            futures = [executor.submit(get_lacan_profile_for_mol, mol) for mol in chunk]
            for future in futures:
                profiles.append(future.result())
                pbar.update(1)
        pbar.close()
    return combine_lacan_profiles(profiles, size)


def get_unique_scaffolds(mol_list, n_workers=None, chunk_size=1000):
    unique_scaffolds = set()
    unique_skeletons = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for i in range(0, len(mol_list), chunk_size):
            chunk = mol_list[i : i + chunk_size]
            futures.append(executor.submit(process_scaffold_chunk, chunk))

        pbar = tqdm(total=len(mol_list), desc="Processing scaffolds")
        for future in as_completed(futures):
            try:
                scaffolds, skeletons = future.result()
                unique_scaffolds.update(scaffolds)
                unique_skeletons.update(skeletons)
                pbar.update(chunk_size)
            except Exception as e:
                logging.error(f"Error processing a chunk: {e}")
                pbar.update(chunk_size)
                continue
        pbar.close()

    return unique_scaffolds, unique_skeletons


def process_scaffold_chunk(mol_list):
    scaffolds = set()
    skeletons = set()
    for mol in mol_list:
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            skeleton = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            skeleton_smiles = Chem.MolToSmiles(skeleton)
            scaffolds.add(scaffold_smiles)
            skeletons.add(skeleton_smiles)
        except:
            continue
    return scaffolds, skeletons


def get_unique_ring_systems(mol_list, n_workers=None, chunk_size=1000):
    ring_finder = RingSystemFinder()
    unique_rings = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for i in range(0, len(mol_list), chunk_size):
            chunk = mol_list[i : i + chunk_size]
            futures.append(executor.submit(process_ring_systems, chunk, ring_finder))

        pbar = tqdm(total=len(mol_list), desc="Processing ring systems")
        for future in as_completed(futures):
            try:
                rings = future.result()
                unique_rings.update(rings)
                pbar.update(chunk_size)
            except Exception as e:
                logging.error(f"Error processing a chunk: {e}")
                pbar.update(chunk_size)
                continue
        pbar.close()

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


def main():
    parser = argparse.ArgumentParser(
        description="Generate ChEMBL filters for molecular generation validation"
    )
    parser.add_argument(
        "--chembl_version",
        type=int,
        default=35,
        help="ChEMBL database version to use (integer, minimum 8, default: 35)",
        choices=range(8, 100),  # Setting upper limit to 100 for future versions
    )
    parser.add_argument(
        "--tsv_file",
        type=str,
        default=None,
        help="Path to TSV file containing molecules as SMILES strings instead of ChEMBL",
    )
    parser.add_argument(
        "--scaffold",
        action="store_true",
        default=True,
        help="Generate scaffold bloom filter (default: True)",
    )
    parser.add_argument(
        "--ring_system",
        action="store_true",
        default=True,
        help="Generate ring system bloom filter (default: True)",
    )
    parser.add_argument(
        "--lacan",
        action="store_true",
        default=True,
        help="Generate LACAN profile (default: True)",
    )
    args = parser.parse_args()

    formatted_version = (
        f"{args.chembl_version:02d}"
        if args.chembl_version < 10
        else str(args.chembl_version)
    )

    if not args.tsv_file:
        logging.info(
            f"Downloading/extracting data from ChEMBL version {formatted_version}"
        )
        path = chembl_downloader.download_extract_sqlite(version=formatted_version)

        if args.chembl_version == 8:
            compounds_table = "compounds"
        else:
            compounds_table = "compound_structures"

        smiles_list = []
        with chembl_downloader.connect(version=formatted_version) as conn:
            cursor = conn.cursor()

            query = f"""
            SELECT canonical_smiles
            FROM {compounds_table}
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
    else:
        logging.info(f"Using TSV file: {args.tsv_file}")
        smiles_list = []
        with open(args.tsv_file, "r") as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                if "smiles" in row and row["smiles"]:
                    smiles_list.append(row["smiles"].strip())

    with ProcessPoolExecutor() as executor:
        mol_list = list(
            tqdm(
                executor.map(smiles_to_mol, smiles_list, chunksize=10000),
                total=len(smiles_list),
                desc="Parsing SMILES",
                smoothing=0,
            )
        )

    if args.scaffold:
        unique_scaffolds, unique_skeletons = get_unique_scaffolds(mol_list)
        create_bloom_filter(unique_scaffolds, "scaffold")
        create_bloom_filter(unique_skeletons, "skeleton")
        logging.info(f"Found {len(unique_scaffolds)} unique scaffolds.")
        logging.info(f"Found {len(unique_skeletons)} unique skeletons.")

    if args.ring_system:
        unique_ring_systems = get_unique_ring_systems(mol_list)
        create_bloom_filter(unique_ring_systems, "ring_system")
        logging.info(f"Found {len(unique_ring_systems)} unique ring systems.")
    if args.lacan:
        lacan_profile = get_lacan_profile_for_mols(mol_list)
        with open("lacan.pkl", "wb") as file:
            pickle.dump(lacan_profile, file)


if __name__ == "__main__":
    main()
