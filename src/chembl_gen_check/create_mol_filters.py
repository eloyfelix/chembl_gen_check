from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from molbloom import CustomFilter
from chembl_gen_check.lacan import mol_to_pairs
from chembl_gen_check.ring_systems import RingSystemFinder
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import chembl_downloader
from tqdm import tqdm
import threading
import queue
import logging
import pickle
import math
import argparse
import csv
import os

RDLogger.DisableLog("rdApp.*")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Sentinel value to signal end of processing
_SENTINEL = None


def batched(iterable, n):
    """Yield successive n-sized chunks from iterable (memory efficient)."""
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == n:
            yield batch
            batch = []
    if batch:
        yield batch


def process_lacan_chunk(mol_chunk):
    """Process a chunk of molecules and return combined LACAN profile for the chunk."""
    idx_counter = Counter()
    pair_counter = Counter()
    total_pairs = 0
    
    for mol in mol_chunk:
        if mol is None:
            continue
        try:
            pairs = mol_to_pairs(mol)
            pair_counter.update(pairs)
            for a, b in pairs:
                idx_counter[a] += 1
                idx_counter[b] += 1
            total_pairs += len(pairs)
        except Exception:
            continue
    
    return idx_counter, pair_counter, total_pairs


def get_lacan_profile_streaming(mol_batches, n_workers=None, chunk_size=1000, pbar=None):
    """
    Stream-process molecule batches for LACAN profile with incremental merging.
    Processes entire chunks in workers instead of per-molecule futures.
    """
    idx_counter = Counter()
    pair_counter = Counter()
    total_setsize = 0
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for mol_batch in mol_batches:
            # Filter None molecules
            valid_mols = [m for m in mol_batch if m is not None]
            if not valid_mols:
                continue
            
            # Submit chunk-level work for better efficiency
            n_chunks = max(1, len(valid_mols) // chunk_size)
            chunk_size_actual = max(1, len(valid_mols) // n_chunks)
            
            futures = []
            for i in range(0, len(valid_mols), chunk_size_actual):
                chunk = valid_mols[i:i + chunk_size_actual]
                futures.append(executor.submit(process_lacan_chunk, chunk))
            
            # Incrementally merge results as they complete
            for future in as_completed(futures):
                try:
                    i_count, p_count, size = future.result()
                    idx_counter.update(i_count)
                    pair_counter.update(p_count)
                    total_setsize += size
                except Exception as e:
                    logging.error(f"Error in LACAN chunk: {e}")
            
            if pbar:
                pbar.update(len(mol_batch))
    
    # Keep only top 1023 idx occurrences
    idx_occurrences = dict(idx_counter.most_common(1023))
    return {"idx": idx_occurrences, "pairs": dict(pair_counter), "setsize": total_setsize}


def get_unique_scaffolds_streaming(mol_batches, n_workers=None, chunk_size=1000, pbar=None):
    """Stream-process molecule batches for scaffolds/skeletons."""
    unique_scaffolds = set()
    unique_skeletons = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for mol_batch in mol_batches:
            valid_mols = [m for m in mol_batch if m is not None]
            if not valid_mols:
                continue
            
            # Dynamic chunk sizing based on batch size
            n_chunks = max(1, len(valid_mols) // chunk_size)
            chunk_size_actual = max(1, len(valid_mols) // n_chunks)
            
            futures = []
            for i in range(0, len(valid_mols), chunk_size_actual):
                chunk = valid_mols[i:i + chunk_size_actual]
                futures.append(executor.submit(process_scaffold_chunk, chunk))

            for future in as_completed(futures):
                try:
                    scaffolds, skeletons = future.result()
                    unique_scaffolds.update(scaffolds)
                    unique_skeletons.update(skeletons)
                except Exception as e:
                    logging.error(f"Error processing scaffold chunk: {e}")
            
            if pbar:
                pbar.update(len(mol_batch))

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
        except Exception:
            continue
    return scaffolds, skeletons


def get_unique_ring_systems_streaming(mol_batches, n_workers=None, chunk_size=1000, pbar=None):
    """Stream-process molecule batches for ring systems."""
    unique_rings = set()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for mol_batch in mol_batches:
            valid_mols = [m for m in mol_batch if m is not None]
            if not valid_mols:
                continue
            
            # Dynamic chunk sizing
            n_chunks = max(1, len(valid_mols) // chunk_size)
            chunk_size_actual = max(1, len(valid_mols) // n_chunks)
            
            futures = []
            for i in range(0, len(valid_mols), chunk_size_actual):
                chunk = valid_mols[i:i + chunk_size_actual]
                futures.append(executor.submit(process_ring_systems_chunk, chunk))

            for future in as_completed(futures):
                try:
                    rings = future.result()
                    unique_rings.update(rings)
                except Exception as e:
                    logging.error(f"Error processing ring chunk: {e}")
            
            if pbar:
                pbar.update(len(mol_batch))

    return unique_rings


def process_ring_systems_chunk(mol_chunk):
    """Process a chunk of molecules for ring systems (creates RingSystemFinder per chunk)."""
    ring_finder = RingSystemFinder()
    rings = set()
    for mol in mol_chunk:
        if mol is None:
            continue
        try:
            found_rings = ring_finder.find_ring_systems(mol, as_mols=False)
            rings.update(found_rings)
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


def parse_smiles_batch(smiles_batch, n_workers=None):
    """Parse a batch of SMILES strings to RDKit Mol objects in parallel."""
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        mols = list(executor.map(smiles_to_mol, smiles_batch, chunksize=1000))
    return mols


def smiles_generator_from_chembl(version, chembl_version_int):
    """Generator that yields SMILES strings from ChEMBL database."""
    if chembl_version_int == 8:
        compounds_table = "compounds"
    else:
        compounds_table = "compound_structures"
    
    with chembl_downloader.connect(version=version) as conn:
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
        for row in cursor:
            yield row[0]


def smiles_generator_from_tsv(tsv_file):
    """Generator that yields SMILES strings from a TSV file."""
    with open(tsv_file, "r") as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if "smiles" in row and row["smiles"]:
                yield row["smiles"].strip()


def smiles_generator_from_parquet(parquet_file, batch_size=50000):
    """Generator that yields SMILES strings from a parquet file using pyarrow streaming."""
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(parquet_file)
    for batch in pf.iter_batches(batch_size=batch_size, columns=['smiles']):
        for smiles in batch.column('smiles').to_pylist():
            if smiles:
                yield smiles


def get_parquet_row_count(parquet_file):
    """Get row count from parquet file metadata (no full scan needed)."""
    import pyarrow.parquet as pq
    return pq.ParquetFile(parquet_file).metadata.num_rows


def download_surechembl(output_path="compounds.parquet"):
    """Download SureChEMBL compounds parquet file."""
    import urllib.request
    url = "https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/latest/compounds.parquet"
    
    if os.path.exists(output_path):
        logging.info(f"SureChEMBL file already exists: {output_path}")
        return output_path
    
    logging.info(f"Downloading SureChEMBL from {url}...")
    
    def report_progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(100, downloaded * 100 / total_size)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)
            print(f"\rDownloading: {mb_downloaded:.1f}/{mb_total:.1f} MB ({percent:.1f}%)", end="", flush=True)
    
    urllib.request.urlretrieve(url, output_path, reporthook=report_progress)
    print()  # newline after progress
    logging.info(f"Downloaded to {output_path}")
    return output_path


def count_lines(filepath):
    """Count lines in a file efficiently."""
    with open(filepath, 'rb') as f:
        line_count = sum(1 for _ in f)
    # Subtract 1 for header, but never return a negative count
    return max(line_count - 1, 0)


class ParallelPipelineProcessor:
    """
    Orchestrates parallel extraction pipelines with memory-efficient batch streaming.
    
    Architecture:
    - Parser thread: reads SMILES → parses to Mol objects → feeds mol_queue
    - Dispatcher thread: takes from mol_queue → distributes to 3 extraction queues
    - 3 extraction threads (scaffold, ring, LACAN): each with ProcessPoolExecutor
    
    This overlaps parsing of batch N+1 with extraction of batch N.
    """
    
    def __init__(self, batch_size=50000, n_workers=None, chunk_size=1000, parse_ahead=2):
        self.batch_size = batch_size
        self.n_workers = n_workers or os.cpu_count()
        self.chunk_size = chunk_size
        self.parse_ahead = parse_ahead  # Number of batches to parse ahead
        
        # Results storage
        self.unique_scaffolds = set()
        self.unique_skeletons = set()
        self.unique_ring_systems = set()
        self.lacan_idx_counter = Counter()
        self.lacan_pair_counter = Counter()
        self.lacan_setsize = 0
        
        # Thread-safe locks for result updates
        self._scaffold_lock = threading.Lock()
        self._ring_lock = threading.Lock()
        self._lacan_lock = threading.Lock()
        
        # Progress tracking
        self._processed_scaffold = 0
        self._processed_ring = 0
        self._processed_lacan = 0
    
    def _parse_producer(self, smiles_generator, mol_queue, pbar_parse):
        """
        Producer thread: reads SMILES batches and parses them to Mol objects.
        Runs ahead of extraction to overlap parsing with processing.
        """
        try:
            for smiles_batch in batched(smiles_generator, self.batch_size):
                mol_batch = parse_smiles_batch(smiles_batch, n_workers=self.n_workers)
                if pbar_parse:
                    pbar_parse.update(len(mol_batch))
                mol_queue.put(mol_batch)
        finally:
            mol_queue.put(_SENTINEL)
    
    def _dispatcher(self, mol_queue, extraction_queues):
        """
        Dispatcher thread: takes parsed mol batches and distributes to all extraction queues.
        This decouples parsing from extraction queue management.
        """
        while True:
            mol_batch = mol_queue.get()
            if mol_batch is _SENTINEL:
                # Send sentinel to all extraction queues
                for q in extraction_queues:
                    q.put(_SENTINEL)
                mol_queue.task_done()
                break
            
            # Distribute to all extraction queues
            for q in extraction_queues:
                q.put(mol_batch)
            mol_queue.task_done()
    
    def _process_scaffolds_worker(self, mol_queue, do_scaffold, pbar):
        """Worker thread for scaffold/skeleton extraction."""
        if not do_scaffold:
            # Drain the queue without processing
            while True:
                item = mol_queue.get()
                if item is _SENTINEL:
                    mol_queue.task_done()
                    break
                mol_queue.task_done()
            return
        
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            while True:
                mol_batch = mol_queue.get()
                if mol_batch is _SENTINEL:
                    mol_queue.task_done()
                    break
                
                valid_mols = [m for m in mol_batch if m is not None]
                if valid_mols:
                    # Submit chunk-level work
                    futures = []
                    for i in range(0, len(valid_mols), self.chunk_size):
                        chunk = valid_mols[i:i + self.chunk_size]
                        futures.append(executor.submit(process_scaffold_chunk, chunk))
                    
                    for future in as_completed(futures):
                        try:
                            scaffolds, skeletons = future.result()
                            with self._scaffold_lock:
                                self.unique_scaffolds.update(scaffolds)
                                self.unique_skeletons.update(skeletons)
                        except Exception as e:
                            logging.error(f"Scaffold error: {e}")
                
                with self._scaffold_lock:
                    self._processed_scaffold += len(mol_batch)
                    if pbar:
                        pbar.update(len(mol_batch))
                
                mol_queue.task_done()
    
    def _process_rings_worker(self, mol_queue, do_rings, pbar):
        """Worker thread for ring system extraction."""
        if not do_rings:
            while True:
                item = mol_queue.get()
                if item is _SENTINEL:
                    mol_queue.task_done()
                    break
                mol_queue.task_done()
            return
        
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            while True:
                mol_batch = mol_queue.get()
                if mol_batch is _SENTINEL:
                    mol_queue.task_done()
                    break
                
                valid_mols = [m for m in mol_batch if m is not None]
                if valid_mols:
                    futures = []
                    for i in range(0, len(valid_mols), self.chunk_size):
                        chunk = valid_mols[i:i + self.chunk_size]
                        futures.append(executor.submit(process_ring_systems_chunk, chunk))
                    
                    for future in as_completed(futures):
                        try:
                            rings = future.result()
                            with self._ring_lock:
                                self.unique_ring_systems.update(rings)
                        except Exception as e:
                            logging.error(f"Ring system error: {e}")
                
                with self._ring_lock:
                    self._processed_ring += len(mol_batch)
                    if pbar:
                        pbar.update(len(mol_batch))
                
                mol_queue.task_done()
    
    def _process_lacan_worker(self, mol_queue, do_lacan, pbar):
        """Worker thread for LACAN profile extraction."""
        if not do_lacan:
            while True:
                item = mol_queue.get()
                if item is _SENTINEL:
                    mol_queue.task_done()
                    break
                mol_queue.task_done()
            return
        
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            while True:
                mol_batch = mol_queue.get()
                if mol_batch is _SENTINEL:
                    mol_queue.task_done()
                    break
                
                valid_mols = [m for m in mol_batch if m is not None]
                if valid_mols:
                    futures = []
                    for i in range(0, len(valid_mols), self.chunk_size):
                        chunk = valid_mols[i:i + self.chunk_size]
                        futures.append(executor.submit(process_lacan_chunk, chunk))
                    
                    for future in as_completed(futures):
                        try:
                            i_count, p_count, size = future.result()
                            with self._lacan_lock:
                                self.lacan_idx_counter.update(i_count)
                                self.lacan_pair_counter.update(p_count)
                                self.lacan_setsize += size
                        except Exception as e:
                            logging.error(f"LACAN error: {e}")
                
                with self._lacan_lock:
                    self._processed_lacan += len(mol_batch)
                    if pbar:
                        pbar.update(len(mol_batch))
                
                mol_queue.task_done()
    
    def process(self, smiles_generator, total_count=None, do_scaffold=True, do_rings=True, do_lacan=True):
        """
        Process molecules through parallel extraction pipelines.
        
        Pipeline architecture:
        - Parser thread → mol_queue → Dispatcher thread → 3 extraction queues → 3 extraction threads
        
        This overlaps:
        - Parsing of batch N+1 while extraction of batch N is running
        - All three extractions (scaffold, ring, LACAN) run concurrently
        
        Args:
            smiles_generator: Iterator yielding SMILES strings
            total_count: Optional total count for progress bar
            do_scaffold: Whether to extract scaffolds/skeletons
            do_rings: Whether to extract ring systems
            do_lacan: Whether to build LACAN profile
        """
        # Queue for parsed molecules (between parser and dispatcher)
        mol_queue = queue.Queue(maxsize=self.parse_ahead)
        
        # Queues for each extraction pipeline
        scaffold_queue = queue.Queue(maxsize=3)
        ring_queue = queue.Queue(maxsize=3)
        lacan_queue = queue.Queue(maxsize=3)
        extraction_queues = [scaffold_queue, ring_queue, lacan_queue]
        
        # Create progress bars
        pbar_parse = tqdm(total=total_count, desc="Parsing SMILES", smoothing=0)
        pbar_scaffold = tqdm(total=total_count, desc="Scaffolds", smoothing=0) if do_scaffold else None
        pbar_ring = tqdm(total=total_count, desc="Ring systems", smoothing=0) if do_rings else None
        pbar_lacan = tqdm(total=total_count, desc="LACAN", smoothing=0) if do_lacan else None
        
        # Start all threads
        threads = [
            # Producer: parses SMILES to molecules
            threading.Thread(target=self._parse_producer, args=(smiles_generator, mol_queue, pbar_parse), name="parser"),
            # Dispatcher: distributes mol batches to extraction queues
            threading.Thread(target=self._dispatcher, args=(mol_queue, extraction_queues), name="dispatcher"),
            # Extraction workers
            threading.Thread(target=self._process_scaffolds_worker, args=(scaffold_queue, do_scaffold, pbar_scaffold), name="scaffold"),
            threading.Thread(target=self._process_rings_worker, args=(ring_queue, do_rings, pbar_ring), name="ring"),
            threading.Thread(target=self._process_lacan_worker, args=(lacan_queue, do_lacan, pbar_lacan), name="lacan"),
        ]
        for t in threads:
            t.start()
        
        # Wait for all threads to complete
        for t in threads:
            t.join()
        
        # Close progress bars
        pbar_parse.close()
        if pbar_scaffold:
            pbar_scaffold.close()
        if pbar_ring:
            pbar_ring.close()
        if pbar_lacan:
            pbar_lacan.close()
    
    def get_lacan_profile(self, size=1024):
        """Get the final LACAN profile dictionary."""
        idx_occurrences = dict(self.lacan_idx_counter.most_common(size - 1))
        return {
            "idx": idx_occurrences,
            "pairs": dict(self.lacan_pair_counter),
            "setsize": self.lacan_setsize
        }


def main():
    parser = argparse.ArgumentParser(
        description="Generate ChEMBL filters for molecular generation validation"
    )
    parser.add_argument(
        "--chembl_version",
        type=int,
        default=36,
        help="ChEMBL database version to use (integer, minimum 8, default: 36)",
        choices=range(8, 50),
    )
    parser.add_argument(
        "--tsv_file",
        type=str,
        default=None,
        help="Path to TSV file containing molecules as SMILES strings instead of ChEMBL",
    )
    parser.add_argument(
        "--surechembl",
        action="store_true",
        default=False,
        help="Use SureChEMBL dataset (downloads ~1GB parquet file if not present)",
    )
    parser.add_argument(
        "--parquet_file",
        type=str,
        default=None,
        help="Path to parquet file containing molecules with 'smiles' column",
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
    parser.add_argument(
        "--batch_size",
        type=int,
        default=50000,
        help="Batch size for streaming processing (default: 50000)",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=1000,
        help="Chunk size for parallel workers (default: 1000)",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: CPU count)",
    )
    parser.add_argument(
        "--parse_ahead",
        type=int,
        default=2,
        help="Number of batches to parse ahead of extraction (default: 2)",
    )
    args = parser.parse_args()

    formatted_version = (
        f"{args.chembl_version:02d}"
        if args.chembl_version < 10
        else str(args.chembl_version)
    )

    # Set up SMILES generator and count
    total_count = None
    
    if args.surechembl:
        # SureChEMBL mode - download parquet and process
        parquet_path = download_surechembl()
        total_count = get_parquet_row_count(parquet_path)
        logging.info(f"SureChEMBL dataset: {total_count:,} compounds")
        smiles_gen = smiles_generator_from_parquet(parquet_path, batch_size=args.batch_size)
    elif args.parquet_file:
        # Custom parquet file
        logging.info(f"Using parquet file: {args.parquet_file}")
        total_count = get_parquet_row_count(args.parquet_file)
        smiles_gen = smiles_generator_from_parquet(args.parquet_file, batch_size=args.batch_size)
    elif args.tsv_file:
        # TSV file
        logging.info(f"Using TSV file: {args.tsv_file}")
        total_count = count_lines(args.tsv_file)
        smiles_gen = smiles_generator_from_tsv(args.tsv_file)
    else:
        # Default: ChEMBL
        logging.info(f"Downloading/extracting data from ChEMBL version {formatted_version}")
        chembl_downloader.download_extract_sqlite(version=formatted_version)
        
        # Get count for progress bar
        with chembl_downloader.connect(version=formatted_version) as conn:
            cursor = conn.cursor()
            compounds_table = "compounds" if args.chembl_version == 8 else "compound_structures"
            cursor.execute(f"""
                SELECT COUNT(*) FROM {compounds_table}
                WHERE molregno IN (SELECT DISTINCT parent_molregno FROM molecule_hierarchy)
                AND canonical_smiles IS NOT NULL
            """)
            total_count = cursor.fetchone()[0]
        
        smiles_gen = smiles_generator_from_chembl(formatted_version, args.chembl_version)

    if total_count == 0:
        logging.error("No SMILES data found. Exiting.")
        exit(1)

    logging.info(f"Processing {total_count} molecules with batch_size={args.batch_size}, chunk_size={args.chunk_size}, parse_ahead={args.parse_ahead}")

    # Process through parallel pipelines
    processor = ParallelPipelineProcessor(
        batch_size=args.batch_size,
        n_workers=args.n_workers,
        chunk_size=args.chunk_size,
        parse_ahead=args.parse_ahead
    )
    
    processor.process(
        smiles_gen,
        total_count=total_count,
        do_scaffold=args.scaffold,
        do_rings=args.ring_system,
        do_lacan=args.lacan
    )

    # Determine output prefix based on data source
    if args.surechembl:
        output_prefix = "surechembl_"
    elif args.parquet_file:
        parquet_stem = os.path.splitext(os.path.basename(args.parquet_file))[0]
        output_prefix = f"{parquet_stem}_"
    elif args.tsv_file:
        tsv_stem = os.path.splitext(os.path.basename(args.tsv_file))[0]
        output_prefix = f"{tsv_stem}_"
    else:
        output_prefix = f"chembl_{formatted_version}_"

    # Create and save bloom filters
    if args.scaffold:
        create_bloom_filter(processor.unique_scaffolds, f"{output_prefix}scaffold")
        create_bloom_filter(processor.unique_skeletons, f"{output_prefix}skeleton")
        logging.info(f"Found {len(processor.unique_scaffolds)} unique scaffolds.")
        logging.info(f"Found {len(processor.unique_skeletons)} unique skeletons.")

    if args.ring_system:
        create_bloom_filter(processor.unique_ring_systems, f"{output_prefix}ring_system")
        logging.info(f"Found {len(processor.unique_ring_systems)} unique ring systems.")

    if args.lacan:
        lacan_profile = processor.get_lacan_profile()
        lacan_filename = f"{output_prefix}lacan.pkl"
        with open(lacan_filename, "wb") as file:
            pickle.dump(lacan_profile, file)
        logging.info(f"LACAN profile saved to {lacan_filename}")


if __name__ == "__main__":
    main()
