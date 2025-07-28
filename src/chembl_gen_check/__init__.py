import warnings

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message="to-Python converter for boost::shared_ptr.*",
)

from importlib.resources import files
from pathlib import Path
from rdkit.Chem.Scaffolds import MurckoScaffold
from .ring_systems import RingSystemFinder
from .lacan import score_mol
from molbloom import BloomFilter
from rdkit.Chem import FilterCatalog
from rdkit import Chem, RDLogger
import pickle


RDLogger.DisableLog("rdApp.*")

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(params.FilterCatalogs.CHEMBL)
sa_catalog = FilterCatalog.FilterCatalog(params)


databases = {
    "chembl": {
        "scaffold": "chembl_scaffold.bloom",
        "skeleton": "chembl_skeleton.bloom",
        "ring_system": "chembl_ring_system.bloom",
        "lacan_profile": "chembl_lacan.pkl",
    },
}


class Checker:
    scaffold_filter = None
    skeleton_filter = None
    ring_sytem_filter = None
    lacan_profile = None

    def __init__(self, db_name: str = "chembl") -> None:
        folder_path = Path(db_name)
        
        # Check if db_name is a folder path
        if folder_path.is_dir():
            self._load_from_folder(folder_path)
        else:
            self._load_from_package(db_name)
    
    def _load_from_folder(self, folder_path: Path) -> None:
        required_files = ["scaffold.bloom", "skeleton.bloom", "ring_system.bloom", "lacan.pkl"]
        missing_files = [f for f in required_files if not (folder_path / f).exists()]
        
        if missing_files:
            raise FileNotFoundError(f"Missing required files in {folder_path}: {', '.join(missing_files)}")
        
        self.scaffold_filter = BloomFilter(str(folder_path / "scaffold.bloom"))
        self.skeleton_filter = BloomFilter(str(folder_path / "skeleton.bloom"))
        self.ring_sytem_filter = BloomFilter(str(folder_path / "ring_system.bloom"))
        
        with open(folder_path / "lacan.pkl", "rb") as file:
            self.lacan_profile = pickle.load(file)
    
    def _load_from_package(self, db_name: str) -> None:
        data_files = files("chembl_gen_check.data")
        db_config = databases[db_name]
        
        self.scaffold_filter = BloomFilter(str(data_files.joinpath(db_config["scaffold"])))
        self.skeleton_filter = BloomFilter(str(data_files.joinpath(db_config["skeleton"])))
        self.ring_sytem_filter = BloomFilter(str(data_files.joinpath(db_config["ring_system"])))
        
        with open(data_files.joinpath(db_config["lacan_profile"]), "rb") as file:
            self.lacan_profile = pickle.load(file)

    def load_smiles(self, smiles) -> None:
        self.mol = Chem.MolFromSmiles(smiles)

    def check_scaffold(self) -> bool:
        if not self.mol:
            return False
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            return scaffold_smiles in self.scaffold_filter
        except:
            return False

    def check_skeleton(self) -> bool:
        if not self.mol:
            return False
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)
            skeleton = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            skeleton_smiles = Chem.MolToSmiles(skeleton)
            return skeleton_smiles in self.skeleton_filter
        except:
            return False

    def check_ring_systems(self) -> bool:
        if not self.mol:
            return False
        try:
            ring_system_finder = RingSystemFinder()
            ring_systems = ring_system_finder.find_ring_systems(self.mol, as_mols=False)
            for rs in ring_systems:
                if rs not in self.ring_sytem_filter:
                    return False
            return True
        except:
            return False

    def check_structural_alerts(self) -> int:
        if not self.mol:
            return 9999  # return a large number to indicate invalid input
        return len(sa_catalog.GetMatches(self.mol))

    def check_lacan(self, t: float = 0.05, include_info: bool = False):
        if not self.mol:
            return (0.0, {"bad_bonds": []}) if include_info else 0.0
        result = score_mol(self.mol, self.lacan_profile, t)
        return result if include_info else result[0]
