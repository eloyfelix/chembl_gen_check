import warnings
# Filter out RDKit converter warnings
warnings.filterwarnings('ignore', category=RuntimeWarning,
                       message='to-Python converter for boost::shared_ptr.*')

from importlib.resources import files
from rdkit.Chem.Scaffolds import MurckoScaffold
from .ring_systems import RingSystemFinder
from .lacan import score_mol
from molbloom import BloomFilter
from rdkit.Chem import FilterCatalog
from rdkit import Chem
import pickle

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(params.FilterCatalogs.CHEMBL)
sa_catalog = FilterCatalog.FilterCatalog(params)


databases = {
    "chembl": {
        "scaffold": "chembl_scaffold.bloom",
        "ring_system": "chembl_ring_system.bloom",
        "lacan_profile": "chembl_lacan.pkl",
    },
}


class GenVal:
    scaffold_filter = None
    ring_sytem_filter = None
    lacal_profile = None

    def __init__(self, db_name: str = "chembl") -> None:
        s_file_path = files("genval.data").joinpath(databases[db_name]["scaffold"])
        self.scaffold_filter = BloomFilter(str(s_file_path))
        rs_file_path = files("genval.data").joinpath(databases[db_name]["ring_system"])
        self.ring_sytem_filter = BloomFilter(str(rs_file_path))
        lp_file_path = files("genval.data").joinpath(
            databases[db_name]["lacan_profile"]
        )
        with open(lp_file_path, "rb") as file:
            self.lacan_profile = pickle.load(file)

    def check_scaffold(self, smiles: str) -> bool:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            return scaffold_smiles in self.scaffold_filter
        except:
            return False

    def check_ring_systems(self, smiles: str) -> bool:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        try:
            ring_system_finder = RingSystemFinder()
            ring_systems = ring_system_finder.find_ring_systems(mol, as_mols=False)
            for rs in ring_systems:
                if rs not in self.ring_sytem_filter:
                    return False
            return True
        except:
            return False

    def check_structural_alerts(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        return len(sa_catalog.GetMatches(mol))

    def check_lacan(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        return score_mol(mol, self.lacan_profile)
