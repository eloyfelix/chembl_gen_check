from importlib.resources import files
from molbloom import BloomFilter
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Chem
from typing import Optional


databases = {
    "chembl_scaffold": "chembl_scaffold.bloom",
    "surechembl_scaffold": "surechembl_scaffold.bloom",
}


class GenVal:

    bf = None

    def __init__(self, filter_name: str = "chembl_scaffold") -> None:
        file_path = files("genval.data").joinpath(databases[filter_name])
        self.bf = BloomFilter(file_path)

    def validate_smiles(self, smiles: str) -> Optional[bool]:
        valid = None
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                valid = scaffold_smiles in self.bf
            except:
                pass
        return valid
