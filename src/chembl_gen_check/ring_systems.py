from rdkit import Chem

# Based on Pat Walters' work:
# Blogpost: https://practicalcheminformatics.blogspot.com/2022/12/identifying-ring-systems-in-molecules.html
# Source: https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py

class RingSystemFinder:
    """A class to identify ring systems """

    def __init__(self):
        """Initialize susbstructure search objects to identify key functionality
        """
        self.ring_db_pat = Chem.MolFromSmarts("[#6R,#16R]=[OR0,SR0,CR0,NR0]")
        self.ring_atom_pat = Chem.MolFromSmarts("[R]")

    def tag_bonds_to_preserve(self, mol):
        """Assign the property "protected" to all ring carbonyls, etc.
        :param mol: input molecule
        :return: None
        """
        for bnd in mol.GetBonds():
            bnd.SetBoolProp("protected", False)
        for match in mol.GetSubstructMatches(self.ring_db_pat):
            bgn, end = match
            bnd = mol.GetBondBetweenAtoms(bgn, end)
            bnd.SetBoolProp("protected", True)

    @staticmethod
    def cleave_linker_bonds(mol):
        """Cleave bonds that are not in rings and not protected
        :param mol: input molecule
        :return: None
        """
        frag_bond_list = []
        for bnd in mol.GetBonds():
            if not bnd.IsInRing() and not bnd.GetBoolProp("protected") and bnd.GetBondType() == Chem.BondType.SINGLE:
                frag_bond_list.append(bnd.GetIdx())

        if len(frag_bond_list):
            frag_mol = Chem.FragmentOnBonds(mol, frag_bond_list)
            Chem.SanitizeMol(frag_mol)
            return frag_mol
        else:
            return mol

    def cleanup_fragments(self, mol, keep_dummy=False):
        """Split a molecule containing multiple ring systems into individual ring systems
        :param keep_dummy: retain dummy atoms
        :param mol: input molecule
        :return: a list of SMILES corresponding to individual ring systems
        """
        frag_list = Chem.GetMolFrags(mol, asMols=True)
        ring_system_list = []
        for frag in frag_list:
            if frag.HasSubstructMatch(self.ring_atom_pat):
                for atm in frag.GetAtoms():
                    if atm.GetAtomicNum() == 0:
                        if keep_dummy:
                            atm.SetProp("atomLabel", "R")
                        else:
                            atm.SetAtomicNum(1)
                        atm.SetIsotope(0)
                # Convert explict Hs to implicit
                frag = Chem.RemoveAllHs(frag)
                frag = self.fix_bond_stereo(frag)
                ring_system_list.append(frag)
        return ring_system_list

    @staticmethod
    def fix_bond_stereo(mol):
        """Loop over double bonds and change stereo specification for double bonds that don't have stereo
        :param mol: input RDKit molecule
        :return: output RDKit molecule
        """
        for bnd in mol.GetBonds():
            if bnd.GetBondType() == Chem.BondType.DOUBLE:
                begin_atm = bnd.GetBeginAtom()
                end_atm = bnd.GetEndAtom()
                # Look for double bond atoms with two attached hydrogens
                if begin_atm.GetDegree() == 1 or end_atm.GetDegree() == 1:
                    bnd.SetStereo(Chem.BondStereo.STEREONONE)
        return mol

    def find_ring_systems(self, mol, keep_dummy=False, as_mols=False):
        """Find the ring systems for an RDKit molecule
        :param as_mols: return results a molecules (otherwise return SMILES)
        :param keep_dummy: retain dummy atoms
        :param mol: input molecule
        :return: a list of SMILES corresponding to individual ring systems
        """
        self.tag_bonds_to_preserve(mol)
        frag_mol = self.cleave_linker_bonds(mol)
        output_list = self.cleanup_fragments(frag_mol, keep_dummy=keep_dummy)
        if not as_mols:
            output_list = [Chem.MolToSmiles(x) for x in output_list]
        return output_list
