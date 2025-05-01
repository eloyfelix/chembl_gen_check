from rdkit.Chem import rdFingerprintGenerator
from rdkit import Chem

MFPGEN = rdFingerprintGenerator.GetMorganGenerator(1)
ao = rdFingerprintGenerator.AdditionalOutput()
ao.AllocateAtomToBits()


def mol_to_pairs(mol):
    """
    function that fractures every bond and reports the two ECFP2
    (including dummy) at the fracture point.
    """
    id_pairs = []
    ri_full = mol.GetRingInfo()
    atom_rings = ri_full.AtomRings()
    for bond in mol.GetBonds():
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()

        # create a new molecule by fragmenting the current bond
        newmol = Chem.FragmentOnBonds(mol, [bond.GetIdx()])
        try:
            if bond.IsInRing():
                Chem.SanitizeMol(
                    newmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                )
                ri = newmol.GetRingInfo()
                for ring in atom_rings:
                    for idx in ring:
                        ri.AddRing((idx,), (0,))
            else:
                Chem.SanitizeMol(newmol)
            # sparse fingerprints for the fractured atoms
            MFPGEN.GetSparseFingerprint(
                newmol, fromAtoms=[begin_atom_idx, end_atom_idx], additionalOutput=ao
            )
            atom_to_bits = ao.GetAtomToBits()
            # take bits for radius 1 for atoms at both ends of the bond
            begin_fp = atom_to_bits[begin_atom_idx][1]
            end_fp = atom_to_bits[end_atom_idx][1]
            id_pairs.append(tuple(sorted([begin_fp, end_fp])))
        except Exception:
            pass
    return id_pairs


def assess_per_bond(mol, profile):
    pairs = mol_to_pairs(mol)
    total = profile["setsize"]
    idx = profile["idx"]
    pair_counts = profile["pairs"]
    results = []
    for pair in pairs:
        o1 = idx.get(pair[0], 0) / total / 2
        o2 = idx.get(pair[1], 0) / total / 2
        expected = o1 * o2
        real = pair_counts.get(pair, 0) / total
        results.append(0 if expected == 0 else real / expected)
    return results


def score_mol(mol, profile, t):
    apb = assess_per_bond(mol, profile)
    if not apb:
        apb = [0]
    min_val = min(apb)
    info = {"bad_bonds": [i for i, b in enumerate(apb) if b < t]}
    score = min(0.5 * (min_val / t) ** 0.5, 1.0)
    return score, info
