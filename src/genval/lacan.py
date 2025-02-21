from rdkit.Chem import rdFingerprintGenerator
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw

MFPGEN = rdFingerprintGenerator.GetMorganGenerator(1)
ao = rdFingerprintGenerator.AdditionalOutput()
ao.AllocateBitInfoMap()
ao.AllocateAtomToBits()


def mol_to_pairs(m):
    """
    function that fractures every bond and reports the two ECFP2
    (including dummy) at the fracture point.
    """
    id_pairs = []
    ri_full = m.GetRingInfo()
    ar = ri_full.AtomRings()
    for b in m.GetBonds():
        bidx = [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
        newmol = Chem.FragmentOnBonds(m, [b.GetIdx()])
        try:
            if b.IsInRing():
                Chem.SanitizeMol(
                    newmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                )
                ri = newmol.GetRingInfo()  # reset ringinfo trick
                for ring in ar:
                    for idx in ring:
                        ri.AddRing((idx,), (0,))
            else:
                Chem.SanitizeMol(newmol)
            MFPGEN.GetSparseFingerprint(newmol, fromAtoms=bidx, additionalOutput=ao)
            id_pairs.append(tuple(sorted([ao.GetAtomToBits()[idx][1] for idx in bidx])))
        except Exception as e:
            pass  # silent for now
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


def score_mol(mol, profile, mode="score", t=0.05):
    apb = assess_per_bond(mol, profile)
    if not apb:
        apb = [0]
    min_val = min(apb)
    info = {"bad_bonds": [i for i, b in enumerate(apb) if b < t]}

    if mode == "threshold":
        score = 0 if min_val < t else 1
    elif mode == "score":
        score = min(0.5 * (min_val / t) ** 0.5, 1.0)
    else:
        print("mode not supported yet, sorry.")
        score = None
    return score, info


def highlight_bonds_svg(mol, bond_indices, size=(300, 300)):
    if not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size[0], size[1])

    colors = {i: (1, 0, 0) for i in bond_indices}

    drawer.DrawMolecule(
        mol,
        highlightAtoms=[],
        highlightBonds=bond_indices,
        highlightAtomColors={},
        highlightBondColors=colors,
    )
    drawer.FinishDrawing()

    return drawer.GetDrawingText()
