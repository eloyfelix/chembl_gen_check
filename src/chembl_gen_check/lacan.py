from rdkit import Chem
from itertools import chain
import hashlib
import copy


def get_neighbors(mol):
    """
    return the idx of each atoms directly bounds neighbots
    """
    return [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]


def get_atom_invariants(mol):
    """
    get ECFP like atom identifiers. these are
    - atom number
    - degree
    - h count
    - formal charge
    - ring type:
    - 0 if acyclic
    - 1 if in a non-aromatic ring
    - 2 if in an aromatic ring
    """
    invs = []
    for a in mol.GetAtoms():
        inv = [
            a.GetAtomicNum(),
            a.GetDegree(),
            a.GetNumExplicitHs() + a.GetNumImplicitHs(),
            a.GetFormalCharge(),
            int(a.IsInRing()) + int(a.GetIsAromatic()),
        ]
        invs.append(inv)
    return invs


def hash_invariants(invs):
    """
    md5 hash folded to 32bit int to have
    reproducible and portable hashes for environments
    """
    h = hashlib.md5()
    h.update(str(invs).encode())
    return int.from_bytes(h.digest()[:4], "big", signed=True)  # 32 bit prefix


def mol_to_pairs(mol):
    """
    function that fractures every bond and reports the two ECFP2like
    identifiers at the fracture point.
    New code calculates them explicitly instead of using fpgenetators
    from rdkit. This is to avoid the time sink of bond fracturing and
    aromatic sanitization.
    """
    if mol:
        nb = get_neighbors(mol)
        a_invs = get_atom_invariants(mol)
        invs = []
        for b in mol.GetBonds():
            b1 = b.GetBeginAtomIdx()
            b2 = b.GetEndAtomIdx()
            nb1 = copy.copy(nb[b1])
            nb2 = copy.copy(nb[b2])
            nb1.remove(b2)
            nb2.remove(b1)
            bt = int(b.GetBondType())
            n_inv1 = [a_invs[b1] + [bt]] + sorted(
                [
                    a_invs[n] + [int(mol.GetBondBetweenAtoms(b1, n).GetBondType())]
                    for n in nb1
                ]
            )
            n_inv2 = [a_invs[b2] + [bt]] + sorted(
                [
                    a_invs[n] + [int(mol.GetBondBetweenAtoms(b2, n).GetBondType())]
                    for n in nb2
                ]
            )
            h1 = hash_invariants(list(chain.from_iterable(n_inv1)))
            h2 = hash_invariants(list(chain.from_iterable(n_inv2)))
            invs.append(tuple(sorted([h1, h2])))
        return invs
    else:
        print("there was a molecule that didn't parse.")
        return []


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


def score_mol(mol, profile, t=0.05, mode="score"):
    apb = assess_per_bond(mol, profile)
    if not apb:
        apb = [0.0]

    protected = {
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.HasProp("_lp") and bond.GetBoolProp("_lp")
    }
    apb_active = [score for i, score in enumerate(apb) if i not in protected]
    if not apb_active:
        return 1.0, {"bad_bonds": []}

    min_val = min(apb_active)
    info = {"bad_bonds": [i for i, b in enumerate(apb) if b < t and i not in protected]}

    if mode == "threshold":
        score = 0.0 if min_val < t else 1.0
    elif mode == "score":
        score = min_val / (1 + min_val)
    else:
        raise ValueError(f"Unsupported LACAN scoring mode: {mode}")

    return score, info
