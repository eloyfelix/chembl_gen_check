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
    - smallest ring atom is in. set to 0 if not in ring
    """
    invs = []
    sssr = Chem.GetSSSR(mol)
    min_ring = {}
    for ring in sssr:
        for a in ring:
            if a not in min_ring:
                min_ring[a] = len(ring)
            else:
                if min_ring[a] > len(ring):
                    min_ring[a] = len(ring)
    for idx, a in enumerate(mol.GetAtoms()):
        inv = []
        inv.append(a.GetAtomicNum())
        inv.append(a.GetDegree())
        inv.append(a.GetNumExplicitHs() + a.GetNumImplicitHs())
        inv.append(a.GetFormalCharge())
        try:
            inv.append(min_ring[idx])
        except:
            inv.append(0)
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


def score_mol(mol, profile, t):
    apb = assess_per_bond(mol, profile)
    if not apb:
        apb = [0]
    min_val = min(apb)
    info = {"bad_bonds": [i for i, b in enumerate(apb) if b < t]}
    score = min(0.5 * (min_val / t) ** 0.5, 1.0)
    return score, info
