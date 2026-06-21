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


def score_mol(mol, profile, t=0.05, mode="threshold"):
    apb = assess_per_bond(mol, profile)

    protected = {
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.HasProp("_lp") and bond.GetBoolProp("_lp")
    }

    # Pair each per-bond score with its originating bond so that protection
    # filtering and bad-bond reporting use real bond indices instead of relying
    # on the score list being in bond-index order. assess_per_bond/mol_to_pairs
    # iterate mol.GetBonds() in order, so zipping against it realigns scores
    # with the bonds they were computed from.
    scored_bonds = [
        (bond.GetIdx(), score)
        for bond, score in zip(mol.GetBonds(), apb)
        if bond.GetIdx() not in protected
    ]

    if not scored_bonds:
        if apb:
            # Every bond is protected: nothing left to penalize.
            return 1.0, {"bad_bonds": []}
        # No bonds at all (e.g. single-atom input): no precedent to establish.
        min_val = 0.0
        bad_bonds = []
    else:
        min_val = min(score for _, score in scored_bonds)
        bad_bonds = [idx for idx, score in scored_bonds if score < t]

    info = {"bad_bonds": bad_bonds}

    if mode == "threshold":
        score = 0.0 if min_val < t else 1.0
    elif mode == "score":
        score = min_val / (1 + min_val)
    else:
        raise ValueError(f"Unsupported LACAN scoring mode: {mode}")

    return score, info
