import os, sys
from collections import Counter
import pytest
from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
from chembl_gen_check import Checker
from chembl_gen_check.lacan import mol_to_pairs


def build_lacan_profile(smiles_list):
    idx_counter = Counter()
    pair_counter = Counter()
    total_pairs = 0

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        pairs = mol_to_pairs(mol)
        pair_counter.update(pairs)
        for a, b in pairs:
            idx_counter[a] += 1
            idx_counter[b] += 1
        total_pairs += len(pairs)

    return {
        "idx": dict(idx_counter),
        "pairs": dict(pair_counter),
        "setsize": total_pairs,
        "atom_invariant_scheme": "ring_type",
        "score_formula": "pmi_over_one_plus_pmi",
    }


# A molecule whose bonds define the reference profile, and one with bond
# environments absent from that profile (so at least one bond scores ~0).
LACAN_PROFILED_SMILES = "COc1ccc2[nH]cc(CNC(C)=O)c2c1"
LACAN_UNUSUAL_SMILES = "FNCCC(c1cc2OCOc2cc1)c1ccccc1"


def _lacan_checker():
    checker = Checker()
    checker.lacan_profile = build_lacan_profile([LACAN_PROFILED_SMILES])
    return checker


def test_check_lacan_threshold_mode_passes_profiled_and_fails_unusual():
    """threshold mode returns a hard 1.0 (pass) / 0.0 (fail) verdict."""
    checker = _lacan_checker()

    checker.load_smiles(LACAN_PROFILED_SMILES)
    good_score, good_info = checker.check_lacan(mode="threshold", include_info=True)

    checker.load_smiles(LACAN_UNUSUAL_SMILES)
    bad_score, bad_info = checker.check_lacan(mode="threshold", include_info=True)

    assert good_score == 1.0
    assert bad_score == 0.0
    assert good_info["bad_bonds"] == []
    assert bad_info["bad_bonds"]


def test_check_lacan_score_mode_ranks_profiled_bonds_higher():
    """score mode returns a continuous min_PMI / (1 + min_PMI) value in [0, 1)."""
    checker = _lacan_checker()

    checker.load_smiles(LACAN_PROFILED_SMILES)
    good_score, good_info = checker.check_lacan(mode="score", include_info=True)

    checker.load_smiles(LACAN_UNUSUAL_SMILES)
    bad_score, bad_info = checker.check_lacan(mode="score", include_info=True)

    assert 0.0 < good_score < 1.0
    assert bad_score == 0.0
    assert good_score > bad_score
    assert good_info["bad_bonds"] == []
    assert bad_info["bad_bonds"]


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2ccccc2c1", True),
        ("C1=Cc2ccccc2[U]=C1", False),
    ],
)
def test_scaffold_filter(smiles, expected):
    checker = Checker()
    checker.load_smiles(smiles)
    assert checker.check_scaffold() is expected, (
        f"Expected check_scaffold to return {expected} for {smiles}"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2ccccc2c1", True),
        ("C1=Cc2ccccc2[U]=C1", True),
    ],
)
def test_skeleton_filter(smiles, expected):
    checker = Checker()
    checker.load_smiles(smiles)
    assert checker.check_skeleton() is expected, (
        f"Expected check_skeleton to return {expected} for {smiles}"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2ccccc2c1", True),
        ("C1=Cc2ccccc2[U]=C1", False),
    ],
)
def test_ring_system_filter(smiles, expected):
    checker = Checker()
    checker.load_smiles(smiles)
    assert checker.check_ring_systems() is expected, (
        f"Expected check_ring_systems to return {expected} for {smiles}"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Nc1ccc(O)cc1", 1),
        ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", 0),
    ],
)
def test_structural_alerts(smiles, expected):
    """Test that check_structural_alerts returns expected counts for given SMILES."""
    checker = Checker()
    checker.load_smiles(smiles)
    result = checker.check_structural_alerts()
    assert result == expected, f"Expected {expected} alerts for {smiles}, got {result}"


def test_check_all_returns_dict_with_expected_keys():
    """Test that check_all() returns a dict with all expected keys and types."""
    checker = Checker()
    checker.load_smiles("c1ccc2ccccc2c1")  # naphthalene
    result = checker.check_all()
    
    assert isinstance(result, dict), "check_all() should return a dict"
    expected_keys = {"scaffold", "skeleton", "ring_systems", "structural_alerts", "lacan"}
    assert set(result.keys()) == expected_keys, f"Expected keys {expected_keys}, got {set(result.keys())}"
    
    # Check types
    assert isinstance(result["scaffold"], bool)
    assert isinstance(result["skeleton"], bool)
    assert isinstance(result["ring_systems"], bool)
    assert isinstance(result["structural_alerts"], int)
    assert isinstance(result["lacan"], float)


def test_check_all_consistency_with_individual_checks():
    """Test that check_all() results match individual check methods."""
    checker = Checker()
    checker.load_smiles("COc1ccc2[nH]cc(CNC(C)=O)c2c1")
    
    all_results = checker.check_all()
    
    # Results should match individual checks (note: scaffold/skeleton cached)
    assert all_results["scaffold"] == checker.check_scaffold()
    assert all_results["skeleton"] == checker.check_skeleton()
    assert all_results["ring_systems"] == checker.check_ring_systems()
    assert all_results["structural_alerts"] == checker.check_structural_alerts()
    assert all_results["lacan"] == checker.check_lacan()


def test_surechembl_checker_loads():
    """Test that Checker with surechembl database loads successfully."""
    checker = Checker(db_name="surechembl")
    checker.load_smiles("c1ccccc1")  # benzene
    
    # Should be able to run all checks without error
    result = checker.check_all()
    assert isinstance(result, dict)
    assert len(result) == 5
