import os, sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
from chembl_gen_check import Checker


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COc1ccc2[nH]cc(CNC(C)=O)c2c1", 1.0),
        ("FNCCC(c1cc2OCOc2cc1)c1ccccc1", 0.0),
    ],
)
def test_check_lacan(smiles, expected):
    checker = Checker()
    checker.load_smiles(smiles)
    assert checker.check_lacan() == expected, (
        f"Expected check_lacan to return {expected} for {smiles}"
    )


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
