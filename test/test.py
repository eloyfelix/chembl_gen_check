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
