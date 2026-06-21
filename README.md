# chembl_gen_check

`chembl_gen_check` is a Python library for rapidly assessing how reasonable a (generated) molecule is. It ships with precomputed databases for both [ChEMBL](https://www.ebi.ac.uk/chembl/) (medicinal chemistry from the scientific literature) and [SureChEMBL](https://www.surechembl.org/) (chemical structures extracted from patent literature), so no downloads or preprocessing are required: just pick one when loading your checker. 

Using lightweight [MolBloom](https://github.com/whitead/molbloom) filters, it verifies whether a compound's scaffolds, generic scaffolds or ring systems already exist in the selected database. It can also flag uncommon bonds via the [LACAN](https://github.com/dehaenw/lacan) algorithm and report structural alerts. Taken together, these checks give a fast read on the plausibility of a molecule's ring systems and scaffolds, and ensure that its atom and bond environments have precedent.

## Installation

```
pip install chembl-gen-check
```

## Usage example

```python
from chembl_gen_check import Checker

checker = Checker("chembl")
#checker = Checker("surechembl")

smiles = "CCN(CC)C(=O)C[C@H]1C[C@@H]1c1ccccc1"
checker.load_smiles(smiles)

# Murcko scaffold found in the loaded database (True/False)
checker.check_scaffold()

# Generic Murcko scaffold found in loaded database (True/False)
checker.check_skeleton()

# All molecule ring systems found in loaded database (True/False)
checker.check_ring_systems()

# Number of structural alerts using the ChEMBL set in RDKit(integer)
checker.check_structural_alerts()

# LACAN hard pass/fail filter (default mode="threshold"): reject if any bond's
# PMI is below the threshold t (default 0.05). Returns 1.0 (pass) or 0.0 (fail).
checker.check_lacan()
```

## How LACAN Works

Reference: Dehaen, W. LACAN. [ChemRxiv preprint](https://chemrxiv.org/doi/full/10.26434/chemrxiv.15001196/v1).

LACAN scores a molecule one bond at a time. Each bond is split into its two
atom environments, and a PMI ratio (pointwise mutual information) is computed
from the reference database:

```
PMI = P(env_a, env_b) / (P(env_a) * P(env_b))
```

that is, how often the two halves are actually bonded together versus how often
they would be by pure chance. `PMI > 1` means the bond is more common than
chance, `PMI ≈ 1` is as expected, and `PMI ≈ 0` flags a junction that is
essentially never seen in the database (a likely artifact).

A molecule is summarized by its **weakest** bond, `min_PMI`:

- **`mode="threshold"`**(default): a single bond is enough to fail the whole molecule —
  if any bond has `PMI < t` (default `0.05`) it returns `0.0` (fail), otherwise
  `1.0` (pass).
- **`mode="score"`** : returns `min_PMI / (1 + min_PMI)`, a value in
  `[0, 1)`. Higher is more reasonable; `0.5` corresponds to `min_PMI = 1`.

Code to extract ring systems adapted from: W Patrick Walters. [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py)

Code to calculate LACAN scores adapted from: Dehaen, W. LACAN. https://github.com/dehaenw/lacan/
