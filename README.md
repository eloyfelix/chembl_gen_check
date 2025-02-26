# GenVal

GenVal is a simple tool to validate molecules produced by generative models. It assesses whether the scaffolds or ring systems of these molecules are present in the ChEMBL database.

Leveraging [molbloom](https://github.com/whitead/molbloom) filters for ChEMBL scaffolds and ring systems, GenVal enables rapid determination of whether a molecule's scaffold or ring systems exist within the database, maintaining a low false positive rate of 0.000025.

```python
from genval import GenVal

gv = GenVal()

smiles = "CC(=O)Oc1ccccc1C(=O)O"

gv.check_scaffold(smiles) # Using scaffolds found in ChEMBL
gv.check_ring_systems(smiles) # Using ring systems found in ChEMBL
gv.check_lacan(smiles) # Profile generated using ChEMBL
gv.check_structural_alerts(smiles) # ChEMBL set
```

Code to extract ring systems adapted from: W Patrick Walters. [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py)

Code to calculate LACAN scores adapted from: Dehaen, W. LACAN. https://github.com/dehaenw/lacan/