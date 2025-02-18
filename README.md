# GenVal

GenVal is a simple tool to validate molecules produced by generative models. It assesses whether the scaffolds or ring systems of these molecules are present in the ChEMBL database.

Leveraging [molbloom](https://github.com/whitead/molbloom) filters for ChEMBL scaffolds and ring systems, GenVal enables rapid determination of whether a molecule's scaffold or ring systems exist within the database, maintaining a low false positive rate of 0.000025.

```python
from genval import GenVal

gv = GenVal()
gv.check_scaffold("CC(=O)Oc1ccccc1C(=O)O")
gv.check_ring_systems("CC(=O)Oc1ccccc1C(=O)O")
```


Code to extract ring systems comes from Patt Walters' [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py)
