# posebench-fast

[![PyPI version](https://badge.fury.io/py/posebench-fast.svg)](https://pypi.org/project/posebench-fast/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Fast docking evaluation metrics: symmetry-corrected RMSD and lightweight PoseBusters filters.

## Installation

```bash
uv pip install posebench-fast
```

## Features

- **Symmetry-corrected RMSD** — accounts for molecular symmetry (benzene, carboxylates, etc.)
- **Fast PoseBusters filters** — 4 key checks in ~10ms instead of full 27-test suite (~1-2s)

## Usage

```python
from posebench_fast import compute_all_isomorphisms, get_symmetry_rmsd_with_isomorphisms

# Symmetry-corrected RMSD
isomorphisms = compute_all_isomorphisms(rdkit_mol)
rmsd = get_symmetry_rmsd_with_isomorphisms(true_coords, pred_coords, isomorphisms)
```

```python
from posebench_fast import check_intermolecular_distance

# Fast filters: not_too_far_away, no_clashes, no_volume_clash, no_internal_clash
results = check_intermolecular_distance(
    mol_orig=rdkit_mol,
    pos_pred=pred_positions,      # (n_samples, n_atoms, 3)
    pos_cond=protein_positions,   # (n_protein_atoms, 3)
    atom_names_pred=lig_atoms,
    atom_names_cond=prot_atoms,
)
```

## Related

- [PoseBench](https://github.com/BioinfoMachineLearning/PoseBench) — full benchmark suite
- [PoseBusters](https://github.com/maabuu/posebusters) — full 27-test validation
- [spyrmsd](https://github.com/RMeli/spyrmsd) — symmetry RMSD algorithms

## License

MIT
