# posebench-fast

Fast docking evaluation metrics for molecular docking benchmarks.

## Features

- **Symmetry-corrected RMSD**: Compute RMSD accounting for molecular symmetry using graph isomorphisms
- **Fast PoseBusters filters**: Quick physical validity checks without full energy evaluation
  - Distance checks (not too far, no clashes)
  - Volume overlap detection
  - Internal geometry validation
- **Metrics aggregation**: Success rates, averages, filtering by scores

## Installation

```bash
uv add posebench-fast
```

Or with pip:
```bash
pip install posebench-fast
```

## Usage

### Symmetry RMSD

```python
from posebench_fast import compute_all_isomorphisms, get_symmetry_rmsd_with_isomorphisms
import numpy as np

# Compute isomorphisms once per molecule
isomorphisms = compute_all_isomorphisms(rdkit_mol)

# Compute symmetry-corrected RMSD
rmsd = get_symmetry_rmsd_with_isomorphisms(true_coords, pred_coords, isomorphisms)
```

### Fast PoseBusters Filters

```python
from posebench_fast import calc_posebusters

results = calc_posebusters(
    pos_pred=ligand_positions,      # (n_samples, n_atoms, 3)
    pos_cond=protein_positions,     # (n_atoms, 3)
    atom_ids_pred=ligand_atom_ids,
    atom_names_cond=protein_atom_names,
    names="sample_id",
    lig_mol_for_posebusters=rdkit_mol
)

# Results contain:
# - not_too_far_away: bool list
# - no_clashes: bool list
# - no_volume_clash: bool list
# - no_internal_clash: bool list
# - is_buried_fraction: float list
```

### Metrics Aggregation

```python
from posebench_fast import get_final_results_for_df

rows, scored_results = get_final_results_for_df(
    full_results,
    score_names=['error_estimate_0'],
    posebusters_filter=True,
    fast_filter=True
)
# Returns DataFrame with: RMSD < 2A, RMSD < 5A, SymRMSD < 2A, etc.
```

## Dependencies

- numpy
- pandas
- spyrmsd
- rdkit
- torch
- posebusters
- tqdm

## License

MIT
