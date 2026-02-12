"""Tests for fast_filters — intermolecular distance checks."""

import numpy as np


def test_check_intermolecular_distance_with_protein_hydrogens():
    """Regression: protein H atoms must be stripped from numpy arrays too.

    When ``ignore_types={"H"}``, the function filters hydrogen atoms from the
    torch tensors (``coords_protein``, ``atoms_protein_all``) but previously
    did NOT filter the corresponding numpy arrays (``atom_names_cond``,
    ``pos_cond``).  A downstream boolean mask (``ids_cond``) produced by the
    filtered tensor was then used to index the *unfiltered* numpy arrays,
    causing a shape mismatch / ``IndexError`` whenever the protein PDB
    contained explicit hydrogen atoms.
    """
    from posebench_fast.filters.fast_filters import check_intermolecular_distance

    rng = np.random.default_rng(42)

    # --- Ligand: 5 heavy atoms, no H ----------------------------------------
    n_lig = 5
    atom_names_pred = np.array(["C", "C", "N", "O", "C"])
    pos_pred = rng.standard_normal((1, n_lig, 3)).astype(np.float32)

    # Minimal RDKit mol for check_geometry (ethanol — 5 heavy atoms after RemoveHs)
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCNOC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol = Chem.RemoveHs(mol)

    # --- Protein: mix of heavy atoms and explicit H --------------------------
    heavy_names = ["C", "N", "O", "C", "N"]
    h_names = ["H", "H", "H"]
    atom_names_cond = np.array(heavy_names + h_names)
    n_prot = len(atom_names_cond)  # 8 total (5 heavy + 3 H)
    # Place protein atoms close to the ligand so it's "not too far away"
    pos_cond = (
        pos_pred[0, 0] + rng.standard_normal((n_prot, 3)).astype(np.float32) * 2.0
    )

    # This call should NOT raise IndexError regardless of H filtering.
    result = check_intermolecular_distance(
        mol_orig=mol,
        pos_pred=pos_pred,
        pos_cond=pos_cond,
        atom_names_pred=atom_names_pred,
        atom_names_cond=atom_names_cond,
        ignore_types={"H"},
        clash_cutoff=0.75,
        clash_cutoff_volume=0.075,
        max_distance=5.0,
    )

    # Basic structure check
    assert "results" in result
    for key in (
        "no_clashes",
        "no_volume_clash",
        "not_too_far_away",
        "no_internal_clash",
    ):
        assert key in result["results"], f"Missing key: {key}"
        assert len(result["results"][key]) == 1  # one prediction


def test_check_intermolecular_distance_no_hydrogens():
    """Sanity check: works correctly when protein has zero H atoms."""
    from posebench_fast.filters.fast_filters import check_intermolecular_distance

    rng = np.random.default_rng(123)

    n_lig = 3
    atom_names_pred = np.array(["C", "N", "O"])
    pos_pred = rng.standard_normal((1, n_lig, 3)).astype(np.float32)

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CNO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol = Chem.RemoveHs(mol)

    # Protein with only heavy atoms
    atom_names_cond = np.array(["C", "N", "O", "S"])
    n_prot = len(atom_names_cond)
    pos_cond = (
        pos_pred[0, 0] + rng.standard_normal((n_prot, 3)).astype(np.float32) * 2.0
    )

    result = check_intermolecular_distance(
        mol_orig=mol,
        pos_pred=pos_pred,
        pos_cond=pos_cond,
        atom_names_pred=atom_names_pred,
        atom_names_cond=atom_names_cond,
        ignore_types={"H"},
    )

    assert "results" in result
    for key in ("no_clashes", "no_volume_clash", "not_too_far_away"):
        assert len(result["results"][key]) == 1
