"""Tests for RMSD computation."""

import numpy as np
import pytest


def test_rmsd_identical_coords():
    """Test RMSD is zero for identical coordinates."""
    from posebench_fast.metrics.rmsd import get_symmetry_rmsd_with_isomorphisms

    rng = np.random.default_rng(42)
    coords = rng.standard_normal((10, 3))
    isomorphisms = [(list(range(10)), list(range(10)))]

    rmsd = get_symmetry_rmsd_with_isomorphisms(coords, coords, isomorphisms)
    assert rmsd == pytest.approx(0.0, abs=1e-10)


def test_rmsd_translated_coords():
    """Test RMSD for translated coordinates."""
    from posebench_fast.metrics.rmsd import get_symmetry_rmsd_with_isomorphisms

    coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    coords2 = coords1 + np.array([1, 0, 0])  # translate by 1 in x
    isomorphisms = [(list(range(3)), list(range(3)))]

    rmsd = get_symmetry_rmsd_with_isomorphisms(coords1, coords2, isomorphisms)
    assert rmsd == pytest.approx(1.0, abs=1e-10)


def test_compute_isomorphisms():
    """Test that isomorphisms are computed without error."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from posebench_fast.metrics.rmsd import compute_all_isomorphisms

    # Simple molecule: methane with 3D coords
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    isomorphisms = compute_all_isomorphisms(mol)
    assert len(isomorphisms) > 0
    assert isinstance(isomorphisms[0], tuple)


def test_timeout_exception():
    """Test TimeoutException can be raised and caught."""
    import time

    from posebench_fast.metrics.rmsd import TimeoutException, time_limit

    with pytest.raises(TimeoutException), time_limit(1):
        time.sleep(2)
