"""Fast PoseBusters filters without full energy evaluation."""

from copy import deepcopy

import numpy as np
import pandas as pd
import torch
from posebusters.modules.intermolecular_distance import _pairwise_distance
from rdkit import Chem
from rdkit.Chem import MolFromSmarts
from rdkit.Chem.rdchem import GetPeriodicTable, Mol
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix
from rdkit.Chem.rdmolops import SanitizeMol
from rdkit.Chem.rdShapeHelpers import ShapeTverskyIndex

_periodic_table = GetPeriodicTable()
# get all atoms from periodic table
atoms_vocab = {
    _periodic_table.GetElementSymbol(i + 1): i
    for i in range(_periodic_table.GetMaxAtomicNumber())
}
vdw_radius = torch.tensor(
    [
        _periodic_table.GetRvdw(_periodic_table.GetElementSymbol(i + 1))
        for i in range(_periodic_table.GetMaxAtomicNumber())
    ]
)

col_lb = "lower_bound"
col_ub = "upper_bound"
col_pe = "percent_error"
col_bpe = "bound_percent_error"
col_bape = "bound_absolute_percent_error"

bound_matrix_params = {
    "set15bounds": True,
    "scaleVDW": True,
    "doTriangleSmoothing": True,
    "useMacrocycle14config": False,
}

col_n_bonds = "number_bonds"
col_shortest_bond = "shortest_bond_relative_length"
col_longest_bond = "longest_bond_relative_length"
col_n_short_bonds = "number_short_outlier_bonds"
col_n_long_bonds = "number_long_outlier_bonds"
col_n_good_bonds = "number_valid_bonds"
col_bonds_result = "bond_lengths_within_bounds"
col_n_angles = "number_angles"
col_extremest_angle = "most_extreme_relative_angle"
col_n_bad_angles = "number_outlier_angles"
col_n_good_angles = "number_valid_angles"
col_angles_result = "bond_angles_within_bounds"
col_n_noncov = "number_noncov_pairs"
col_closest_noncov = "shortest_noncovalent_relative_distance"
col_n_clashes = "number_clashes"
col_n_good_noncov = "number_valid_noncov_pairs"
col_clash_result = "no_internal_clash"

_empty_results = {
    col_n_bonds: np.nan,
    col_shortest_bond: np.nan,
    col_longest_bond: np.nan,
    col_n_short_bonds: np.nan,
    col_n_long_bonds: np.nan,
    col_bonds_result: np.nan,
    col_n_angles: np.nan,
    col_extremest_angle: np.nan,
    col_n_bad_angles: np.nan,
    col_angles_result: np.nan,
    col_n_noncov: np.nan,
    col_closest_noncov: np.nan,
    col_n_clashes: np.nan,
    col_clash_result: np.nan,
}


# Allowable features for atom types (simplified version)
ALLOWABLE_ATOM_TYPES = [
    1,
    5,
    6,
    7,
    8,
    9,
    14,
    15,
    16,
    17,
    23,
    26,
    27,
    29,
    30,
    33,
    34,
    35,
    44,
    51,
    53,
    78,
]


def symmetrize_conjugated_terminal_bonds(df: pd.DataFrame, mol: Mol) -> pd.DataFrame:
    """
    Symmetrize the lower and upper bounds of the conjugated terminal bonds.

    Args:
        df: Dataframe with the bond geometry information and bounds.
        mol: RDKit molecule object

    Returns:
        Dataframe with symmetrized bounds for conjugated terminal bonds.
    """

    def _sort_bond_ids(bond_ids: tuple) -> tuple:
        return tuple(tuple(sorted(_)) for _ in bond_ids)

    def _get_terminal_group_matches(_mol: Mol) -> tuple:
        qsmarts = "[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]"
        qsmarts = MolFromSmarts(qsmarts)
        matches = _mol.GetSubstructMatches(qsmarts)
        return _sort_bond_ids(matches)

    df["atom_types_sorted"] = df["atom_types"].apply(
        lambda a: tuple(sorted(a.split("--")))
    )
    matches = _get_terminal_group_matches(mol)
    matched = df[df["atom_pair"].isin(matches)].copy()
    grouped = matched.groupby("atom_types_sorted").agg(
        {"lower_bound": np.amin, "upper_bound": np.amax}
    )
    index_orig = matched.index
    matched = matched.set_index("atom_types_sorted")
    matched.update(grouped)
    matched = matched.set_index(index_orig)
    df.update(matched)
    return df.drop(columns=["atom_types_sorted"])


def _get_bond_atom_indices(mol: Mol) -> list:
    bonds = []
    for bond in mol.GetBonds():
        bond_tuple = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        bond_tuple = _sort_bond(bond_tuple)
        bonds.append(bond_tuple)
    return bonds


def _get_angle_atom_indices(bonds: list) -> list:
    """Check all combinations of bonds to generate list of molecule angles."""
    angles = []
    bonds = list(bonds)
    for i in range(len(bonds)):
        for j in range(i + 1, len(bonds)):
            angle = _two_bonds_to_angle(bonds[i], bonds[j])
            if angle is not None:
                angles.append(angle)
    return angles


def _two_bonds_to_angle(bond1: tuple, bond2: tuple):
    set1 = set(bond1)
    set2 = set(bond2)
    all_atoms = set1 | set2
    if len(all_atoms) != 3:
        return None
    shared_atom = set1 & set2
    other_atoms = all_atoms - shared_atom
    return (min(other_atoms), shared_atom.pop(), max(other_atoms))


def _sort_bond(bond: tuple) -> tuple:
    return (min(bond), max(bond))


def _has_hydrogen(mol: Mol, idcs) -> bool:
    return any(_is_hydrogen(mol, idx) for idx in idcs)


def _is_hydrogen(mol: Mol, idx: int) -> bool:
    return mol.GetAtomWithIdx(int(idx)).GetAtomicNum() == 1


def mol_from_symbols_and_npcoords(symbols, coords_np: np.ndarray):
    """Create RDKit mol from symbols and coordinates.

    Args:
        symbols: List of atom symbols
        coords_np: Coordinates array shape (N, 3)

    Returns:
        RDKit Mol with conformer
    """
    assert coords_np.shape == (len(symbols), 3)
    rw = Chem.RWMol()
    for sym in symbols:
        a = Chem.Atom(sym)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        rw.AddAtom(a)
    m = rw.GetMol()
    conf = Chem.Conformer(len(symbols))
    conf.SetPositions(coords_np.astype(np.float64, copy=False))
    m.AddConformer(conf, assignId=True)
    return m


def check_intermolecular_distance(
    mol_orig,
    pos_pred,
    pos_cond,
    atom_names_pred,
    atom_names_cond,
    radius_type: str = "vdw",
    radius_scale: float = 1.0,
    clash_cutoff: float = 0.75,
    clash_cutoff_volume: float = 0.075,
    ignore_types: set | None = None,
    max_distance: float = 5.0,
    search_distance: float = 6.0,
    vdw_scale: float = 0.8,
):
    """Check that predicted molecule is not too close and not too far away from conditioning molecule.

    Args:
        mol_orig: Original ligand molecule
        pos_pred: Predicted ligand positions (n_preds, n_atoms, 3)
        pos_cond: Conditioning (protein) positions (n_atoms, 3)
        atom_names_pred: Ligand atom names
        atom_names_cond: Protein atom names
        radius_type: Type of atomic radius ("vdw" or "covalent")
        radius_scale: Scaling factor for radii
        clash_cutoff: Threshold for clash detection
        clash_cutoff_volume: Threshold for volume overlap
        ignore_types: Atom types to ignore
        max_distance: Maximum allowed distance
        search_distance: Search distance for nearby atoms
        vdw_scale: VDW radius scale for volume calculation

    Returns:
        Dictionary with filter results
    """
    if ignore_types is None:
        ignore_types = {"H"}
    device = "cuda" if torch.cuda.is_available() else "cpu"
    coords_ligand = torch.tensor(pos_pred, device=device).float()
    coords_protein = torch.tensor(pos_cond, device=device).float()

    atoms_ligand = torch.tensor(
        [atoms_vocab[atom] for atom in atom_names_pred], device=device
    ).long()
    atoms_protein_all = torch.tensor(
        [atoms_vocab[atom] for atom in atom_names_cond], device=device
    ).long()

    mask = atoms_ligand != atoms_vocab["H"]
    coords_ligand = coords_ligand[:, mask, :]
    atoms_ligand = atoms_ligand[mask]
    if ignore_types:
        mask = atoms_protein_all != atoms_vocab["H"]
        coords_protein = coords_protein[mask, :]
        atoms_protein_all = atoms_protein_all[mask]
        # Keep numpy arrays in sync with the filtered tensors so that boolean
        # masks derived from the tensors (e.g. ids_cond) can safely index them
        # in the ShapeTverskyIndex volume-overlap calculation below.
        mask_np = mask.cpu().numpy()
        atom_names_cond = atom_names_cond[mask_np]
        pos_cond = pos_cond[mask_np]

    radius_ligand = vdw_radius.to(device)[atoms_ligand]
    radius_protein_all = vdw_radius.to(device)[atoms_protein_all]

    distances_all = (coords_ligand[:, :, None] - coords_protein[None, None, :]).norm(
        dim=-1
    )
    distances = distances_all
    radius_protein = radius_protein_all

    is_buried_fraction = (distances < 5).any(dim=-1).sum(dim=-1) / distances.size(1)

    radius_sum = radius_ligand[None, :, None] + radius_protein[None, None, :]
    distance = distances
    sum_radii_scaled = radius_sum * radius_scale
    relative_distance = distance / sum_radii_scaled
    clash = relative_distance < clash_cutoff

    candidates = distance < (
        (radius_ligand[None, :, None] + radius_protein_all[None, None, :]) * vdw_scale
        + 2 * 3 * 0.25
    )
    ids_conds = candidates.any(dim=1).cpu().numpy()
    overlap = []
    for i in range(coords_ligand.size(0)):
        ids_cond = ids_conds[i]
        overlap.append(
            ShapeTverskyIndex(
                mol_from_symbols_and_npcoords(atom_names_pred, pos_pred[i]),
                mol_from_symbols_and_npcoords(
                    atom_names_cond[ids_cond], pos_cond[ids_cond]
                ),
                alpha=1,
                beta=0,
                vdwScale=vdw_scale,
            )
            < clash_cutoff_volume
        )

    results = {
        "not_too_far_away": (
            distance.reshape(distance.size(0), -1).min(dim=-1)[0] <= max_distance
        ).tolist(),
        "no_clashes": torch.logical_not(clash.any(dim=(1, 2))).tolist(),
        "no_volume_clash": overlap,
        "is_buried_fraction": is_buried_fraction.tolist(),
        "no_internal_clash": check_geometry(
            mol_orig,
            coords_ligand,
            threshold_bad_bond_length=0.25,
            threshold_clash=0.3,
            threshold_bad_angle=0.25,
            bound_matrix_params=bound_matrix_params,
            ignore_hydrogens=True,
            sanitize=True,
            symmetrize_conjugated_terminal_groups=True,
        ),
    }
    return {"results": results}


def check_volume_overlap(
    pos_pred,
    pos_cond,
    atom_names_pred,
    atom_names_cond,
    clash_cutoff: float = 0.05,
    vdw_scale: float = 0.8,
    ignore_types: set | None = None,
    search_distance: float = 6.0,
):
    """Check volume overlap between ligand and protein.

    Args:
        pos_pred: Predicted ligand positions
        pos_cond: Protein positions
        atom_names_pred: Ligand atom names
        atom_names_cond: Protein atom names
        clash_cutoff: Maximum allowed volume overlap fraction
        vdw_scale: VDW radius scale
        ignore_types: Atom types to ignore
        search_distance: Search distance for nearby atoms

    Returns:
        Dictionary with volume overlap results
    """
    if ignore_types is None:
        ignore_types = {"H"}
    keep_mask = atom_names_cond != "H"
    pos_cond = pos_cond[keep_mask]
    atom_names_cond = atom_names_cond[keep_mask]
    if len(pos_cond) == 0:
        return {"results": {"volume_overlap": np.nan, "no_volume_clash": True}}

    distances = _pairwise_distance(pos_pred, pos_cond)
    keep_mask = distances.min(axis=0) <= search_distance * vdw_scale
    pos_cond = pos_cond[keep_mask]
    atom_names_cond = atom_names_cond[keep_mask]
    if len(pos_cond) == 0:
        return {"results": {"volume_overlap": np.nan, "no_volume_clash": True}}

    ignore_hydrogens = "H" in ignore_types
    overlap = ShapeTverskyIndex(
        mol_from_symbols_and_npcoords(atom_names_pred, pos_pred),
        mol_from_symbols_and_npcoords(atom_names_cond, pos_cond),
        alpha=1,
        beta=0,
        vdwScale=vdw_scale,
        ignoreHs=ignore_hydrogens,
    )

    results = {
        "volume_overlap": overlap,
        "no_volume_clash": overlap <= clash_cutoff,
    }

    return {"results": results}


def check_geometry(
    mol_orig,
    pos_preds,
    threshold_bad_bond_length: float = 0.25,
    threshold_clash: float = 0.3,
    threshold_bad_angle: float = 0.25,
    bound_matrix_params=bound_matrix_params,
    ignore_hydrogens: bool = True,
    sanitize: bool = True,
    symmetrize_conjugated_terminal_groups: bool = True,
):
    """Use RDKit distance geometry bounds to check the geometry of a molecule.

    Args:
        mol_orig: Original molecule
        pos_preds: Predicted positions tensor
        threshold_bad_bond_length: Bond length threshold
        threshold_clash: Clash threshold
        threshold_bad_angle: Angle threshold
        bound_matrix_params: Parameters for GetMoleculeBoundsMatrix
        ignore_hydrogens: Whether to ignore hydrogens
        sanitize: Whether to sanitize molecule
        symmetrize_conjugated_terminal_groups: Whether to symmetrize terminal groups

    Returns:
        List of booleans indicating valid geometry
    """
    mol_pred = deepcopy(mol_orig)
    mol_pred.GetConformer().SetPositions(pos_preds[0].cpu().numpy().astype(np.float64))
    assert mol_pred.GetNumConformers() == 1, "Molecule must have exactly one conformer"
    mol = deepcopy(mol_pred)
    results = _empty_results.copy()
    if mol.GetNumConformers() == 0:
        print("Molecule does not have a conformer.")
        return {"results": results}
    if mol.GetNumAtoms() == 1:
        print(f"Molecule has only {mol.GetNumAtoms()} atoms.")
        results[col_angles_result] = True
        results[col_bonds_result] = True
        results[col_clash_result] = True
        return {"results": results}
    try:
        if sanitize:
            flags = SanitizeMol(mol)
            assert flags == 0, f"Sanitization failed with flags {flags}"
    except Exception:
        return {"results": results}
    bond_set = sorted(_get_bond_atom_indices(mol))
    angles = sorted(_get_angle_atom_indices(bond_set))
    angle_set = {(a[0], a[2]): a for a in angles}
    if len(bond_set) == 0:
        print("Molecule has no bonds.")

    bounds = GetMoleculeBoundsMatrix(mol, **bound_matrix_params)
    lower_triangle_idcs = np.tril_indices(mol.GetNumAtoms(), k=-1)
    upper_triangle_idcs = (lower_triangle_idcs[1], lower_triangle_idcs[0])
    df_12 = pd.DataFrame()
    df_12["atom_pair"] = list(zip(*upper_triangle_idcs, strict=False))
    df_12["atom_types"] = [
        "--".join(tuple(mol.GetAtomWithIdx(int(j)).GetSymbol() for j in i))
        for i in df_12["atom_pair"]
    ]
    df_12["angle"] = df_12["atom_pair"].apply(lambda x: angle_set.get(x))
    df_12["has_hydrogen"] = [_has_hydrogen(mol, i) for i in df_12["atom_pair"]]
    df_12["is_bond"] = [i in bond_set for i in df_12["atom_pair"]]
    df_12["is_angle"] = df_12["angle"].apply(lambda x: x is not None)
    df_12[col_lb] = bounds[lower_triangle_idcs]
    df_12[col_ub] = bounds[upper_triangle_idcs]
    if symmetrize_conjugated_terminal_groups:
        df_12 = symmetrize_conjugated_terminal_bonds(df_12, mol)
    distances_all = (pos_preds[:, :, None] - pos_preds[:, None, :]).norm(dim=-1)[
        :, lower_triangle_idcs[0], lower_triangle_idcs[1]
    ]
    distances_valid = distances_all[:, (~df_12["is_bond"] & ~df_12["is_angle"]).values]
    lower_bounds_valid = torch.tensor(
        df_12[col_lb][~df_12["is_bond"] & ~df_12["is_angle"]].values,
        device=distances_valid.device,
    )
    df_clash = torch.where(
        distances_valid >= lower_bounds_valid[None],
        0,
        (distances_valid - lower_bounds_valid[None]) / lower_bounds_valid[None],
    )
    col_n_clashes_count = (df_clash < -threshold_clash).sum(dim=-1)
    col_n_good_noncov_count = len(df_clash) - col_n_clashes_count
    res = (col_n_good_noncov_count == len(df_clash)).tolist()
    return res


def calc_posebusters(
    pos_pred, pos_cond, atom_ids_pred, atom_names_cond, names, lig_mol_for_posebusters
):
    """Calculate fast PoseBusters filters.

    Args:
        pos_pred: Predicted ligand positions
        pos_cond: Protein positions
        atom_ids_pred: Ligand atom type IDs
        atom_names_cond: Protein atom names
        names: Sample names (for error logging)
        lig_mol_for_posebusters: Ligand molecule for PoseBusters

    Returns:
        Dictionary with filter results or None on error
    """
    if 22 in atom_ids_pred:
        with open("error.txt", "a") as f:
            f.write(f"Error in {names}\n")
            f.write("22 (misc) in atom_ids_pred\n")
        return None
    atom_names_pred = np.array(
        [
            _periodic_table.GetElementSymbol(ALLOWABLE_ATOM_TYPES[atom_id])
            for atom_id in atom_ids_pred
            if atom_id >= 0
        ],
        dtype=object,
    )

    posebusters_results = {}
    try:
        assert len(pos_pred[0]) == len(atom_names_pred), (
            f"len(pos_pred[i]) = {len(pos_pred[0])} != len(atom_names_pred[i]) = {len(atom_names_pred)}"
        )
        assert len(pos_cond) == len(atom_names_cond), (
            f"len(pos_cond[i]) = {len(pos_cond[0])} != len(atom_names_cond[i]) = {len(atom_names_cond)}"
        )
    except Exception as e:
        print(f"Error in {names}")
        print(e)
        with open("error.txt", "a") as f:
            f.write(f"Error in {names}\n")
            f.write(
                f"len(pos_pred[i]) = {len(pos_pred[0])} != len(atom_names_pred[i]) = {len(atom_names_pred)}\n"
            )
        return None
    res1 = check_intermolecular_distance(
        lig_mol_for_posebusters, pos_pred, pos_cond, atom_names_pred, atom_names_cond
    )
    res = {**res1["results"]}
    for key in res:
        posebusters_results[key] = res[key]
    return posebusters_results
