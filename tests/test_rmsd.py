import copy
import math

import numpy as np
import pytest
from Bio.PDB import Atom, Residue

from scripts.core import squared_diffs_between_sidechain_heavy_atoms
from scripts.plddt import per_residue_plddt
from scripts.rmsd import global_rmsd, per_residue_rmsd, per_residue_sidechain_rmsd


def test_per_residue_rmsd_self_alignment_is_zero(sample_chain):
    pred_chain = copy.deepcopy(sample_chain)

    result = per_residue_rmsd("toy", "AGS", "AGS", sample_chain, pred_chain)

    assert result.keys() == {0, 1, 2}
    assert all(np.isclose(value, 0.0) for value in result.values())


def test_global_rmsd_self_alignment_is_zero(sample_chain):
    pred_chain = copy.deepcopy(sample_chain)

    assert np.isclose(global_rmsd("toy", "AGS", "AGS", sample_chain, pred_chain), 0.0)


def test_per_residue_plddt_reads_ca_bfactors(make_chain):
    ref_chain = make_chain("AGS")
    pred_chain = make_chain("AGS", ca_bfactors=[91.5, 82.0, 70.25])

    result = per_residue_plddt("toy", "AGS", "AGS", ref_chain, pred_chain)

    assert result == {0: 91.5, 1: 82.0, 2: 70.25}


# ── Unit tests: squared_diffs_between_sidechain_heavy_atoms ──────────────────

def _make_residue(resname: str, seq_id: int, atoms: dict) -> Residue.Residue:
    """Build a minimal residue with the given atom name → coordinate mapping."""
    res = Residue.Residue((" ", seq_id, " "), resname, " ")
    for serial, (name, coord) in enumerate(atoms.items(), start=1):
        element = name[0] if name[0].isalpha() and name[0] != "H" else "H"
        res.add(Atom.Atom(name, np.array(coord, dtype=float), 0.0, 1.0, " ", name, serial, element))
    return res


def test_sc_diffs_backbone_only_returns_empty():
    """Residues with only backbone atoms (N, CA, C, O) yield an empty list.

    This mirrors the ESM3-Open / Ember3D case where backbone-only PDBs are
    provided: no sidechain heavy atoms → nothing to compare → [].
    """
    backbone = {"N": [0, 0, 0], "CA": [1, 0, 0], "C": [2, 0, 0], "O": [2, 1, 0]}
    r_ref = _make_residue("ALA", 1, backbone)
    r_pred = _make_residue("ALA", 1, backbone)

    assert squared_diffs_between_sidechain_heavy_atoms(r_ref, r_pred) == []


def test_sc_diffs_glycine_has_no_sidechain():
    """GLY has no CB, so sidechain diffs are always empty."""
    backbone = {"N": [0, 0, 0], "CA": [1, 0, 0], "C": [2, 0, 0], "O": [2, 1, 0]}
    gly_ref = _make_residue("GLY", 1, backbone)
    gly_pred = _make_residue("GLY", 1, backbone)

    assert squared_diffs_between_sidechain_heavy_atoms(gly_ref, gly_pred) == []


def test_sc_diffs_cb_displacement_gives_expected_squared_diff():
    """Single CB displaced by (3, 4, 0) → exactly one squared diff equal to 25."""
    r_ref = _make_residue("ALA", 1, {"CA": [1, 0, 0], "CB": [0, 0, 0]})
    r_pred = _make_residue("ALA", 1, {"CA": [1, 0, 0], "CB": [3, 4, 0]})

    diffs = squared_diffs_between_sidechain_heavy_atoms(r_ref, r_pred)

    assert len(diffs) == 1
    assert math.isclose(diffs[0], 25.0)


def test_sc_diffs_excludes_backbone_atoms():
    """Backbone atoms (N, C, O, OXT) must not appear in the squared diffs even
    when their coordinates differ between ref and pred."""
    r_ref = _make_residue("ALA", 1, {
        "N": [0, 0, 0], "CA": [1, 0, 0], "C": [2, 0, 0], "O": [2, 1, 0], "CB": [1, 1, 0],
    })
    r_pred = _make_residue("ALA", 1, {
        "N": [9, 9, 9], "CA": [1, 0, 0], "C": [9, 9, 9], "O": [9, 9, 9], "CB": [1, 1, 0],
    })

    diffs = squared_diffs_between_sidechain_heavy_atoms(r_ref, r_pred)

    # Only CB is sidechain; CB positions are identical → single zero diff.
    assert len(diffs) == 1
    assert math.isclose(diffs[0], 0.0)


# ── Integration tests: per_residue_sidechain_rmsd ────────────────────────────

def test_per_residue_sc_rmsd_backbone_only_all_nan(sample_chain):
    """Chains with no sidechain atoms → every position is NaN.

    sample_chain contains only backbone atoms (N, CA, C, O), matching the
    structure of ESM3-Open and Ember3D predictions.
    """
    result = per_residue_sidechain_rmsd(
        "toy", "AGS", "AGS", sample_chain, copy.deepcopy(sample_chain)
    )

    assert all(math.isnan(v) for v in result.values())


def test_per_residue_sc_rmsd_self_is_zero(make_chain_with_cb):
    """Identical chains with sidechain atoms → scRMSD is 0.0 at every position."""
    ref = make_chain_with_cb("AAS")
    pred = copy.deepcopy(ref)

    result = per_residue_sidechain_rmsd("toy", "AAS", "AAS", ref, pred)

    assert result.keys() == {0, 1, 2}
    assert all(math.isclose(v, 0.0) for v in result.values())


def test_per_residue_sc_rmsd_glycine_is_nan(make_chain_with_cb):
    """GLY (no CB) produces NaN; flanking residues with a CB give 0.0."""
    ref = make_chain_with_cb("AGS")
    pred = copy.deepcopy(ref)

    result = per_residue_sidechain_rmsd("toy", "AGS", "AGS", ref, pred)

    assert math.isnan(result[1]), "GLY position should be NaN"
    assert math.isclose(result[0], 0.0), "ALA position should be 0.0"
    assert math.isclose(result[2], 0.0), "SER position should be 0.0"


def test_per_residue_sc_rmsd_known_cb_displacement(make_chain_with_cb):
    """CB displaced by (3, 4, 0) from the aligned reference → scRMSD = 5.0.

    Backbone positions (including CA) are identical in ref and pred, so the
    superposition is the identity transform and the CB displacement is preserved
    exactly after alignment.
    """
    ref  = make_chain_with_cb("A", cb_coords=[(0.0, 0.0, 0.0)])
    pred = make_chain_with_cb("A", cb_coords=[(3.0, 4.0, 0.0)])

    result = per_residue_sidechain_rmsd("toy", "A", "A", ref, pred)

    assert math.isclose(result[0], 5.0, abs_tol=1e-6)


def test_per_residue_sc_rmsd_backbone_atoms_not_included(make_chain_with_cb):
    """Moving backbone atoms (N, C, O) in pred does not affect scRMSD.

    The superposition is performed on CA atoms only.  After alignment, the
    sidechain RMSD is computed solely from sidechain heavy atoms, so displacing
    N/C/O must leave the result unchanged.
    """
    ref  = make_chain_with_cb("AAS")
    pred = copy.deepcopy(ref)

    # Displace N, C, O in pred by a large offset; leave CA and CB in place.
    for res in pred.get_residues():
        for atom in res.get_atoms():
            if atom.get_name() in ("N", "C", "O"):
                atom.set_coord(atom.get_coord() + np.array([10.0, 10.0, 10.0]))

    result = per_residue_sidechain_rmsd("toy", "AAS", "AAS", ref, pred)

    # CB positions are identical → scRMSD should be 0.0 at every position.
    assert all(math.isclose(v, 0.0) for v in result.values())
