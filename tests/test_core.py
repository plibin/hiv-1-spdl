import copy

import numpy as np
import pytest
from Bio.PDB import Atom, Residue

from scripts.core import _find_alignment_start, aa_residues, aa_seq, squared_diffs_between_residues, stat_per_residue


def test_aa_seq(sample_chain):
    assert aa_seq(list(sample_chain)) == "AGS"


def test_aa_seq_supports_override_residue_names():
    mse = Residue.Residue(("H_MSE", 1, " "), "MSE", " ")
    mse.add(Atom.Atom("N", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 1, "N"))
    mse.add(Atom.Atom("CA", np.array([1.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 2, "C"))
    mse.add(Atom.Atom("C", np.array([2.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 3, "C"))
    assert aa_seq([mse]) == "M"


def test_aa_residues(sample_chain):
    residues = aa_residues(sample_chain)
    assert len(residues) == 3
    assert [residue.get_resname() for residue in residues] == ["ALA", "GLY", "SER"]


def test_find_alignment_start_identifies_first_shared_non_gap():
    assert _find_alignment_start("--AGSTV", "QQAGSTV") == (2, 7)


def test_find_alignment_start_rejects_length_mismatch():
    with pytest.raises(Exception, match="Alignment sizes differ"):
        _find_alignment_start("AGS", "AGST")


def test_find_alignment_start_rejects_no_shared_positions():
    with pytest.raises(Exception, match="Could not find alignment start point"):
        _find_alignment_start("A--", "-A-")


def test_squared_diffs_no_common_atoms():
    r1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    r1.add(Atom.Atom("N", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 1, "N"))

    r2 = Residue.Residue((" ", 1, " "), "ALA", " ")
    r2.add(Atom.Atom("C", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 2, "C"))

    with pytest.raises(ValueError, match="No common atoms"):
        squared_diffs_between_residues(r1, r2)


def test_stat_per_residue_does_not_mutate_input_chain(sample_chain):
    def dummy_stat(_ref, _pred):
        return 0.0

    pred_chain_shifted = copy.deepcopy(sample_chain)
    original_coords = [copy.deepcopy(atom.get_coord()) for atom in pred_chain_shifted.get_atoms()]
    for atom in pred_chain_shifted.get_atoms():
        atom.set_coord(atom.get_coord() + np.array([1.0, 1.0, 1.0]))

    shifted_coords = [copy.deepcopy(atom.get_coord()) for atom in pred_chain_shifted.get_atoms()]
    stats = stat_per_residue("toy", "AGS", "AGS", sample_chain, pred_chain_shifted, dummy_stat)

    assert stats == {0: 0.0, 1: 0.0, 2: 0.0}
    for atom, expected in zip(pred_chain_shifted.get_atoms(), shifted_coords):
        assert np.allclose(atom.get_coord(), expected)
    assert any(not np.allclose(shifted, original) for shifted, original in zip(shifted_coords, original_coords))


def test_stat_per_residue_rejects_gap_inside_start_motif(make_chain):
    chain = make_chain("AGSTV")

    with pytest.raises(Exception, match=r"Motif AG-ST has a gap"):
        stat_per_residue("toy", "AG-STV", "AG-STV", chain, copy.deepcopy(chain), lambda *_: 0.0)


def test_stat_per_residue_rejects_ambiguous_motif(make_chain):
    chain = make_chain("AGSTVAGSTV")

    with pytest.raises(Exception, match="Ambiguous motif"):
        stat_per_residue("toy", "AGSTVAGSTV", "AGSTVAGSTV", chain, copy.deepcopy(chain), lambda *_: 0.0)


