import pytest
from Bio.PDB import Residue, Atom
import numpy as np
from scripts.tm import global_tm, _dist


def test_dist_ca_only(sample_chain):
    """Test that _dist uses CA atoms."""
    r1 = list(sample_chain)[0]  # ALA
    r2 = list(sample_chain)[1]  # GLY

    # r1 CA is at [1, 0, 0]
    # r2 CA is at [4, 0, 0]
    # Distance should be 3.0

    d = _dist(r1, r2)
    assert d == 3.0


def test_dist_missing_ca():
    """Test that _dist raises ValueError if CA is missing."""
    r1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    r1.add(Atom.Atom("N", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 1, "N"))

    r2 = Residue.Residue((" ", 2, " "), "GLY", " ")
    r2.add(Atom.Atom("CA", np.array([1.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 2, "C"))

    with pytest.raises(ValueError, match="Missing CA"):
        _dist(r1, r2)


def test_global_tm_small_ln(sample_chain):
    tm = global_tm("toy", "AGS", "AGS", sample_chain, sample_chain)
    assert np.isclose(tm, 1.0)


def test_global_tm_returns_bounded_score(sample_chain):
    tm = global_tm("toy", "AGS", "AGS", sample_chain, sample_chain)
    assert 0.0 <= tm <= 1.0


def test_global_tm_rejects_empty_reference_alignment(sample_chain):
    with pytest.raises(ValueError, match="Reference alignment contains no residues"):
        global_tm("toy", "---", "---", sample_chain, sample_chain)
