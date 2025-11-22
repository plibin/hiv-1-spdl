import pytest
from Bio.PDB import Chain, Residue, Atom
import numpy as np
from scripts.tm import global_tm, _dist

def test_dist_ca_only(sample_chain):
    """Test that _dist uses CA atoms."""
    r1 = list(sample_chain)[0] # ALA
    r2 = list(sample_chain)[1] # GLY
    
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
    """Test global_tm with small L_N."""
    # L_N = 3, so d_0 should be 0.5
    # We pass the same chain, so distance is 0 (after superposition)
    # TM score should be 1.0 (sum(1/(1+0)) / L_N * L_aligned?)
    # Wait, global_tm returns sum(terms) / L_N
    # If perfect match, terms are all 1.0.
    # If we align 3 residues, sum is 3.0.
    # If L_N is 3, result is 1.0.
    
    tm = global_tm(0, 3, sample_chain, sample_chain, L_N=3)
    assert np.isclose(tm, 1.0)

def test_global_tm_complex_check(sample_chain):
    """Test that global_tm handles L_N < 15 without error."""
    # This would crash if we didn't handle the negative root
    tm = global_tm(0, 3, sample_chain, sample_chain, L_N=10)
    # d_0 = 0.5
    # Perfect match -> sum = 3.0
    # Result = 3.0 / 10 = 0.3
    assert np.isclose(tm, 0.3)
