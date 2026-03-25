import pytest
from Bio.PDB import Chain, Residue, Atom
import numpy as np
from Bio.PDB.Polypeptide import is_aa as bio_is_aa
from scripts.core import aa_seq, aa_residues, stat_per_residue, squared_diffs_between_residues

def test_is_aa(sample_chain):
    """Test is_aa function."""
    for residue in sample_chain:
        assert bio_is_aa(residue)

    # Create a water residue
    water = Residue.Residue(("W", 1, " "), "HOH", " ")
    assert not bio_is_aa(water)

def test_aa_seq(sample_chain):
    """Test aa_seq function."""
    assert aa_seq(list(sample_chain)) == "AGS"

def test_aa_residues(sample_chain):
    """Test aa_residues function."""
    residues = aa_residues(sample_chain)
    assert len(residues) == 3
    assert residues[0].get_resname() == "ALA"
    assert residues[1].get_resname() == "GLY"
    assert residues[2].get_resname() == "SER"

def test_is_aa_overrides():
    """Test is_aa with overridden residues."""
    # MSE is mapped to MET, so it should be considered an AA if we treat it as such
    # But Bio.PDB might read it as HETATM. Let's simulate a HETATM MSE.
    mse = Residue.Residue(("H_MSE", 1, " "), "MSE", " ")
    assert bio_is_aa(mse)

def test_squared_diffs_no_common_atoms():
    """Test that ValueError is raised when no common atoms exist."""
    r1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    r1.add(Atom.Atom("N", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 1, "N"))
    
    r2 = Residue.Residue((" ", 1, " "), "ALA", " ")
    # No common atoms
    r2.add(Atom.Atom("C", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 2, "C"))
    
    with pytest.raises(ValueError, match="No common atoms"):
        squared_diffs_between_residues(r1, r2)

def test_stat_per_residue_mutation(sample_chain):
    """Test that stat_per_residue does not mutate the original chain."""
    # Create a dummy stat function
    def dummy_stat(r1, r2):
        return 0.0
    
    # Use the same chain for ref and pred
    ref_chain = sample_chain
    pred_chain = sample_chain # This is the one that shouldn't be mutated in place if we passed it directly
    
    # But wait, sample_chain is a fixture, so it's a fresh object.
    # Let's make a copy to compare later.
    import copy
    original_coords = [a.get_coord() for r in pred_chain for a in r]
    original_coords = copy.deepcopy(original_coords)
    
    # We need to pass a pred_chain that is actually different so superposition does something
    # But here we just want to check if the *object* passed is mutated.
    # If we pass the same object as ref and pred, and it aligns to itself, it won't move.
    # So let's create a shifted pred chain.
    pred_chain_shifted = copy.deepcopy(sample_chain)
    for atom in pred_chain_shifted.get_atoms():
        atom.set_coord(atom.get_coord() + np.array([1.0, 1.0, 1.0]))
        
    stat_per_residue(0, 3, ref_chain, pred_chain_shifted, dummy_stat)
    
    # Check if pred_chain_shifted coords are still the same (shifted)
    # Because stat_per_residue should work on a copy.
    for i, atom in enumerate(pred_chain_shifted.get_atoms()):
        assert np.allclose(atom.get_coord(), original_coords[i] + np.array([1.0, 1.0, 1.0]))

def test_stat_per_residue_bounds(sample_chain):
    """Test out of bounds check."""
    def dummy_stat(r1, r2):
        return 0.0
    
    with pytest.raises(IndexError, match="extends past the end"):
        stat_per_residue(0, 10, sample_chain, sample_chain, dummy_stat)

