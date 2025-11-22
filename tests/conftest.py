import pytest
from Bio.PDB import Structure, Model, Chain, Residue, Atom
import numpy as np

@pytest.fixture
def sample_chain():
    """Creates a sample chain with 3 residues (ALA, GLY, SER)."""
    chain = Chain.Chain("A")
    
    # Residue 1: ALA
    r1 = Residue.Residue((" ", 1, " "), "ALA", " ")
    r1.add(Atom.Atom("N", np.array([0.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 1, "N"))
    r1.add(Atom.Atom("CA", np.array([1.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 2, "C"))
    r1.add(Atom.Atom("C", np.array([2.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 3, "C"))
    r1.add(Atom.Atom("O", np.array([2.0, 1.0, 0.0]), 1.0, 0.0, " ", "O", 4, "O"))
    chain.add(r1)

    # Residue 2: GLY
    r2 = Residue.Residue((" ", 2, " "), "GLY", " ")
    r2.add(Atom.Atom("N", np.array([3.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 5, "N"))
    r2.add(Atom.Atom("CA", np.array([4.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 6, "C"))
    r2.add(Atom.Atom("C", np.array([5.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 7, "C"))
    r2.add(Atom.Atom("O", np.array([5.0, 1.0, 0.0]), 1.0, 0.0, " ", "O", 8, "O"))
    chain.add(r2)

    # Residue 3: SER
    r3 = Residue.Residue((" ", 3, " "), "SER", " ")
    r3.add(Atom.Atom("N", np.array([6.0, 0.0, 0.0]), 1.0, 0.0, " ", "N", 9, "N"))
    r3.add(Atom.Atom("CA", np.array([7.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 10, "C"))
    r3.add(Atom.Atom("C", np.array([8.0, 0.0, 0.0]), 1.0, 0.0, " ", "C", 11, "C"))
    r3.add(Atom.Atom("O", np.array([8.0, 1.0, 0.0]), 1.0, 0.0, " ", "O", 12, "O"))
    chain.add(r3)
    
    return chain

@pytest.fixture
def sample_structure(sample_chain):
    """Creates a sample structure containing the sample chain."""
    model = Model.Model(0)
    model.add(sample_chain)
    structure = Structure.Structure("test")
    structure.add(model)
    return structure
