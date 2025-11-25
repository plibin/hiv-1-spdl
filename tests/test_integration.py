import os

import pytest
from Bio import PDB

from scripts.tm import global_tm
from scripts.rmsd import global_rmsd
from scripts.plddt import per_residue_plddt

# Path to structures used in tests
data_path = os.path.join(os.path.dirname(__file__), "data")
pdb_pred = os.path.join(data_path, "7leg_af2_pred.pdb")
pdb_ref = os.path.join(data_path, "7leg_ref.pdb")


def load_ca_chain(pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_path)
    # Use first chain
    chain = next(structure.get_chains())
    return chain


@pytest.fixture
def ref_chain():
    return load_ca_chain(pdb_ref)


@pytest.fixture
def pred_chain():
    return load_ca_chain(pdb_pred)


@pytest.fixture
def length(ref_chain, pred_chain):
    residues_ref = [r for r in ref_chain if r.has_id('CA')]
    residues_pred = [r for r in pred_chain if r.has_id('CA')]
    return min(len(residues_ref), len(residues_pred))


def test_tm_integration(length, ref_chain, pred_chain):
    tm_score = global_tm(0, length, ref_chain, pred_chain, length)
    assert 0 <= tm_score <= 1, "TM-score must be between 0 and 1"


def test_rmsd_integration(length, ref_chain, pred_chain):
    rmsd_value = global_rmsd(0, length, ref_chain, pred_chain)
    assert rmsd_value >= 0, "RMSD must be non-negative"


def test_plddt_integration(length, ref_chain, pred_chain):
    plddt_scores = per_residue_plddt(0, length, ref_chain, pred_chain)
    assert len(plddt_scores) == length, "pLDDT scores length mismatch"
    for score in plddt_scores.values():
        assert 0 <= score <= 100, "pLDDT score must be between 0 and 100"


if __name__ == "__main__":
    test_tm_integration(length, ref_chain, pred_chain)
    test_rmsd_integration(length, ref_chain, pred_chain)
    test_plddt_integration(length, ref_chain, pred_chain)
