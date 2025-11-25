import os

from Bio import PDB

from scripts.tm import global_tm

# Path to structures used in tests
data_path = os.path.join(os.path.dirname(__file__), "data")
pdb_pred = os.path.join(data_path, "7leg_af2_pred.pdb")
pdb_ref = os.path.join(data_path, "7leh_ref.pdb")


def load_ca_chain(pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_path)
    # Use first chain
    chain = next(structure.get_chains())
    return chain


def test_tm_integration():
    ref_chain = load_ca_chain(pdb_ref)
    pred_chain = load_ca_chain(pdb_pred)
    residues_ref = [r for r in ref_chain if r.has_id('CA')]
    residues_pred = [r for r in pred_chain if r.has_id('CA')]
    length = min(len(residues_ref), len(residues_pred))
    score = global_tm(0, length, ref_chain, pred_chain, length)
    assert 0 <= score <= 1, "TM-score must be between 0 and 1"


if __name__ == "__main__":
    test_tm_integration()
