import os
import requests
from Bio import PDB
from scripts.tm import global_tm


def download_file(url, fn):
    if not os.path.exists(fn):
        print(f"Downloading {fn}...")
        resp = requests.get(url)
        resp.raise_for_status()
        with open(fn, "wb") as f:
            f.write(resp.content)
    else:
        print(f"File {fn} already exists.")


# Human Serum Albumin: experimental and AlphaFold-predicted structures
pdb_gt_url = "https://files.rcsb.org/download/1AO6.pdb"
pdb_pred_url = "https://alphafold.ebi.ac.uk/files/AF-P02768-F1-model_v4.pdb"
pdb_true = "1AO6.pdb"
pdb_pred = "AF-P02768-2-F1-model_v6.pdb"


def load_ca_chain(pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_path)
    # Use first chain
    chain = next(structure.get_chains())
    return chain


def test_tm_human_albumin():
    download_file(pdb_gt_url, pdb_true)
    download_file(pdb_pred_url, pdb_pred)
    ref_chain = load_ca_chain(pdb_true)
    pred_chain = load_ca_chain(pdb_pred)
    residues_ref = [r for r in ref_chain if r.has_id('CA')]
    residues_pred = [r for r in pred_chain if r.has_id('CA')]
    length = min(len(residues_ref), len(residues_pred))
    print(f"Comparing {length} residues.")
    score = global_tm(0, length, ref_chain, pred_chain, length)
    print(f"TM-score: Human serum albumin X-ray (1AO6) vs AlphaFold2: {score:.4f}")
    assert 0 <= score <= 1, "TM-score must be between 0 and 1"


if __name__ == "__main__":
    test_tm_human_albumin()
