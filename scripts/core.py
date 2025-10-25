from __future__ import annotations
from typing import Iterable, List, Sequence, Tuple
import numpy as np
import pandas as pd
from Bio.PDB import Chain, Residue, Atom, Superimposer
from Bio import pairwise2
from Bio.Data.IUPACData import protein_letters_3to1

AA_OVERRIDES = {
    "MSE": "MET",  # Selenomethionine → MET
    "SEC": "CYS",  # Selenocysteine → CYS (approximate)
    "PYL": "LYS",  # Pyrrolysine → LYS (approximate)
}

def _res_name(r: Residue.Residue) -> str:
    name = r.get_resname().strip().upper()
    name = AA_OVERRIDES.get(name, name).capitalize()
    return protein_letters_3to1.get(name)

def _seq_from_chain(chain: Chain.Chain) -> str:
    AAs = _aa_residues(chain)
        
    return "".join(_res_name(aa) for aa in AAs)

def _aa_residues(chain: Chain.Chain) -> List[Residue.Residue]:
    AAs = []
    for r in chain.get_residues():
        if r.id[0] != " ":  # Not a standard amino acid
            continue
        AAs.append(r)
    return AAs

#TODO: is this used?
def _resid_tuple(r: Residue.Residue) -> Tuple[str, int, str]:
    het, resseq, icode = r.get_id()
    return het, int(resseq), icode or ""

def _ca_atoms(residues: Sequence[Residue.Residue]) -> List[Atom.Atom]:
    #r["CA"]: Retrieve the CA atom from residue r
    return [r["CA"] for r in residues if "CA" in r]

def _superpose_on_ca(ref_chain: Chain.Chain, pred_chain: Chain.Chain):
    """Compute superposition on matched CA atoms and apply to all atoms in pred chain.
       Changes the coordinates of pred_chain in place.
    """
    ref_seq = _seq_from_chain(ref_chain)
    pred_seq = _seq_from_chain(pred_chain)
    print(ref_seq)
    print(pred_seq)
    if ref_seq == pred_seq:
        pairs = [(i, i) for i in range(len(ref_seq))]
    else:
        raise RuntimeError("Ref and pred seq are different!")

    ref_res = list(ref_chain.get_residues())
    pred_res = list(pred_chain.get_residues())
    
    ref_cas = []
    pred_cas = []
    for ri, pj in pairs:
        if "CA" in ref_res[ri] and "CA" in pred_res[pj]:
            ref_cas.append(ref_res[ri]["CA"]) 
            pred_cas.append(pred_res[pj]["CA"]) 

    if not ref_cas or not pred_cas:
        #TODO: throw exception? (raise a runtime error)
        return

    sup = Superimposer()
    sup.set_atoms(ref_cas, pred_cas)
    # apply transform to every atom in pred residues
    all_pred_atoms = [a for r in pred_res for a in r.get_atoms()]
    sup.apply(all_pred_atoms)
