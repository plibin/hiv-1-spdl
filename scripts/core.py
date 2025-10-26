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

def _res_aa_letter(r: Residue.Residue) -> str:
    name = r.get_resname().strip().upper()
    name = AA_OVERRIDES.get(name, name).capitalize()
    return protein_letters_3to1.get(name)

def is_aa(r: Residue.Residue) -> bool:
    return r.id[0] == " "

def aa_seq(residues: List[Residue.Residue]) -> str:
    for r in residues:
        if not is_aa(r):
            raise RuntimeError("Not a amino acide residue!")
    return "".join(_res_aa_letter(aa) for aa in residues)

def aa_residues(chain: Chain.Chain) -> List[Residue.Residue]:
    return [r for r in chain.get_residues() if is_aa(r)]

def _ca_atoms(residues: Sequence[Residue.Residue]) -> List[Atom.Atom]:
    #r["CA"]: Retrieve the CA atom from residue r
    return [r["CA"] for r in residues if "CA" in r]

#TODO: specify signature of state
def stat_per_residue(start: int, end: int,
                     ref_chain: Chain.Chain, pred_chain: Chain.Chain,
                     stat: Callable[Residue.Residue, Residue.Residue, float]) -> dict[int, float]:
    #Steps:
    #1. Filter to amino-acid residues only (skip waters/ligands/etc.).
    #2. We assume the amino acid sequences are aligned in the considered range (start-end),
    #     and check this (error is raised if this is not met)
    #     (we try to predict structures to compare them with their ground truth)
    #3. Superpose pred onto ref using CA pairs (in-place), *only* between start-end.
    #     TODO: is this correct? (that means, outside the considered range, the atoms are not aligned!)
    #4. Apply a function stat to each pair of residues of chain ref and pred, between start-end.
    #5. Return a dict of position with their corresponding score as computed by stat().
    # Since it is possible that some atoms are not aligned,
    # it is important to keep this workflow together to ensure consistency!
    
    #1. Filter to amino-acid residues only.
    ref_res = aa_residues(ref_chain)
    pred_res = aa_residues(pred_chain)

    #2. Check if the amino acid sequences are aligned (i.e., equal) in the range (start-end).
    ref_seq = aa_seq(ref_res)[start:end]
    pred_seq = aa_seq(pred_res)[start:end]
    pairs = None
    if ref_seq == pred_seq:
        if len(ref_res) == 0:
            raise RuntimeError("Empty residue lists!")
        pairs = [(i, i) for i in range(start, end)]
    else:
        raise RuntimeError("Ref and pred seq are different between start-end!")

    #3. Superpose pred onto ref using CA pairs (in-place), *only* between start-end.
    ref_cas = []
    pred_cas = []
    for ri, pj in pairs:
        if "CA" in ref_res[ri] and "CA" in pred_res[pj] :
            ref_cas.append(ref_res[ri]["CA"]) 
            pred_cas.append(pred_res[pj]["CA"])
        else :
            raise RuntimeError("No CA atom in amino acid at pair (" + str(ri) + "," + str(pj) + ")")

    sup = Superimposer()
    sup.set_atoms(ref_cas, pred_cas)
    # TODO?: apply transform to every atom in pred residues
    all_pred_atoms = [a for r in pred_res for a in r.get_atoms()]
    #TODO: not sure if this is necessary, or even allowed?
    sup.apply(all_pred_atoms)

    #4. Apply a function stat to each pair of residues of chain ref and pred, between start-end.
    stats = {}
    for pos in range(start, end):
        r_ref = ref_res[pos]
        r_pred = pred_res[pos]
        stats[pos] = stat(r_ref, r_pred)

    #5. Return a dict of position with their corresponding score as computed by stat().
    return stats
