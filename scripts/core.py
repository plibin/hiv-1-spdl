from __future__ import annotations

import numpy as np
import copy
from typing import List, Sequence, Callable
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import Chain, Residue, Atom, Superimposer
from Bio.PDB.Polypeptide import is_aa as bio_is_aa

# Mapping for non-standard amino acids
AA_OVERRIDES = {
    "MSE": "MET",  # Selenomethionine -> MET
    "SEC": "CYS",  # Selenocysteine -> CYS (approximate)
    "PYL": "LYS",  # Pyrrolysine -> LYS (approximate)
    "CSO": "CYS",  # TODO
}


def squared_diffs_between_residues(ref_r: Residue.Residue, pred_r: Residue.Residue) -> list[float]:
    # Since we assume alignment, the overlap should be perfect,
    # but to be robust for missing artefacts in the ground truth PDB,
    # we follow this approach.
    names = set(a.get_name() for a in ref_r) & set(a.get_name() for a in pred_r)
    if not names:
        raise ValueError(f"No common atoms found between residues {ref_r} and {pred_r}")

    squared_diffs = []
    for n in names:
        v = ref_r[n].get_coord() - pred_r[n].get_coord()
        squared_diffs.append(np.dot(v, v))
    return squared_diffs


def _res_aa_letter(r: Residue.Residue) -> str:
    name = r.get_resname().strip().upper()
    name = AA_OVERRIDES.get(name, name).capitalize()
    one = protein_letters_3to1.get(name)
    if one is None: raise KeyError(f"Unknown residue name: {name}")
    return one


#TODO: no need for this mapping, use the Bio's is_aa directly
def is_aa(r: Residue.Residue) -> bool:
    return bio_is_aa(r)


def aa_seq(residues: List[Residue.Residue]) -> str:
    for r in residues:
        if not is_aa(r):
            raise RuntimeError("Not an amino acid residue!")
    return "".join(_res_aa_letter(aa) for aa in residues)


def aa_residues(chain: Chain.Chain) -> List[Residue.Residue]:
    return [r for r in chain.get_residues() if is_aa(r)]


def _ca_atoms(residues: Sequence[Residue.Residue]) -> List[Atom.Atom]:
    # r["CA"]: Retrieve the CA atom from residue r
    return [r["CA"] for r in residues if "CA" in r]


def stat_per_residue(start: int, end: int,
                     ref_chain: Chain.Chain, pred_chain: Chain.Chain,
                     stat: Callable[[Residue.Residue, Residue.Residue], float]) -> dict[int, float]:

    # Deep copy the predicted chain to avoid in-place mutation side effects
    pred_chain = copy.deepcopy(pred_chain)

    #Steps:
    #1. Filter to amino-acid residues only (skip waters/ligands/etc.).
    #2. We assume the amino acid sequences are aligned in the considered range (start-end),
    #     (we try to predict structures to compare them with their ground truth),
    #     and check this (error is raised if this is not met).
    #3. Superpose pred onto ref using CA pairs (in-place), *only* between start-end.
    #   Note: This means atoms outside the range are not necessarily aligned, which is intended for local analysis.
    #4. Apply a function stat() to each pair of residues of chain ref and pred, between start-end.
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

    # Assuming aligned sequences, we can just pair by index. Not trivial.
    if ref_seq == pred_seq:
        if len(ref_res) == 0:
            raise RuntimeError("Empty residue lists!")
        pairs = [(i, i) for i in range(start, end)]
    else:
        raise RuntimeError("Ref and pred seq are different between start-end!")
    # Check for out of bounds to avoid silent truncation
    if end > len(ref_res) or end > len(pred_res):
        raise IndexError(f"Range [{start}:{end}] extends past the end of chains (lengths: ref={len(ref_res)}, pred={len(pred_res)})")

    #3. Superpose pred onto ref using CA pairs (in-place), *only* between start-end.
    #  Note: This means atoms outside the range are not necessarily aligned, which is intended for local analysis.
    #  Note: Other methods might exist that use all backbone heavy atoms.
    # https://blog.matteoferla.com/2021/07/per-residue-rmsd.html?m=0 --> local superpositioning and calculation based on CA
    # https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut2.html ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4321859/ --> calculating based on CA
    ref_cas = []
    pred_cas = []
    for ri, pj in pairs:
        if "CA" in ref_res[ri] and "CA" in pred_res[pj]:
            ref_cas.append(ref_res[ri]["CA"])
            pred_cas.append(pred_res[pj]["CA"])
        else:
            raise RuntimeError("No CA atom in amino acid at pair (" + str(ri) + "," + str(pj) + ")")

    sup = Superimposer()
    sup.set_atoms(ref_cas, pred_cas)
    # Apply transform to every atom in pred residues to ensure the whole chain moves
    all_pred_atoms = [a for r in pred_res for a in r.get_atoms()]
    sup.apply(all_pred_atoms)

    #4. Apply a function stat to each pair of residues of chain ref and pred, between start-end.
    stats = {}
    for pos in range(start, end):
        r_ref = ref_res[pos]
        r_pred = pred_res[pos]
        stats[pos] = stat(r_ref, r_pred)

    #5. Return a dict of position with their corresponding score as computed by stat().
    return stats
