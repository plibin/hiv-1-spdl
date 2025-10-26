from Bio.PDB import Chain, Residue, Atom, Superimposer
import pandas as pd
import core
import numpy as np

def _squared_diffs_between_residues(ref_r: Residue.Residue, pred_r: Residue.Residue) -> float:
    names = set(a.get_name() for a in ref_r) & set(a.get_name() for a in pred_r)
    if not names:
        return np.nan
    squared_diffs = []
    for n in names:
        v = ref_r[n].get_coord() - pred_r[n].get_coord()
        squared_diffs.append(np.dot(v, v))
    return squared_diffs

def _rmsd_between_residues(r_ref: Residue.Residue, r_pred: Residue.Residue):
    squared_diffs = _squared_diffs_between_residues(r_ref, r_pred)
    return float(np.sqrt(np.mean(squared_diffs)))

#Return a dict with the positions and corresponding RMSDs
def per_residue_rmsd(start: int, end: int, ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(start, end, ref_chain, pred_chain, _rmsd_between_residues)

def global_rmsd(start: int, end: int, ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> float:
    squared_diffs = core.stat_per_residue(start, end, ref_chain, pred_chain, _squared_diffs_between_residues)

    return float(np.sqrt(np.mean(squared_diffs)))
