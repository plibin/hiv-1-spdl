from Bio.PDB import Residue, Chain
import core
import numpy as np

def _rmsd_between_residues(r_ref: Residue.Residue, r_pred: Residue.Residue):
    squared_diffs = core.squared_diffs_between_residues(r_ref, r_pred)
    return float(np.sqrt(np.mean(squared_diffs)))

#Return a dict with the positions and corresponding RMSDs
def per_residue_rmsd(start: int, end: int, ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(start, end, ref_chain, pred_chain, _rmsd_between_residues)

#TODO: GPT says: we superpose on CA atoms, but compute distances on all atoms -> for TM, it should be on CA atoms by definition, but if so, should we not do it here as well? anyhow, we need to be *very* clear about the choices in the methods section 
def global_rmsd(start: int, end: int, ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> float:
    squared_diffs_per_pos = core.stat_per_residue(start, end, ref_chain, pred_chain, core.squared_diffs_between_residues)
    global_diffs = []
    for pos in squared_diffs_per_pos.keys():
        global_diffs.extend(squared_diffs_per_pos[pos])
    return float(np.sqrt(np.mean(global_diffs)))
