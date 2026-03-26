from Bio.PDB import Residue, Chain
from Bio.SeqRecord import SeqRecord
import core
import numpy as np


def _rmsd_between_residues(r_ref: Residue.Residue, r_pred: Residue.Residue):
    squared_diffs = core.squared_diffs_between_residues(r_ref, r_pred)
    return float(np.sqrt(np.mean(squared_diffs)))


# Return a dict with the positions and corresponding RMSDs
def per_residue_rmsd(id_: str,
                     ref_align, pred_align,
                     ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(id_, ref_align, pred_align, ref_chain, pred_chain, _rmsd_between_residues)


# Note: This computes the RMSD of all common atoms after superposition on CA atoms.
# This metric evaluates how well the side-chains align given the backbone alignment.
# The global RMSD computed here is the root-mean-square of the per-residue mean squared deviations,
# which weights each residue equally (regardless of atom count). When using all atoms, residues with
# more atoms contribute more to the total RMSD, which can bias results!!
# (http://pongor.itk.ppke.hu/library/Group-Publications/papers/142.pdf ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4321859/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC1471868/)
def global_rmsd(id_: str,
                ref_align, pred_align,
                ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> float:
    squared_diffs_per_pos = core.stat_per_residue(id_, ref_align, pred_align, ref_chain, pred_chain,
                                                  core.squared_diffs_between_residues)
    global_diffs = []
    for pos in squared_diffs_per_pos.keys():
        global_diffs.extend(squared_diffs_per_pos[pos])
    return float(np.sqrt(np.mean(global_diffs)))
