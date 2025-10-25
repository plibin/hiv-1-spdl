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

def _rmsd_between_residues(r_ref, r_pred):
    squared_diffs = _squared_diffs_between_residues(r_ref, r_pred)
    return float(np.sqrt(np.mean(squared_diffs)))

def per_residue_rmsd(ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> pd.DataFrame:
    """Compute per-residue RMSD for *aligned* residue pairs between `ref_chain` and `pred_chain`.

    Returns a DataFrame with columns pos and RMSD.
    """
    core._superpose_on_ca(ref_chain, pred_chain)

    #we cannot directly use the residues,
    #as there might be other atoms (especially in the ground truth) as well,
    #and we are only interested in positional differences (between the AAs)
    ref_res = core._aa_residues(ref_chain)
    pred_res = core._aa_residues(pred_chain)
    
    rows = []
    for pos in range(len(ref_res)):
        r_ref = ref_res[pos]
        r_pred = pred_res[pos]

        rows.append(
            {
                "pos": pos + 1,
                "RMSD": _rmsd_between_residues(r_ref, r_pred)
            }
        )
    return pd.DataFrame(rows)

def global_rmsd(ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> float:
    core._superpose_on_ca(ref_chain, pred_chain)

    ref_res = list(ref_chain.get_residues())
    pred_res = list(pred_chain.get_residues())

    squared_diffs = []
    for pos in range(1, len(ref_res) + 1):
        r_ref = ref_res[pos]
        r_pred = pred_res[pos]
        squared_diffs.append(_squared_diffs_between_residues(r_ref, r_pred))
        
    return float(np.sqrt(np.mean(squared_diffs)))
