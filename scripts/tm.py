import numpy as np
from Bio.PDB import Residue, Chain

import core


# TM-score is defined on CA atoms (https://doi.org/10.1093/nar/gki524 Methods Section)
def _dist(r_ref: Residue.Residue, r_pred: Residue.Residue) -> float:
    if "CA" not in r_ref or "CA" not in r_pred:
        raise ValueError(f"Missing CA atom in residues {r_ref} or {r_pred}")
    diff = r_ref["CA"].get_coord() - r_pred["CA"].get_coord()
    return float(np.linalg.norm(diff))


def global_tm(id_: str,
              ref_align, pred_align,
              ref_chain: Chain.Chain,
              pred_chain: Chain.Chain) -> float:
    pos_to_d_i = core.stat_per_residue(id_, ref_align, pred_align, ref_chain, pred_chain, _dist)
    L_N = len(pos_to_d_i)

    # As specified in the original TM-score paper: 10.1002/prot.20264,
    # d_i  is the distance between the i-th pair of aligned residues,
    # d_i_list is the list of distances over all considered residues,
    # d0 is a scale to normalize the match difference,
    # L_N is the length of the native structure,
    # L_T is the length of the aligned residues to the template structure,
    # (here that is the length of d_i_list).

    # Handle small proteins to avoid complex numbers
    if L_N < 15:
        d_0 = 0.5
    else:
        d_0 = 1.24 * (L_N - 15) ** (1.0 / 3.0) - 1.8
        if d_0 < 0.5:
            d_0 = 0.5

    terms = [1.0 / (1.0 + (d_i / d_0) ** 2) for d_i in pos_to_d_i.values()]

    return sum(terms) / float(L_N)
