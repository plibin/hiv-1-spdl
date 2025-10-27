import numpy as np
from math import sqrt
from Bio.PDB import Residue, Chain
import core

#TODO: GPT told me TM should be computed on CA distances, can we verify this?
def _dist(r_ref: Residue.Residue, r_pred: Residue.Residue) -> float:
    sq_diff = core.squared_diffs_between_residues(r_ref, r_pred)
    #Euclidean norm
    return float(np.sqrt(np.mean(sq_diff)))

def global_tm(start: int, end: int,
              ref_chain: Chain.Chain,
              pred_chain: Chain.Chain,
              L_N: int) -> float:
    d_i_list = core.stat_per_residue(start, end, ref_chain, pred_chain, _dist)

    # As specified in the original TM-score paper: 10.1002/prot.20264,
    # d_i  is the distance between the i-th pair of aligned residues,
    # d_i_list is the list of distances over all considered residues,
    # d0 is a scale to normalize the match difference,
    # L_N is the length of the native structure,
    # L_T is the length of the aligned residues to the template structure,
    # (here that is the length of d_i_list).
    d_0 = 1.24 * (L_N - 15)**(1.0/3.0) - 1.8
    if d_0 < 0.5:
        d_0 = 0.5

    #TODO: GPT said: If L_N is small (<15), (L_N - 15)**(1.0/3.0) becomes a fractional power of a negative number, which Python treats as complex. That will break later math, but not in a super obvious way until runtime. I would propose to raise if this is the case (this should never happen...)
    #TODO: and some more: If you ask for --stat tm on a short segment (like 10 residues), you’ll hit the complex-number issue in the d_0 calculation before your own if d_0 < 0.5 clamp can “rescue” it. That’s an actual runtime failure.

    terms = [1.0 / (1.0 + (d_i / d_0)**2) for d_i in d_i_list]
    
    return sum(terms)/float(L_N)
