from Bio.PDB import Chain, Residue, Atom, Superimposer
from . import core

def _plldt(ref_r: Residue.Residue, pred_r: Residue.Residue) -> float:
    # pLDDT is typically stored in the B-factor field. 
    # AlphaFold stores the same value for all atoms in a residue, so checking CA is sufficient.
    return float(pred_r["CA"].get_bfactor())

def per_residue_plddt(start: int, end: int,
                      ref_chain: Chain.Chain,
                      pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(start, end, ref_chain, pred_chain, _plldt)
