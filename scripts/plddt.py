from Bio.PDB import Chain, Residue, Atom, Superimposer
import core

def _plldt(ref_r: Residue.Residue, pred_r: Residue.Residue) -> float:
    #TODO: why the CA atom?
    return float(pred_r["CA"].get_bfactor())

def per_residue_plddt(start: int, end: int,
                      ref_chain: Chain.Chain,
                      pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(start, end, ref_chain, pred_chain, _plldt)
