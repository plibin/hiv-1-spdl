from Bio.PDB import Chain, Residue

import core


def _plldt(ref_r: Residue.Residue, pred_r: Residue.Residue) -> float:
    # pDDT is typically stored in the B-factor field.
    # AlphaFold stores the same value for all atoms in a residue, so checking CA is sufficient.
    # Code also should work for ESMFold predictions, only there the pLDDT values might be in eather [0,100] or [0,1] range.
    # (https://biomodel.uah.es/Jmol/fold_AI/index.html, https://310.ai/docs/folding/esmfold)
    return float(pred_r["CA"].get_bfactor())


def per_residue_plddt(start: int, end: int,
                      ref_chain: Chain.Chain,
                      pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(start, end, ref_chain, pred_chain, _plldt)
