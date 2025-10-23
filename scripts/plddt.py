def per_residue_plddt(pred_chain: Chain.Chain) -> pd.DataFrame:
    """Extract per-residue pLDDT from CA B-factors of a *predicted* chain.

    Returns DataFrame with columns:
    - pred_resseq: residue number from the chain
    - pLDDT: float (taken from CA B-factor)
    """
    rows = []
    for r in pred_chain.get_residues():
        if "CA" not in r:
            continue
        _, resseq, _ = _resid_tuple(r)
        rows.append({"pred_resseq": resseq, "pLDDT": float(r["CA"].get_bfactor())})
    return pd.DataFrame(rows)
