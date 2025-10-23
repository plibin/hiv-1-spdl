def _rmsd_between_residues(ref_r: Residue.Residue, pred_r: Residue.Residue) -> float:
    names = set(a.get_name() for a in ref_r) & set(a.get_name() for a in pred_r)
    if not names:
        return np.nan
    diffs_sq = []
    for n in names:
        v = ref_r[n].get_coord() - pred_r[n].get_coord()
        diffs_sq.append(np.dot(v, v))
    return float(np.sqrt(np.mean(diffs_sq)))


def _align_pairs(ref: Chain.Chain, pred: Chain.Chain) -> List[Tuple[int, int]]:
    """Return aligned *indices into residue lists* (ref_i, pred_j) excluding gaps.

    Uses a simple global alignment (match=2, mismatch=-1, gap_open=-5, gap_ext=-0.5).
    """
    A = _seq_from_chain(ref)
    B = _seq_from_chain(pred)
    aln = pairwise2.align.globalms(A, B, 2, -1, -5, -0.5, one_alignment_only=True)[0]
    pairs: List[Tuple[int, int]] = []
    i = j = 0
    for a, b in zip(aln.seqA, aln.seqB):
        if a != "-" and b != "-":
            pairs.append((i, j))
        if a != "-":
            i += 1
        if b != "-":
            j += 1
    return pairs


def _superpose_on_ca(ref_res: Sequence[Residue.Residue], pred_res: Sequence[Residue.Residue], pairs: List[Tuple[int, int]]):
    """Compute superposition on matched CA atoms and APPLY to all atoms in pred chain.

    NOTE: This *modifies coordinates* of the provided `pred_res` atoms in-place.
    """
    ref_cas = []
    pred_cas = []
    for ri, pj in pairs:
        if "CA" in ref_res[ri] and "CA" in pred_res[pj]:
            ref_cas.append(ref_res[ri]["CA"]) 
            pred_cas.append(pred_res[pj]["CA"]) 
    if not ref_cas or not pred_cas:
        return
    sup = Superimposer()
    sup.set_atoms(ref_cas, pred_cas)
    # apply transform to every atom in pred residues
    all_pred_atoms: List[Atom.Atom] = [a for r in pred_res for a in r.get_atoms()]
    sup.apply(all_pred_atoms)


# --- public API --------------------------------------------------------------

def per_residue_rmsd(ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> pd.DataFrame:
    """Compute per-residue RMSD for *aligned* residue pairs between `ref_chain` and `pred_chain`.

    Returns a DataFrame with columns:
    - aln_idx: 1..N along the alignment (no gaps)
    - ref_resseq: residue number in the reference chain (resseq)
    - pred_resseq: residue number in the predicted chain (resseq)
    - ref_aa / pred_aa: one-letter codes
    - RMSD: Å, using shared atom names per residue

    Notes:
    - The function aligns sequences, superposes on matched CA atoms, then computes per-residue RMSD.
    - Coordinates of `pred_chain` are modified in-place by the superposition.
    """
    ref_res = list(ref_chain.get_residues())
    pred_res = list(pred_chain.get_residues())
    pairs = _align_pairs(ref_chain, pred_chain)
    _superpose_on_ca(ref_res, pred_res, pairs)

    rows = []
    for k, (ri, pj) in enumerate(pairs, start=1):
        r_ref = ref_res[ri]
        r_pred = pred_res[pj]
        het_r, resseq_r, ic_r = _resid_tuple(r_ref)
        het_p, resseq_p, ic_p = _resid_tuple(r_pred)
        rows.append(
            {
                "aln_idx": k,
                "ref_resseq": resseq_r,
                "pred_resseq": resseq_p,
                "ref_aa": _res_name_1(r_ref),
                "pred_aa": _res_name_1(r_pred),
                "RMSD": _rmsd_between_residues(r_ref, r_pred),
            }
        )
    return pd.DataFrame(rows)

def global_rmsd(ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> float:
    """Global RMSD after CA superposition on aligned pairs (Å)."""
    ref_res = list(ref_chain.get_residues())
    pred_res = list(pred_chain.get_residues())
    pairs = _align_pairs(ref_chain, pred_chain)
    _superpose_on_ca(ref_res, pred_res, pairs)

    # Compute RMSD on the matched CA atoms used for superposition
    ref_cas, pred_cas = [], []
    for ri, pj in pairs:
        if "CA" in ref_res[ri] and "CA" in pred_res[pj]:
            ref_cas.append(ref_res[ri]["CA"])
            pred_cas.append(pred_res[pj]["CA"])
    if not ref_cas:
        return float("nan")
    diffs = [(a.get_coord() - b.get_coord()) for a, b in zip(ref_cas, pred_cas)]
    return float(np.sqrt(np.mean([np.dot(d, d) for d in diffs])))
