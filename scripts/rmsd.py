from typing import List

import numpy as np
from Bio.PDB import Residue, Chain

import core


def _rmsd_between_residues(r_ref: Residue.Residue, r_pred: Residue.Residue):
    squared_diffs = core.squared_diffs_between_residues(r_ref, r_pred)
    return float(np.sqrt(np.mean(squared_diffs)))


def _sidechain_heavy_atom_rmsd_between_residues(
        r_ref: Residue.Residue,
        r_pred: Residue.Residue,
) -> float:
    squared_diffs = core.squared_diffs_between_sidechain_heavy_atoms(r_ref, r_pred)

    # Some residues have no comparable side-chain heavy atoms; represent these
    # positions as missing values in the positional side-chain RMSD output.
    if not squared_diffs:
        return float("nan")

    return float(np.sqrt(np.mean(squared_diffs)))


# Return a dict with the positions and corresponding RMSDs
def per_residue_rmsd(id_: str,
                     ref_align, pred_align,
                     ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(id_, ref_align, pred_align, ref_chain, pred_chain, _rmsd_between_residues)


def per_residue_sidechain_rmsd(id_: str,
                               ref_align, pred_align,
                               ref_chain: Chain.Chain, pred_chain: Chain.Chain) -> dict[int, float]:
    return core.stat_per_residue(id_, ref_align, pred_align, ref_chain, pred_chain,
                                 _sidechain_heavy_atom_rmsd_between_residues)


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


def _all_atom_rmsd_stat(
        ref_selection: List[Residue.Residue],
        pred_selection: List[Residue.Residue],
) -> float:
    """Compute all-heavy-atom RMSD over aligned residue selections.

    Uses the intersection of atom names per residue pair (same strategy as
    squared_diffs_between_residues) so that atoms unresolved in the reference
    or absent from the prediction are silently excluded.  Hydrogens in the
    prediction (e.g. AlphaFold2) are automatically excluded because they are
    absent from the experimental reference, so the intersection never contains
    them.

    Returns NaN when the prediction contains no sidechain heavy atoms,
    signalling that the all-atom metric is not applicable (backbone-only model).
    """
    has_sidechain = any(
        core._is_sidechain_heavy_atom(a)
        for r in pred_selection
        for a in r.get_atoms()
        if not core._is_hydrogen(a)
    )
    if not has_sidechain:
        print("Warning: prediction has no sidechain heavy atoms; all-atom RMSD not applicable.", file=__import__("sys").stderr)
        return float("nan")

    all_squared_diffs: list[float] = []
    for ref_r, pred_r in zip(ref_selection, pred_selection):
        all_squared_diffs.extend(core.squared_diffs_between_residues(ref_r, pred_r))

    if not all_squared_diffs:
        return float("nan")
    return float(np.sqrt(np.mean(all_squared_diffs)))


def global_all_atom_rmsd(
        id_: str,
        ref_align, pred_align,
        ref_chain: Chain.Chain,
        pred_chain: Chain.Chain,
) -> float:
    """All-heavy-atom RMSD after optimal superposition on all heavy atoms.

    Unlike global_rmsd — which uses Cα atoms for superposition then evaluates
    RMSD over all common atoms — this function derives the optimal rigid-body
    transform by minimising the RMSD over the full set of matched heavy atoms
    (backbone + sidechain), giving the best achievable all-atom alignment.

    Because the optimal all-atom superposition minimises the all-atom RMSD,
    this value is always ≤ global_rmsd for full-atom models.

    Returns NaN for backbone-only predictions (ESM3-Open, Ember3D) where the
    metric is not applicable.
    """
    return core.stat_global(
        id_, ref_align, pred_align, ref_chain, pred_chain,
        stat=_all_atom_rmsd_stat,
        superposition_atoms=core._heavy_atom_pairs,
    )
