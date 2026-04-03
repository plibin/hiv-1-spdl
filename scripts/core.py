from __future__ import annotations

import copy
import sys
from typing import List, Sequence, Callable

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import Chain, Residue, Atom, Superimposer
from Bio.PDB.Polypeptide import is_aa as bio_is_aa

from utils import count_overlapping

# Mapping for non-standard amino acids
AA_OVERRIDES = {
    "MSE": "MET",  # Selenomethionine -> MET
    "SEC": "CYS",  # Selenocysteine -> CYS (approximate)
    "PYL": "LYS",  # Pyrrolysine -> LYS (approximate)
    # Note on CSO and CAF: CSO is S-hydroxycysteine and CAF is a cross-linked cysteine adduct.
    # Treating both as plain CYS is defensible as a structural approximation but should be noted explicitly in the paper.
    "CSO": "CYS",
    "CAF": "CYS",
    # Note on CSD: CSD is structurally distinctive, with long sulfur–carbon and sulfur–oxygen bonds and tetrahedral geometry around sulfur due to its lone pair,
    # and it shows a greater preference for α-helix than unmodified CYS. The modification alters local backbone preferences. For the current CA-based RMSD workflow, mapping CSD → CYS is correct.
    "CSD": "CYS",
    # Note on CAS: The backbone is identical to CYS, so for Cα-based RMSD, mapping CAS to CYS is structurally safe as a fallback.
    # However, CAS is rare, context-dependent, and might warrant runtime waring instead of silent override.
    "CAS": "CYS",
}


def squared_diffs_between_residues(ref_r: Residue.Residue, pred_r: Residue.Residue) -> list[float]:
    # Since we assume alignment, the overlap should be perfect,
    # but to be robust for missing artefacts in the ground truth PDB,
    # we follow this approach.
    names = set(a.get_name() for a in ref_r) & set(a.get_name() for a in pred_r)
    if not names:
        raise ValueError(f"No common atoms found between residues {ref_r} and {pred_r}")

    squared_diffs = []
    for n in names:
        v = ref_r[n].get_coord() - pred_r[n].get_coord()
        squared_diffs.append(np.dot(v, v))
    return squared_diffs


def _res_aa_letter(r: Residue.Residue) -> str:
    name = r.get_resname().strip().upper()
    name = AA_OVERRIDES.get(name, name).capitalize()
    one = protein_letters_3to1.get(name)
    if one is None: raise KeyError(f"Unknown residue name: {name}")
    return one


def aa_seq(residues: List[Residue.Residue]) -> str:
    for r in residues:
        if not bio_is_aa(r):
            raise RuntimeError("Not an amino acid residue!")
    return "".join(_res_aa_letter(aa) for aa in residues)


def aa_residues(chain: Chain.Chain) -> List[Residue.Residue]:
    return [r for r in chain.get_residues() if bio_is_aa(r)]


def _ca_atoms(residues: Sequence[Residue.Residue]) -> List[Atom.Atom]:
    # r["CA"]: Retrieve the CA atom from residue r
    return [r["CA"] for r in residues if "CA" in r]


def first_motif_res(id_, chain, motif):
    residues = aa_residues(chain)

    m = len(motif)
    for i in range(len(residues) - m + 1):
        window = residues[i:i + m]
        window_str = "".join(_res_aa_letter(r) for r in window)
        if window_str == motif:
            start = window[0]
            return start

    raise RuntimeError(f"{id_}: No start for motif {motif} in chain")


def next_res(chain, current_res):
    residues = aa_residues(chain)
    for i, r in enumerate(residues):
        if r is current_res:
            return residues[i + 1] if i + 1 < len(residues) else None
    raise RuntimeError("Current residue not found in residue list")


def _find_alignment_start(ref_align, pred_align):
    align_len = len(pred_align)
    if align_len != len(ref_align):
        raise Exception("Alignment sizes differ!")

    for i in range(align_len):
        if pred_align[i] != '-' and ref_align[i] != '-':
            return i, align_len

    raise Exception("Could not find alignment start point!")


def _resolve_start_motif(ref_align, pred_align, align_start: int, motif_size: int = 5) -> str:
    motif = pred_align[align_start:align_start + motif_size]
    if '-' in motif:
        raise Exception(f"Motif {motif} has a gap!")
    if ref_align[align_start:align_start + motif_size] != motif:
        raise Exception("Query and PDB have different starting motifs!")
    if count_overlapping(pred_align, motif) > 1 or count_overlapping(ref_align, motif) > 1:
        raise Exception("Ambiguous motif!")
    return motif


def stat_per_residue(id_: str,
                     ref_align, pred_align,
                     ref_chain: Chain.Chain, pred_chain: Chain.Chain,
                     stat: Callable[[Residue.Residue, Residue.Residue], float]) -> dict[int, float]:
    # Deep copy the predicted chain to avoid in-place mutation side effects
    pred_chain = copy.deepcopy(pred_chain)

    #Steps:
    #1. Filter to amino-acid residues only (skip waters/ligands/etc.).
    #2. We depend on a correct amino acid alignment of a/ the query sequences and b/ the sequences extracted from the PDBs. Based on this alignment, we superpose pred onto ref using CA pairs (in-place), *only* when there is an amino acid in both the pred and the ref, as their might be gaps.
    #3. We apply a function stat() to each pair of superposed residues of chain ref and pred.
    #4. We return a dict of position (as in the alignment!) with their corresponding score as computed by stat().

    #1. Filter to amino-acid residues only.
    ref_residues = aa_residues(ref_chain)
    pred_residues = aa_residues(pred_chain)

    align_start, align_len = _find_alignment_start(ref_align, pred_align)

    # Find and validate a unique motif at alignment start to anchor chain numbering.
    motif = _resolve_start_motif(ref_align, pred_align, align_start)

    # We use residues to navigate through the chains (jumping from one res to the next one),
    # since the residue index does not necessarily follow the gaps in the alignment. 
    ref_res = first_motif_res(id_, ref_chain, motif)
    pred_res = first_motif_res(id_, pred_chain, motif)

    ref_selection = []
    pred_selection = []
    positions = []

    #TODO!!!: given that this has been a particular hard to get right,
    #perhaps a sanity check is warranted, write the  
    for i in range(align_start, align_len):
        # Both the pred and ref have a gap in the alignment: safe to skip.
        if pred_align[i] == '-' and ref_align[i] == '-':
            continue

        # Unexpected gap in the pred (i.e., query): raise error
        if pred_align[i] == '-' and ref_align[i] != '-':
            raise RuntimeError(f"Unexpected gap in pred at pos {i} for id {id_}")

        # It can happen that the ref PDB has gaps, as it is the result of a wet lab experiment.
        # We ignore the position, but first move to the next residue for the pred.
        if pred_align[i] != '-' and ref_align[i] == '-':
            pred_res = next_res(pred_chain, pred_res)
            continue

        if pred_align[i] != ref_align[i]:
            print(f"{id_}: Query {pred_align[i]} and PDB {ref_align[i]} don't match in the alignment at pos {i}",
                  file=sys.stderr)

        if ref_res is None or pred_res is None:
            raise RuntimeError("Missing residues")

        ref_res_aa = _res_aa_letter(ref_res)
        pred_res_aa = _res_aa_letter(pred_res)
        if ref_res_aa != pred_res_aa:
            print(f"{id_}: ref_aa {ref_res_aa} and pred_aa {pred_res_aa} don't match in the PDB at pos {i}",
                  file=sys.stderr)

        if "CA" in ref_res and "CA" in pred_res:
            ref_selection.append(ref_res)
            pred_selection.append(pred_res)
            positions.append(i)

            pred_res = next_res(pred_chain, pred_res)
            ref_res = next_res(ref_chain, ref_res)
        else:
            raise RuntimeError("No CA atom in amino acid!")

    #3. Superpose pred onto ref using CA pairs (in-place), *only* in the region where the amino acids overlap.
    #  Note: This means atoms outside the range are not necessarily aligned, which is intended for local analysis.
    #  Note: Other methods might exist that use all backbone heavy atoms.
    # https://blog.matteoferla.com/2021/07/per-residue-rmsd.html?m=0 --> local superpositioning and calculation based on CA
    # https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut2.html ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4321859/ --> calculating based on CA
    sup = Superimposer()

    #Extract CA atoms
    ref_cas = [x["CA"] for x in ref_selection]
    pred_cas = [x["CA"] for x in pred_selection]

    sup.set_atoms(ref_cas, pred_cas)
    # Apply transform to every atom in pred residues to ensure the whole chain moves
    all_pred_atoms = [a for r in pred_residues for a in r.get_atoms()]
    sup.apply(all_pred_atoms)

    #4. Apply a function stat to each pair of residues of chain ref and pred, between start-end.
    stats = {}
    for ref, pred, pos in zip(ref_selection, pred_selection, positions):
        stats[pos] = stat(ref, pred)

    #5. Return a dict of position with their corresponding score as computed by stat().
    return stats
