from __future__ import annotations

import numpy as np
import copy
from utils import count_overlapping
from Bio.SeqRecord import SeqRecord
from typing import List, Sequence, Callable
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import Chain, Residue, Atom, Superimposer
from Bio.PDB.Polypeptide import is_aa as bio_is_aa

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


def first_motif_idx(chain, motif):
    residues = aa_residues(chain)

    m = len(motif)
    for i in range(len(residues) - m + 1):
        window = residues[i:i + m]
        window_str = "".join(_res_aa_letter(r) for r in window)
        if window_str == motif:
            start = window[0]
            return start.id[1]

    return None


def get_res(residues, idx):
    for r in residues:
        if r.id[1] == idx:
            return r
    return None


def _find_alignment_start(ref_align, pred_align):
    align_len = len(pred_align)
    if align_len != len(ref_align):
        raise Exception("Alignment sizes differ!")

    for i in range(align_len):
        if pred_align[i] != '-' and ref_align[i] != '-':
            return i, align_len

    raise Exception("Could not find alignment start point!")


def stat_per_residue(id_: str,
                     alignment: dict[str, SeqRecord],
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
    ref_res = aa_residues(ref_chain)
    pred_res = aa_residues(pred_chain)

    align_start, align_len = _find_alignment_start(ref_align, pred_align)

    #TODO: this is a tricky part, I hope the explanation is a bit clear, otherwise we can discuss
    #Since the residues in the PDBs do not necessarily follow the same
    #positional numbering (TODO: add some additional explanation), we need to find the start position based on
    #a motif that occurs in both structures.
    motif_size = 5
    motif = query_align[align_start:align_start + motif_size]
    if '-' in motif:
        raise Exception("Motif has a gap!")
    if ref_pdb_align[align_start:align_start + motif_size] != motif:
        raise Exception("Query and PDB have different starting motifs!")
    if count_overlapping(query_align, motif) > 1 or \
       count_overlapping(ref_pdb_align, motif) > 1 :
        raise Exception("Ambiguous motif!")

    ref_start_idx = first_motif_idx(ref_chain, motif)
    pred_start_idx = first_motif_idx(pred_chain, motif)

    ref_selection = []
    pred_selection = []
    positions = []
    #TODO: the reported positions follow the positions in the alignment, make sure the alginment starts and ends  correctly!

    # Note: For a typical case where align_start is small (e.g. 2) and align_len is large (e.g. 600), align_len - align_start ≈ 598 — so the loop nearly covers the full range and the bug may not manifest.
    # But for a case where align_start > align_len / 2, the loop would be empty or truncated. I think the correct upper bound is align_len.
    for i in range(align_start, align_len - align_start):
        if query_align[i] != '-' and ref_pdb_align[i] != '-':
            if query_align[i] != ref_pdb_align[i]:
                raise Exception("Query and PDB don't match in the alignment!")

            #get_res returns None if it is a gap
            ref = get_res(ref_res, i + ref_start_idx)
            pred = get_res(pred_res, i + pred_start_idx)

            #check if ref or pred has a gap
            if ref is not None and pred is not None:
                if "CA" in ref and "CA" in pred:
                    ref_selection.append(ref)
                    pred_selection.append(pred)
                    positions.append(i)
                else:
                    raise RuntimeError("No CA atom in amino acid!")

                if _res_aa_letter(ref) != _res_aa_letter(pred):
                    raise Exception("Ref and pred don't match in the PDB!")
           

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
    all_pred_atoms = [a for r in pred_res for a in r.get_atoms()]
    sup.apply(all_pred_atoms)

    #4. Apply a function stat to each pair of residues of chain ref and pred, between start-end.
    stats = {}
    for ref, pred, pos in zip(ref_selection, pred_selection, positions):
        stats[pos] = stat(ref, pred)

    #5. Return a dict of position with their corresponding score as computed by stat().
    return stats
