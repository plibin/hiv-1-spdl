from __future__ import annotations
from typing import Iterable, List, Sequence, Tuple
import numpy as np
import pandas as pd
from Bio.PDB import Chain, Residue, Atom, Superimposer
from Bio import pairwise2
from Bio.Data.IUPACData import protein_letters_3to1

# --- helpers -----------------------------------------------------------------

AA_OVERRIDES = {
    "MSE": "MET",  # Selenomethionine → MET
    "SEC": "CYS",  # Selenocysteine → CYS (approximate)
    "PYL": "LYS",  # Pyrrolysine → LYS (approximate)
}


def _res_name_1(r: Residue.Residue) -> str:
    name3 = (r.get_resname() or "").upper()
    name3 = AA_OVERRIDES.get(name3, name3).capitalize()
    return protein_letters_3to1.get(name3, "X")


def _seq_from_chain(chain: Chain.Chain) -> str:
    return "".join(_res_name_1(r) for r in chain.get_residues())


def _resid_tuple(r: Residue.Residue) -> Tuple[str, int, str]:
    het, resseq, icode = r.get_id()
    return het, int(resseq), icode or ""


def _ca_atoms(residues: Sequence[Residue.Residue]) -> List[Atom.Atom]:
    return [r["CA"] for r in residues if "CA" in r]
