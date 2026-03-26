import sys
from pathlib import Path
from typing import Optional

import pytest
from Bio.PDB import Chain, Residue, Atom
from Bio.Data.IUPACData import protein_letters_1to3
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _three_letter(aa: str) -> str:
    return protein_letters_1to3[aa].upper()


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def data_root(repo_root: Path) -> Path:
    return repo_root / "data"


@pytest.fixture
def make_chain():
    def _make_chain(
        sequence: str,
        *,
        chain_id: str = "A",
        translation: tuple[float, float, float] = (0.0, 0.0, 0.0),
        ca_bfactors: Optional[list[float]] = None,
    ) -> Chain.Chain:
        chain = Chain.Chain(chain_id)
        offset = np.array(translation, dtype=float)
        atom_serial = 1

        for idx, aa in enumerate(sequence, start=1):
            residue = Residue.Residue((" ", idx, " "), _three_letter(aa), " ")
            base = offset + np.array([(idx - 1) * 3.0, 0.0, 0.0], dtype=float)
            ca_bfactor = 0.0 if ca_bfactors is None else ca_bfactors[idx - 1]
            residue.add(Atom.Atom("N", base + np.array([0.0, 0.0, 0.0]), 0.0, 1.0, " ", "N", atom_serial, "N"))
            atom_serial += 1
            residue.add(Atom.Atom("CA", base + np.array([1.0, 0.0, 0.0]), ca_bfactor, 1.0, " ", "C", atom_serial, "C"))
            atom_serial += 1
            residue.add(Atom.Atom("C", base + np.array([2.0, 0.0, 0.0]), 0.0, 1.0, " ", "C", atom_serial, "C"))
            atom_serial += 1
            residue.add(Atom.Atom("O", base + np.array([2.0, 1.0, 0.0]), 0.0, 1.0, " ", "O", atom_serial, "O"))
            atom_serial += 1
            chain.add(residue)

        return chain

    return _make_chain


@pytest.fixture
def sample_chain(make_chain) -> Chain.Chain:
    return make_chain("AGS")
