from __future__ import annotations
from pathlib import Path
from typing import Optional, Union
from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain

PathLike = Union[str, Path]

def parse_structure(path: PathLike) -> Structure.Structure:
    """Load a PDB/mmCIF file to a Biopython Structure.

    Chooses parser from file suffix. Quiet by default.
    """
    p = Path(path)
    if p.suffix.lower() == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        parser = MMCIFParser(QUIET=True)
    return parser.get_structure(p.stem, str(p))


def _first_model(structure: Structure.Structure) -> Model.Model:
    return next(structure.get_models())


def _first_chain(model: Model.Model) -> Chain.Chain:
    return next(model.get_chains())


def select_chain(
    obj: Union[Structure.Structure, Model.Model, Chain.Chain],
    chain_id: Optional[str] = None,
) -> Chain.Chain:
    """Return a Chain from Structure/Model/Chain.

    - If `obj` is a Structure: take first model, then matching chain_id or first chain.
    - If `obj` is a Model: take matching chain_id or first chain.
    - If `obj` is a Chain: return it (ignoring chain_id).
    """
    if isinstance(obj, Chain.Chain):
        return obj
    if isinstance(obj, Structure.Structure):
        model = _first_model(obj)
    elif isinstance(obj, Model.Model):
        model = obj
    else:
        raise TypeError(f"Unsupported object type: {type(obj)}")

    if chain_id is None:
        return _first_chain(model)

    for ch in model.get_chains():
        if ch.id == chain_id:
            return ch
    raise ValueError(f"Chain '{chain_id}' not found; available: {[c.id for c in model.get_chains()]}.")
