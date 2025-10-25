from __future__ import annotations
from pathlib import Path
from typing import Optional, Union
from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain

def parse_structure_pdb(path: Path) -> Structure.Structure:
    #TODO: check the extension
    parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))

def parse_structure_mmcif(path: Path) -> Structure.Structure:
    #TODO: check the extension
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))

def _first_model(structure: Structure.Structure) -> Model.Model:
    return next(structure.get_models())

def _first_chain(model: Model.Model) -> Chain.Chain:
    return next(model.get_chains())

#load the references, for one protein
def load_refs(base_path, protein):
    refs = {}
    refs_path = base_path / protein / "refs"
    for p in refs_path.iterdir():
        if p.is_file() and p.suffix.lower() == ".pdb":
            refs[p.stem] = parse_structure_pdb(p)

    return refs

#load the predictions, for one protein, for one algorithm
def load_preds(refs, base_path, protein, algorithm):
    preds = {}
    for ref in refs.keys():
        algo_path = base_path / protein / algorithm 
        matches = []
        for p in algo_path.iterdir():
            if p.stem.lower().startswith(ref.lower()):
                matches.append(p)
        if len(matches) == 0 :
            #TODO: should be a msg to std.err (because otherwise it will crash for AlphaFold3)?
            raise RuntimeError("No prediction found for ref " + ref + " for algorithm " + algorithm)
        elif len(matches) > 1 :
            #TODO: should be a msg to std.err (because otherwise it will crash for AlphaFold3)?
            raise RuntimeError("More then one prediction found for ref " + ref + " for algorithm " + algorithm)
        else :
            preds[ref] = parse_structure_pdb(p)

    return preds
    

#TODO: do we need this? remove?
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
