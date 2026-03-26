from __future__ import annotations

import sys
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain


def parse_structure_pdb(path: Path) -> Structure.Structure:
    if path.suffix.lower() != ".pdb":
        raise ValueError(f"Expected .pdb file, got {path}")
    parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def parse_structure_mmcif(path: Path) -> Structure.Structure:
    if path.suffix.lower() not in [".cif", ".mmcif"]:
        raise ValueError(f"Expected .cif or .mmcif file, got {path}")
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def _first_model(structure: Structure.Structure) -> Model.Model:
    return next(structure.get_models())


def _first_chain(model: Model.Model) -> Chain.Chain:
    return next(model.get_chains())


# load the references, for one protein
def load_refs(base_path: str, protein: str) -> dict[str, Structure.Structure]:
    refs = {}
    refs_path = base_path / protein / "refs"
    for p in refs_path.iterdir():
        if p.is_file():
            if p.suffix.lower() == ".pdb":
                refs[p.stem.upper()] = parse_structure_pdb(p)
            # elif p.suffix.lower() in [".cif", ".mmcif"]:
            #     refs[p.stem] = parse_structure_mmcif(p)

    return refs


# load the predictions, for one protein, for one algorithm
def load_preds(refs: dict[str, Structure.Structure], base_path: str, protein: str, algorithm: str) -> dict[str, Structure.Structure]:
    preds = {}
    algo_path = base_path / protein / algorithm
    for ref in refs.keys():
        matches = []
        for p in algo_path.iterdir():
            if p.stem.lower().startswith(ref.lower()) and p.suffix.lower() == ".pdb":
                matches.append(p)
        if len(matches) == 0:
            print("No prediction found for ref " + ref + " for algorithm " + algorithm, file=sys.stderr)
        elif len(matches) > 1:
            raise RuntimeError("More than one prediction found for ref " + ref + " for algorithm " + algorithm)
        else:
            preds[ref] = parse_structure_pdb(matches[0])

    return preds
