from __future__ import annotations

import numpy as np
import copy
from typing import List, Sequence, Callable
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import Chain, Residue, Atom, Superimposer
import core

import argparse
import sys
import pandas as pd
import config as config
import myio as io
import rmsd
import plddt
import cli
from pathlib import Path

#TODO: add a few lines of doc
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p", required=True)

    args = parser.parse_args()
    
    base_path = Path(args.base_path)

    refs = io.load_refs(base_path, args.protein)

    for id_ in refs:
        ref_chain = io._first_chain(io._first_model(refs[id_]))
        ref_res = core.aa_residues(ref_chain)
        print(">" + id_.upper() + "_pdb")
        print(core.aa_seq(ref_res))

if __name__ == "__main__":
    main()

