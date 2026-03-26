from __future__ import annotations

import argparse
from pathlib import Path

import core
import myio as io


#TODO: add a few lines of doc -> get a fasta of the AAs that are actually in the PDB
def main():
    parser = argparse.ArgumentParser(description="CLI for extracting FASTA sequences from PDB files.")
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
