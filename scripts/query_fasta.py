from __future__ import annotations

import argparse
from pathlib import Path

import sys

import core
import myio as io


#Extract the query, from the PDB (TODO: make very clear this is different from pdb_fasta)
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p", required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)

    refs = io.load_refs(base_path, args.protein)

    for id_ in refs:
        ref_chain = io._first_chain(io._first_model(refs[id_]))
        chain_id = ref_chain.get_id()
        if chain_id in refs[id_].seqres:
            seq = refs[id_].seqres[chain_id]
        else:
            print(f"WARNING: no SEQRES for {id_} chain {chain_id}, falling back to ATOM sequence", file=sys.stderr)
            seq = core.aa_seq(core.aa_residues(ref_chain))
        print(">" + id_.upper() + "")
        print(seq)

if __name__ == "__main__":
    main()

