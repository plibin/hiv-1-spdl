from __future__ import annotations

import argparse
from pathlib import Path

import sys

import core
import myio as io

from Bio import SeqIO

def load_overrides(path: str | None) -> dict[str, str]:
    if path is None:
        return {}
    return {
        rec.id.upper(): str(rec.seq).upper()
        for rec in SeqIO.parse(path, "fasta")
    }

#Extract the query, from the PDB (TODO: make very clear this is different from pdb_fasta)
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p", required=True)
    parser.add_argument("--overrides_fasta", "-o", default=None)

    args = parser.parse_args()

    overrides = {}
    if args.overrides_fasta is not None:
        overrides = load_overrides(args.overrides_fasta)

    base_path = Path(args.base_path)

    refs = io.load_refs(base_path, args.protein)

    for id_ in refs:
        ref_chain = io._first_chain(io._first_model(refs[id_]))
        chain_id = ref_chain.get_id()
        if chain_id in refs[id_].seqres:
            seq = refs[id_].seqres[chain_id]
        else:
            if id_ not in overrides:
                raise ValueError(f"No SEQRES for {id_} chain {chain_id} and no override provided")
            
            #"No SEQRES for this id's chain, look one up in the overrides FASTA
            seq = overrides[id_]
        print(">" + id_.upper() + "")
        print(seq)

if __name__ == "__main__":
    main()

