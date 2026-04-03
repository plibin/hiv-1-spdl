from __future__ import annotations

import argparse
from pathlib import Path

import sys

import core
import myio as io

"""Build `query.fasta` from the reference structures for one protein set.

This script is the counterpart to `pdb_fasta.py`, but extracts a *different*
sequence from the same reference structures:

- `pdb_fasta.py`   → uses the **experimentally observed** ATOM sequence: only
  residues that are actually present (have coordinates) in the structure.
- `query_fasta.py` → uses the **SEQRES** record: the full, intended sequence of
  the protein as deposited by the authors, which may include residues not
  resolved in the ATOM sequence.

Both outputs are combined in `alignment.fasta` and passed to a multiple-sequence
aligner.  The resulting alignment is used downstream to map residue positions
between the reference structures and the model predictions, so it is important
that both sequences originate from the same PDB entry.

What this script does (and why):
1) Load reference structures under `<base-path>/<protein>/refs`.
   Why: the same set of reference structures used as ground truth throughout
   the pipeline.
2) Select the first model and first chain from each structure.
   Why: the downstream pipeline compares one chain per structure.
3) Read the SEQRES sequence for that chain; fall back to the ATOM sequence
   (with a warning) if no SEQRES record is present.
   Why: SEQRES represents the full biological sequence the structure was meant
   to capture, making it a better query for aligning model predictions that
   were generated from a complete sequence rather than a partial structure.
4) Print each sequence as FASTA with header `>PDB_ID` (no `_pdb` suffix).
   Why: plain `>PDB_ID` headers identify the query sequences in
   `alignment.fasta`, distinguishing them from the `>PDB_ID_pdb` headers
   produced by `pdb_fasta.py`.
"""


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
