from __future__ import annotations

import argparse
from pathlib import Path

import core
import myio as io


"""Build `pdb.fasta` from the reference structures for one protein set.

What this script does (and why):
1) Load reference structures under `<base-path>/<protein>/refs`.
   Why: these are the experimentally observed structures used as ground truth.
2) Select the first model and first chain from each structure.
   Why: the downstream pipeline compares one chain per structure.
3) Keep only amino-acid residues and convert them to a one-letter sequence.
   Why: reference PDB files can contain waters/ligands and may miss residues; we
   want the actual amino-acid sequence present in the structure.
4) Print each sequence as FASTA with header `>PDB_ID_pdb`.
   Why: `>PDB_ID_pdb` headers are later paired with query-sequence headers in `alignment.fasta`.
"""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from reference structures for one protein dataset."
    )
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
