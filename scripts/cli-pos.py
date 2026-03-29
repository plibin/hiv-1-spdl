import argparse
import sys
from pathlib import Path

import pandas as pd
import config as config
import myio as io
import rmsd
import plddt

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="CLI for positional statistics.")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p", required=True)
    parser.add_argument("--stat", "-s", choices=["rmsd", "plddt"], required=True)
    parser.add_argument("--alignment", "-a", required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)

    refs = io.load_refs(base_path, args.protein)

    aligned_seqs = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))

    dfs = []
    for algorithm in config.algorithms():
        preds = io.load_preds(refs, base_path, args.protein, algorithm)
        for ref in refs.keys():
            if ref not in preds:
                print(f"Skipping {algorithm}/{ref}: no prediction found", file=sys.stderr)
                continue
            if ref not in aligned_seqs or ref + "_pdb" not in aligned_seqs:
                print(f"Skipping {algorithm}/{ref}: missing alignment record", file=sys.stderr)
                continue

            r = refs[ref]
            p = preds[ref]

            r_chain = io._first_chain(io._first_model(r))
            p_chain = io._first_chain(io._first_model(p))

            r_align = aligned_seqs[ref + "_pdb"].seq
            p_align = aligned_seqs[ref].seq

            rows = None
            if args.stat == "rmsd":
                pos_to_rmsd = rmsd.per_residue_rmsd(ref, r_align, p_align, r_chain, p_chain)
                rows = [{"pos": pos + 1, "RMSD": rmsd_value} for pos, rmsd_value in pos_to_rmsd.items()]
            elif args.stat == "plddt":
                pos_to_plddt = plddt.per_residue_plddt(ref, r_align, p_align, r_chain, p_chain)
                rows = [{"pos": pos + 1, "pLDDT": plddt_value} for pos, plddt_value in pos_to_plddt.items()]

            df = pd.DataFrame(rows)
            df["Algorithm"] = algorithm
            df["ref"] = ref
            dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
