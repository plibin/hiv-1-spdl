import argparse
import sys
import pandas as pd
import config as config
import myio as io
import rmsd
import tm
from pathlib import Path

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="CLI for global statistics.")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p", required=True)
    parser.add_argument("--stat", "-s", choices=["rmsd", "tm"], required=True)
    parser.add_argument("--alignment", "-a", required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)

    refs = io.load_refs(base_path, args.protein)

    aligned_seqs = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))

    rows = []
    for algorithm in config.algorithms():
        preds, _ = io.load_preds(refs, base_path, args.protein, algorithm)
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

            row = {}
            if args.stat == "rmsd":
                g_rmsd = rmsd.global_rmsd(ref, r_align, p_align, r_chain, p_chain)
                row["RMSD"] = g_rmsd
            elif args.stat == "tm":
                g_tm = tm.global_tm(ref, r_align, p_align, r_chain, p_chain)
                row["TM"] = g_tm

            row["Algorithm"] = algorithm
            row["ref"] = ref

            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
