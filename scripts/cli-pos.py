import argparse
import sys
import pandas as pd
import config as config
import myio as io
import rmsd
import plddt
import cli
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="CLI for positional statistics.")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p",  required=True)
    parser.add_argument("--stat", "-s", choices=["rmsd","plddt"], required=True)
    parser.add_argument("--range", type=cli.parse_range, required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)
    
    refs = io.load_refs(base_path, args.protein)

    #TODO: GPT thought the indexing was confusing (I tend to agree), here we assume a 0-based index, which is then implicitely converted to a one-based index, perhaps using a one-based index in the interface would make more sense?
    start, end = args.range

    dfs = []
    for algorithm in config.algorithms():
        preds = io.load_preds(refs, base_path, args.protein, algorithm)
        for ref in refs.keys():
            r = refs[ref]
            p = preds[ref]

            r_chain = io._first_chain(io._first_model(r))
            p_chain = io._first_chain(io._first_model(p))

            df = None
            rows = None
            if args.stat == "rmsd":
                pos_to_rmsd = rmsd.per_residue_rmsd(start, end, r_chain, p_chain)
                rows = [{"pos": pos + 1, "RMSD": rmsd} for pos, rmsd in pos_to_rmsd.items()]
            elif args.stat == "plddt":
                pos_to_plddt = plddt.per_residue_plddt(start, end, r_chain, p_chain)
                rows = [{"pos": pos + 1, "pLDDT": plddt} for pos, plddt in pos_to_plddt.items()]

            df = pd.DataFrame(rows)
            df["Algorithm"] = algorithm
            df["ref"] = ref
            dfs.append(df)
                
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    main()
