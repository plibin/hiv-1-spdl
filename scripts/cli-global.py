import argparse
import sys
import pandas as pd
import config as config
import myio as io
import rmsd
import tm
import cli
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="CLI for positional statistics.")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p",  required=True)
    parser.add_argument("--stat", "-s", choices=["rmsd","tm"], required=True)
    parser.add_argument("--range", type=cli.parse_range, required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)
    
    refs = io.load_refs(base_path, args.protein)

    start, end = args.range

    rows = []
    for algorithm in config.algorithms():
        preds = io.load_preds(refs, base_path, args.protein, algorithm)
        for ref in refs.keys():
            r = refs[ref]
            p = preds[ref]

            r_chain = io._first_chain(io._first_model(r))
            p_chain = io._first_chain(io._first_model(p))

            row = {}
            if args.stat == "rmsd":
                g_rmsd = rmsd.global_rmsd(start, end, r_chain, p_chain)
                row["RMSD"] = g_rmsd
            elif args.stat == "tm":
                g_tm = tm.global_tm(start, end, r_chain, p_chain)
                row["TM"] = g_tm

            row["Algorithm"] = algorithm
            row["ref"] = ref
            
            rows.append(row)
                
    df = pd.DataFrame(row)
    df.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    main()
