import argparse
import sys
import pandas as pd
import config as config
import myio as io
import rmsd
from pathlib import Path

def parse_range(s):
    try:
        start, end = map(int, s.split('-'))
        return start, end
    except Exception:
        raise argparse.ArgumentTypeError("Range must be in START-END format (e.g. 5-10)")

def main():
    parser = argparse.ArgumentParser(description="CLI for positional statistics.")
    parser.add_argument("--base-path", "-b", required=True)
    parser.add_argument("--protein", "-p",  required=True)
    parser.add_argument("--stat", "-s", choices=["rmsd","plddt"], required=True)
    parser.add_argument("--range", type=parse_range, required=True)

    args = parser.parse_args()

    base_path = Path(args.base_path)
    
    refs = io.load_refs(base_path, args.protein)

    start, end = args.range
    
    for algorithm in config.algorithms():
        print("algo:" + algorithm)
        preds = io.load_preds(refs, base_path, args.protein, algorithm)
        for ref in refs.keys():
            print("ref:" + ref)
            r = refs[ref]
            p = preds[ref]
            if args.stat == "rmsd":
                r_chain = io._first_chain(io._first_model(r))
                p_chain = io._first_chain(io._first_model(p))
                df = rmsd.per_residue_rmsd(start, end, r, p)
                #TODO: add algo + structure name to the df
                

    # Export DataFrame as CSV to stdout
    #df.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    main()
