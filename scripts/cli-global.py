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
            p = preds[ref] #TODO: there is not a prediction for every ref (AlphaFold3 has fewer predictions, as some of the sequences were part of the training set!!!), we should check for this (also in cli-pos)

            r_chain = io._first_chain(io._first_model(r))
            p_chain = io._first_chain(io._first_model(p))

            row = {}
            if args.stat == "rmsd":
                g_rmsd = rmsd.global_rmsd(start, end, r_chain, p_chain)
                row["RMSD"] = g_rmsd
            elif args.stat == "tm":
                #TODO: given our problem set this choice for L_N seems to make sense?
                #TODO: GPT said: The original TM-score divides by L_N, the length of the native structure (the reference chain length), not just the number of residues in the evaluated window. You currently pass in L_N = end - start from cli-global.py, with a comment. This is a variant — it makes the TM-score “local-TM over the window.” That’s not inherently wrong if that’s what you want, but it’s not standard TM-score. Just be aware that these values won’t be directly comparable to literature TM-scores unless you use the full native length. -> I would say this makes most sense here, but we need to specify it in the manuscript...
                L_N = end - start
                g_tm = tm.global_tm(start, end, r_chain, p_chain, L_N)
                row["TM"] = g_tm

            row["Algorithm"] = algorithm
            row["ref"] = ref
            
            rows.append(row)
                
    df = pd.DataFrame(rows)
    df.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    main()
