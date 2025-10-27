import seaborn as sns
import matplotlib.pyplot as plt

import argparse
import pandas as pd
from pathlib import Path
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_path", type=Path, help="Path to the CSV file")
    parser.add_argument("--type", choices=["rmsd","plddt", "grmsd", "tm"], required=True)
    args = parser.parse_args()

    print (args.csv_path)
    df = pd.read_csv(args.csv_path)
    if args.type == "rmsd":
        sns.lineplot(data=df, x="pos", y="RMSD", hue="Algorithm")
    elif args.type == "plddt":
        sns.lineplot(data=df, x="pos", y="pLDDT", hue="Algorithm")
    elif args.type == "grmsd":
        sns.boxplot(data=df, x="Algorithm", y="RMSD")
    elif args.type == "tm":
        sns.boxplot(data=df, x="Algorithm", y="TM")
        
    plt.savefig(args.csv_path.stem + ".png", format="png", dpi=300)

if __name__ == "__main__":
    main()

