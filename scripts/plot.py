import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
from matplotlib.pyplot import ylabel

SECONDARY_STRUCTURES = {
    "PR": [(13, 19, 'sheet'), (25, 28, 'helix'), (45, 55, 'sheet'), (59, 64, 'sheet'),
           (67, 72, 'sheet'), (76, 81, 'sheet'), (84, 90, 'sheet'), (95, 99, 'helix')],
    "IN": [(55, 63, 'sheet'), (70, 82, 'helix'),
           (101, 112, 'helix'), (121, 133, 'sheet'), (138, 147, 'sheet')],
    "RT": [(14, 24, 'helix'), (38, 47, 'helix'), (57, 66, 'sheet'), (71, 81, 'sheet'),
           (94, 105, 'helix'), (119, 128, 'helix'), (145, 154, 'sheet'), (174, 183, 'helix'),
           (186, 191, 'sheet'), (194, 200, 'sheet'), (204, 209, 'helix'),
           (212, 220, 'sheet'), (229, 238, 'helix'), (241, 246, 'sheet'),
           (250, 259, 'helix'), (265, 269, 'sheet'), (273, 278, 'helix'),
           (283, 295, 'helix'), (298, 304, 'sheet'), (312, 321, 'helix'),
           (333, 339, 'helix'), (342, 349, 'sheet'), (351, 355, 'sheet'),
           (361, 368, 'helix'), (371, 378, 'helix'), (382, 390, 'helix'),
           (399, 408, 'helix'), (410, 416, 'sheet'), (423, 430, 'sheet'),
           (437, 448, 'helix'), (454, 461, 'sheet'), (468, 479, 'helix'),
           (487, 494, 'helix'), (500, 506, 'sheet'), (512, 519, 'helix'),
           (523, 529, 'helix'), (534, 543, 'sheet'), (547, 553, 'sheet')],
}

# Pastel colors for secondary-structure background regions.
HELIX_BG_COLOR = "#cfe8ff"   # pastel light blue
SHEET_BG_COLOR = "#ffd6d6"   # pastel light red


def _infer_protein_from_positions(df: pd.DataFrame) -> str:
    max_pos = float(df["pos"].max())
    if max_pos <= 110:
        return "PR"
    if max_pos <= 250:
        return "IN"
    return "RT"


def _plot_secondary_structure_background(df: pd.DataFrame, protein: str) -> None:
    if protein not in SECONDARY_STRUCTURES:
        return

    ax = plt.gca()
    for start, end, kind in SECONDARY_STRUCTURES[protein]:
        color = HELIX_BG_COLOR if kind == "helix" else SHEET_BG_COLOR
        ax.axvspan(start, end, ymin=0, ymax=1, color=color, alpha=0.55, zorder=0)


def _add_secondary_structure_legend(ax) -> None:
    handles, labels = ax.get_legend_handles_labels()
    sec_handles = [
        Patch(facecolor=HELIX_BG_COLOR, edgecolor="none", alpha=0.55, label="Alpha helix"),
        Patch(facecolor=SHEET_BG_COLOR, edgecolor="none", alpha=0.55, label="Beta sheet"),
    ]
    ax.legend(handles + sec_handles, labels + ["Alpha helix", "Beta sheet"])


def plot_rmsd(df, protein: str, plot_sec_struct: bool = True):
    if plot_sec_struct:
        _plot_secondary_structure_background(df, protein)
    ax = sns.lineplot(data=df, x="pos", y="RMSD", hue="Algorithm")
    if plot_sec_struct:
        _add_secondary_structure_legend(ax)
    

def plot_plddt(df, protein: str, plot_sec_struct: bool = True):
    if plot_sec_struct:
        _plot_secondary_structure_background(df, protein)

    # ESM3 exports pLDDT on a 0-1 scale; convert to 0-100 for plotting parity.
    df = df.copy()
    esm_mask = df["Algorithm"].isin(["ESM3-Open", "ESM3-Large"])
    df.loc[esm_mask, "pLDDT"] = df.loc[esm_mask, "pLDDT"] * 100

    ax = sns.lineplot(data=df, x="pos", y="pLDDT", hue="Algorithm")
    if plot_sec_struct:
        _add_secondary_structure_legend(ax)
    

def plot_grmsd(df):
    sns.boxplot(data=df, x="Algorithm", y="RMSD")
    

def plot_tm(df):
    sns.boxplot(data=df, x="Algorithm", y="TM")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_path", type=Path, help="Path to the CSV file", required=True)
    parser.add_argument("--type", choices=["rmsd", "plddt", "grmsd", "tm"], required=True,
                        help="Type of plot to generate")
    parser.add_argument("--protein", choices=["PR", "IN", "RT"], required=True,
                        help="Protein for secondary-structure overlays in line plots")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    if args.type == "rmsd":
        plot_rmsd(df, plot_sec_struct=True, protein=args.protein)
    elif args.type == "plddt":
        plot_plddt(df, plot_sec_struct=True, protein=args.protein)
    elif args.type == "grmsd":
        plot_grmsd(df)
    elif args.type == "tm":
        plot_tm(df)

    plt.savefig(args.csv_path.stem + ".png", format="png", dpi=300)


if __name__ == "__main__":
    main()
