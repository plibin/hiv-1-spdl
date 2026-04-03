import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Patch

SECONDARY_STRUCTURES = {
    "PR": [(13, 19, 'sheet'), (25, 28, 'helix'), (45, 55, 'sheet'), (59, 64, 'sheet'),
           (67, 72, 'sheet'), (76, 81, 'sheet'), (84, 90, 'sheet'), (95, 99, 'helix')],
    # Consensus secondary structures for IN derived from HELIX/SHEET records of 25 chain-A reference PDBs.
    "IN": [(60, 68, 'sheet'),  # β1
           (71, 78, 'sheet'),  # β2
           (84, 89, 'sheet'),  # β3
           (93, 108, 'helix'),  # α1
           (112, 115, 'sheet'),  # β4
           (118, 122, 'helix'),  # 3₁₀ turn
           (123, 134, 'helix'),  # α2
           (136, 139, 'sheet'),  # β5
           (153, 166, 'helix'),  # α3
           (167, 169, 'helix'),  # 3₁₀ turn
           (171, 186, 'helix'),  # α4
           (195, 208, 'helix')],  # α5
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
HELIX_BG_COLOR = "#cfe8ff"  # pastel light blue
SHEET_BG_COLOR = "#ffd6d6"  # pastel light red


class _SubtitleHandler(HandlerBase):
    """Invisible legend handle used to render subtitle-only rows."""

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        handlebox.set_visible(False)
        return handlebox


def _plot_secondary_structure_background(ax, protein: str) -> None:
    if protein not in SECONDARY_STRUCTURES:
        return

    for start, end, kind in SECONDARY_STRUCTURES[protein]:
        color = HELIX_BG_COLOR if kind == "helix" else SHEET_BG_COLOR
        ax.axvspan(start, end, ymin=0, ymax=1, color=color, alpha=0.55, zorder=0)


def _build_positional_legend(ax, with_sec_struct: bool) -> None:
    """Build a right-side legend spanning the full plot height.

    Entries are grouped under bold subtitles ("Algorithms" and optionally
    "Secondary structures") and vertically centred within the legend frame.
    """
    algo_handles, algo_labels = ax.get_legend_handles_labels()

    handles, labels = [], []
    subtitle_handles = []

    # --- Algorithms section ---
    algo_subtitle = Patch(facecolor="none", edgecolor="none")
    subtitle_handles.append(algo_subtitle)
    handles.append(algo_subtitle)
    labels.append("Algorithms")

    handles.extend(algo_handles)
    labels.extend(algo_labels)

    # --- Secondary structures section ---
    if with_sec_struct:
        sec_subtitle = Patch(facecolor="none", edgecolor="none")
        subtitle_handles.append(sec_subtitle)
        handles.append(sec_subtitle)
        labels.append("Structures")

        handles.append(Patch(facecolor=HELIX_BG_COLOR, edgecolor="none", alpha=0.55))
        labels.append("Alpha helix")
        handles.append(Patch(facecolor=SHEET_BG_COLOR, edgecolor="none", alpha=0.55))
        labels.append("Beta sheet")

    handler_map = {h: _SubtitleHandler() for h in subtitle_handles}

    leg = ax.legend(
        handles, labels,
        bbox_to_anchor=(1.02, 0.0, 0.0, 1.0),
        loc="center left",
        borderaxespad=0,
        handler_map=handler_map,
        frameon=True,
        fontsize=13,
    )

    # Centre entries vertically inside the legend box.
    leg._legend_box.align = "center"

    # Bold the subtitle labels.
    for text in leg.get_texts():
        if text.get_text() in ("Algorithms", "Structures"):
            text.set_fontweight("bold")
            text.set_fontsize(14)


def _format_pos_ax(ax, ylabel):
    ax.set_xlabel("Position", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="both", labelsize=13)
    ax.grid(True, linestyle="--", alpha=0.4, zorder=0)


def plot_rmsd(df, protein: str, plot_sec_struct: bool = True, figsize: tuple = (15, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    if plot_sec_struct:
        _plot_secondary_structure_background(ax, protein)
    sns.lineplot(data=df, x="pos", y="RMSD", hue="Algorithm")
    _build_positional_legend(ax, with_sec_struct=plot_sec_struct)
    _format_pos_ax(ax, "RMSD")
    fig.tight_layout()


def plot_plddt(df, protein: str, plot_sec_struct: bool = True, figsize: tuple = (15, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    if plot_sec_struct:
        _plot_secondary_structure_background(ax, protein)

    # ESM3 exports pLDDT on a 0-1 scale; convert to 0-100 for plotting parity.
    df = df.copy()
    esm_mask = df["Algorithm"].isin(["ESM3-Open", "ESM3-Large"])
    df.loc[esm_mask, "pLDDT"] = df.loc[esm_mask, "pLDDT"] * 100

    sns.lineplot(data=df, x="pos", y="pLDDT", hue="Algorithm")
    _build_positional_legend(ax, with_sec_struct=plot_sec_struct)
    _format_pos_ax(ax, "pLDDT")
    fig.tight_layout()


def plot_grmsd(df):
    fig, ax = plt.subplots()
    sns.boxplot(data=df, x="Algorithm", y="RMSD", ax=ax)
    ax.tick_params(axis="x", rotation=25)


def plot_tm(df):
    fig, ax = plt.subplots()
    sns.boxplot(data=df, x="Algorithm", y="TM", ax=ax)
    ax.tick_params(axis="x", rotation=25)


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

    plt.savefig(args.csv_path.stem + ".png", format="png", dpi=500, bbox_inches="tight")


if __name__ == "__main__":
    main()
