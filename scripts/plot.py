import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from utils import linear_regression_fit, binned_mean_regression_fit


def _weighted_binned_correlation(x, y, n_bins: int = 4):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return np.nan

    xf, yf = x[mask], y[mask]
    bin_edges = np.linspace(xf.min(), xf.max(), n_bins + 1)
    bin_indices = np.digitize(xf, bin_edges, right=True)

    mean_x, mean_y, weights = [], [], []
    for b in range(1, n_bins + 1):
        sel = bin_indices == b
        if sel.any():
            mean_x.append(xf[sel].mean())
            mean_y.append(yf[sel].mean())
            weights.append(sel.sum())

    if len(mean_x) < 2:
        return np.nan

    mean_x = np.asarray(mean_x)
    mean_y = np.asarray(mean_y)
    weights = np.asarray(weights, dtype=float)

    mx = np.average(mean_x, weights=weights)
    my = np.average(mean_y, weights=weights)
    cov = np.average((mean_x - mx) * (mean_y - my), weights=weights)
    vx = np.average((mean_x - mx) ** 2, weights=weights)
    vy = np.average((mean_y - my) ** 2, weights=weights)
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)

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

# The four algorithms shown in the 2×2 pLDDT-vs-RMSD correlation grid.
CORRELATION_ALGORITHMS = ["AlphaFold2", "AlphaFold3", "ESMFold", "Ember3D"]


# Per-algorithm scatter colors for the 2×2 correlation grid.
CORRELATION_COLORS = {
    "AlphaFold2": "#1f77b4",   # muted blue
    "AlphaFold3": "#ff7f0e",   # orange
    "ESMFold":    "#2ca02c",   # green
    "Ember3D":    "#9467bd",   # purple
}


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

    # Exclude ESM3-Open and ESM3-Large from the positional pLDDT plot.
    df = df[~df["Algorithm"].isin(["ESM3-Open", "ESM3-Large"])].copy()

    sns.lineplot(data=df, x="pos", y="pLDDT", hue="Algorithm")
    _build_positional_legend(ax, with_sec_struct=plot_sec_struct)
    _format_pos_ax(ax, "pLDDT")
    fig.tight_layout()


def plot_correlation(df_rmsd: pd.DataFrame, df_plddt: pd.DataFrame, figsize: tuple = (12, 10)):
    # Merge on shared keys.
    df = pd.merge(df_rmsd, df_plddt, on=["pos", "Algorithm", "ref"], how="inner")

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    for ax, algo in zip(axes.flatten(), CORRELATION_ALGORITHMS):
        sub = df[df["Algorithm"] == algo]
        color = CORRELATION_COLORS.get(algo, "#888888")
        ax.scatter(sub["pLDDT"], sub["RMSD"], alpha=0.3, s=10, edgecolors="none", color=color)

        # Linear regression fit line.
        if len(sub) > 1:
            fit = linear_regression_fit(sub["pLDDT"].values, sub["RMSD"].values)
            if fit is not None:
                x_line, y_line, _ = fit
                ax.plot(x_line, y_line, color="red", linewidth=1.5)

            # Binned-mean regression: 4 equal-width pLDDT bins → mean points → fit.
            binned_fit = binned_mean_regression_fit(sub["pLDDT"].values, sub["RMSD"].values, n_bins=4)
            if binned_fit is not None:
                x_bin, y_bin, _, mean_x, mean_y = binned_fit
                ax.plot(x_bin, y_bin, color="black", linewidth=1.5, linestyle="--")
                ax.scatter(mean_x, mean_y, color="black", s=40, zorder=5)

            r = np.corrcoef(sub["pLDDT"].values, sub["RMSD"].values)[0, 1]
            r2 = r ** 2 if np.isfinite(r) else np.nan
            rw = _weighted_binned_correlation(sub["pLDDT"].values, sub["RMSD"].values, n_bins=4)
            ax.text(
                0.03, 0.97,
                f"$R^2$ = {r2:.3f}\n$r_w$ = {rw:.3f}",
                transform=ax.transAxes,
                va="top", ha="left", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="0.8"),
            )

        ax.set_title(algo, fontsize=14, fontweight="bold")
        ax.set_xlabel("pLDDT", fontsize=12)
        ax.set_ylabel("RMSD", fontsize=12)
        ax.tick_params(axis="both", labelsize=11)
        ax.grid(True, linestyle="--", alpha=0.4)

    # Shared legend centred below the bottom row of subplots.
    legend_handles = [
        Line2D([], [], color="red", linewidth=1.5, label="Linear fit"),
        Line2D([], [], color="black", linewidth=1.5, linestyle="--", label="Binned-mean fit"),
        Line2D([], [], marker="o", color="black", linewidth=0, markersize=6, label="Bin means"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               fontsize=12, frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0, 0.04, 1, 1])


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
    parser.add_argument("--csv_path2", type=Path, default=None,
                        help="Second CSV path (pLDDT CSV for --type correlation)")
    parser.add_argument("--type",
                        choices=["rmsd", "plddt", "grmsd", "tm", "correlation"],
                        required=True,
                        help="Type of plot to generate")
    parser.add_argument("--protein", choices=["PR", "IN", "RT"], required=True,
                        help="Protein for secondary-structure overlays in line plots")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output filename (default: <csv_stem>.png)")
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
    elif args.type == "correlation":
        if args.csv_path2 is None:
            parser.error("--csv_path2 (pLDDT CSV) is required for --type correlation")
        df_plddt = pd.read_csv(args.csv_path2)
        plot_correlation(df, df_plddt)

    out = args.output if args.output else (args.csv_path.stem + ".png")
    plt.savefig(out, format="png", dpi=500, bbox_inches="tight")


if __name__ == "__main__":
    main()
