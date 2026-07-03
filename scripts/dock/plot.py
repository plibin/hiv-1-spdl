import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#TODO: we have this in config (in the main repo), no?
ALGORITHM_COLORS = {
    "Reference": "#7f7f7f",  # grey
    "AlphaFold2": "#1f77b4",  # muted blue
    "AlphaFold3": "#ff7f0e",  # orange
    "ESMFold": "#2ca02c",  # green
    "ESM3-Open": "#d62728",  # red
    "ESM3-Large": "#e377c2",  # pink
    "Ember3D": "#9467bd",  # purple
    "AF2-overlap-AF3": "#f1c40f",  # yellow
}

def load_csv(csv_file, algorithm):
    df = pd.read_csv(csv_file)
    df["Algorithm"] = algorithm
    return df

def plot_binding_affinities(df, png_fn):
    fig, ax = plt.subplots(figsize=(10, 6))
    order = [a for a in ALGORITHM_COLORS if a in df["Algorithm"].unique()]
    palette = {a: ALGORITHM_COLORS[a] for a in order}

    sns.boxplot(data=df, x='Algorithm', y='score', order=order, palette=palette, ax=ax, hue='Algorithm', legend=False)

    ax.set_ylabel('Mean Binding Affinity (kcal/mol)', fontsize=14)
    ax.set_xlabel('', fontsize=14)
    ax.tick_params(axis="x", rotation=25, labelsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    plt.savefig(png_fn, dpi=500, bbox_inches="tight")
    plt.close()

def squared_diffs_of_means(refs, df):
    ref_means = (
        refs
        .groupby("name", as_index=False)["score"]
        .mean()
        .rename(columns={"score": "ref_mean_score"})
    )

    df_means = (
        df
        .groupby(["name", "Algorithm"], as_index=False)["score"]
        .mean()
        .rename(columns={"score": "mean_score"})
    )

    df_means = df_means.merge(ref_means, on="name", how="inner")

    df_means["squared_diff"] = (
        df_means["mean_score"] - df_means["ref_mean_score"]
    ) ** 2

    return df_means

def plot_squared_errors(df, png_fn):
    fig, ax = plt.subplots(figsize=(10, 6))

    order = [a for a in ALGORITHM_COLORS if a in df["Algorithm"].unique()]
    palette = {a: ALGORITHM_COLORS[a] for a in order}

    sns.boxplot(
        data=df,
        x="Algorithm",
        y="squared_diff",
        order=order,
        palette=palette,
        ax=ax,
        hue="Algorithm",
        legend=False,
    )

    ax.set_ylabel("Squared Difference of Mean Affinity (kcal/mol)$^2$", fontsize=14)
    ax.set_xlabel("", fontsize=14)
    ax.tick_params(axis="x", rotation=25, labelsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    plt.savefig(png_fn, dpi=500, bbox_inches="tight")
    plt.close()

def make_af2_overlap_with_af3(af2_df, af3_df):
    af3_names = set(af3_df["name"])

    af2_overlap = af2_df[af2_df["name"].isin(af3_names)].copy()
    af2_overlap["Algorithm"] = "AF2-overlap-AF3"

    return af2_overlap

if __name__ == "__main__":
    refs = load_csv("refs.docking.csv", "Reference")
    af2 = load_csv("AF2.docking.csv", "AlphaFold2")
    af3 = load_csv("AF3.docking.csv", "AlphaFold3")
    af2_overlap = make_af2_overlap_with_af3(af2, af3)

    data_df = pd.concat([refs, af2, af3], ignore_index=True)
    data_df_overlap = pd.concat([refs, af2, af3, af2_overlap], ignore_index=True)

    plot_binding_affinities(data_df, 'binding_affinities_boxplot.png')
    plot_binding_affinities(data_df_overlap, 'binding_affinities_overlap_boxplot.png')

    squared_df = pd.concat(
        [
            squared_diffs_of_means(refs, af2),
            squared_diffs_of_means(refs, af3),
        ],
        ignore_index=True,
    )

    plot_squared_errors(squared_df, 'binding_affinities_squared_error_boxplot.png')
