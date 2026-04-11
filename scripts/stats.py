#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

VALUE_COLUMNS = {
    "rmsd": "RMSD",
    "plddt": "pLDDT",
    "grmsd": "RMSD",
    "tm": "TM",
}

def prepare_data(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    cols = ["Algorithm", value_col]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}. Available: {list(df.columns)}")

    out = df[cols].copy()
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")
    return out.dropna(subset=cols)

def mean_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    return (
        df.groupby("Algorithm", as_index=False)[value_col]
        .mean()
        .rename(columns={value_col: f"Mean {value_col}"})
        .sort_values("Algorithm")
    )

def anova_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    groups = [g[value_col].to_numpy() for _, g in df.groupby("Algorithm") if len(g)]
    if len(groups) < 2:
        return pd.DataFrame([{"F": pd.NA, "p_value": pd.NA}])

    f_stat, p_value = f_oneway(*groups)
    return pd.DataFrame([{"F": f_stat, "p_value": p_value}])

def tukey_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    if df["Algorithm"].nunique() < 2:
        return pd.DataFrame(columns=["group1", "group2", "meandiff", "p-adj", "lower", "upper", "reject"])

    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df["Algorithm"], alpha=0.05)
    return pd.DataFrame(tukey._results_table.data[1:], columns=tukey._results_table.data[0])

def save_tables(df: pd.DataFrame, value_col: str, output_prefix: Path) -> None:
    mean_table(df, value_col).to_csv(output_prefix.with_name(output_prefix.name + "-means.csv"), index=False)
    anova_table(df, value_col).to_csv(output_prefix.with_name(output_prefix.name + "-anova.csv"), index=False)
    tukey_table(df, value_col).to_csv(output_prefix.with_name(output_prefix.name + "-tukey.csv"), index=False)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_path", type=Path, required=True)
    parser.add_argument("--type", choices=VALUE_COLUMNS, required=True)
    parser.add_argument("--output_prefix", type=Path, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    value_col = VALUE_COLUMNS[args.type]
    df = prepare_data(df, value_col)

    output_prefix = args.output_prefix or args.csv_path.with_suffix("")
    save_tables(df, value_col, output_prefix)

if __name__ == "__main__":
    main()
