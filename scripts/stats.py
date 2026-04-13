#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scikit_posthocs as sp
from scipy.stats import f_oneway, kruskal, shapiro, levene
from statsmodels.stats.multicomp import pairwise_tukeyhsd

VALUE_COLUMNS = {
    "rmsd": "RMSD",
    "plddt": "pLDDT",
    "grmsd": "RMSD",
    "tm": "TM",
}

ALPHA = 0.05
_SHAPIRO_MAX_N = 5000

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


def _check_normality(groups: dict[str, np.ndarray], alpha: float = ALPHA) -> bool:
    """Return True when every group passes Shapiro-Wilk at the given α.

    Groups with fewer than 3 observations are skipped (too few to test).
    Groups with constant values (zero variance) are treated as non-normal.
    Groups larger than _SHAPIRO_MAX_N are subsampled to avoid Shapiro-Wilk's
    over-sensitivity at very high n.
    """
    rng = np.random.default_rng(42)
    for vals in groups.values():
        n = len(vals)
        if n < 3:
            continue
        if vals.std() < 1e-12:
            return False
        test_vals = vals if n <= _SHAPIRO_MAX_N else rng.choice(vals, _SHAPIRO_MAX_N, replace=False)
        _, p = shapiro(test_vals)
        if p < alpha:
            return False
    return True


def _check_equal_variance(groups: dict[str, np.ndarray], alpha: float = ALPHA) -> bool:
    """Return True when Levene's test (median variant) does not reject H₀.

    Groups with fewer than 2 observations or zero variance are excluded
    before running the test.  Returns True when fewer than 2 testable
    groups remain (nothing to compare).
    """
    arrays = [v for v in groups.values() if len(v) >= 2 and v.std() > 1e-12]
    if len(arrays) < 2:
        return True
    _, p = levene(*arrays)
    return p >= alpha


def _eta_squared(groups: dict[str, np.ndarray]) -> float:
    """η² = SS_between / SS_total  (effect size for one-way ANOVA).

    A p-value only tells us whether group differences are *detectable*; with
    large positional datasets (thousands of rows per algorithm) even trivial
    differences become highly significant.  η² quantifies *how much* of the
    total variance is explained by algorithm choice, giving a practical
    measure of importance:

        η² ≈ 0.01  →  small effect
        η² ≈ 0.06  →  medium effect
        η² ≈ 0.14+ →  large effect [1]

        [1] "Statistical Power Analysis for the Behavioral Sciences" by Jacob Cohen, 2nd edition, 1988, Lawrence Erlbaum Associates

    This lets us distinguish "statistically significant but negligible"
    results from genuinely meaningful performance gaps between algorithms.
    """
    all_vals = np.concatenate(list(groups.values()))
    grand_mean = all_vals.mean()
    ss_total = np.sum((all_vals - grand_mean) ** 2)
    if ss_total == 0:
        return np.nan
    ss_between = sum(len(v) * (v.mean() - grand_mean) ** 2 for v in groups.values())
    return ss_between / ss_total


def _epsilon_squared(h_stat: float, n_total: int, k: int) -> float:
    """ε² = (H − k + 1) / (N − k)  (effect size for Kruskal-Wallis).

    Analogous to η² but designed for the non-parametric setting.  Because
    Kruskal-Wallis operates on ranks rather than raw values, a dedicated
    effect-size measure is needed.  ε² ranges from 0 (no effect) to 1
    (perfect separation) and uses the same interpretive thresholds as η²:

        ε² ≈ 0.01  →  small effect
        ε² ≈ 0.06  →  medium effect
        ε² ≈ 0.14+ →  large effect [1]

    Reporting ε² alongside the Kruskal-Wallis H statistic is essential for
    our datasets: positional CSVs contain thousands of residue-level rows,
    so the test will almost always yield p ≈ 0.  The effect size reveals
    whether the detected difference is scientifically meaningful or merely
    an artifact of the large sample size.
    """
    denom = n_total - k
    if denom <= 0:
        return np.nan
    return (h_stat - k + 1) / denom


def anova_table(df: pd.DataFrame, value_col: str, normality: bool, equal_var: bool) -> pd.DataFrame:
    groups = {name: g[value_col].to_numpy() for name, g in df.groupby("Algorithm") if len(g)}

    if len(groups) < 2:
        return pd.DataFrame([{"test": pd.NA, "statistic": pd.NA, "p_value": pd.NA, "effect_size": pd.NA, "normality": pd.NA, "equal_variance": pd.NA}])

    arrays = list(groups.values())

    if normality and equal_var:
        stat, p_value = f_oneway(*arrays)
        test = "One-way ANOVA"
        effect = _eta_squared(groups)
    else:
        stat, p_value = kruskal(*arrays)
        test = "Kruskal-Wallis"
        n_total = sum(len(v) for v in arrays)
        effect = _epsilon_squared(stat, n_total, len(arrays))

    return pd.DataFrame([{
        "test": test,
        "statistic": stat,
        "p_value": p_value,
        "effect_size": effect,
        "normality": normality,
        "equal_variance": equal_var,
    }])


def posthoc_table(df: pd.DataFrame, value_col: str, normality: bool, equal_var: bool) -> pd.DataFrame:
    if df["Algorithm"].nunique() < 2:
        return pd.DataFrame(columns=["test", "group1", "group2", "meandiff", "p-adj", "reject",])

    if normality and equal_var:
        return _tukey_table(df, value_col)
    else:
        return _dunn_table(df, value_col)


def _tukey_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df["Algorithm"], alpha=ALPHA)
    out = pd.DataFrame(tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    out.insert(0, "test", "Tukey HSD")
    return out


def _dunn_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Pairwise post-hoc comparisons using Dunn's test (Bonferroni-corrected).

    Dunn's test is the non-parametric counterpart of Tukey HSD and is used
    when the ANOVA assumptions (normality and/or equal variance) are violated,
    which is the case for all of our current datasets.  It compares every pair
    of algorithms by examining the differences in their mean ranks (the same
    ranks used by Kruskal-Wallis), making it a natural follow-up to a
    significant Kruskal-Wallis result.

    Bonferroni correction is applied to control the family-wise error rate
    across all (k choose 2) = 15 pairwise comparisons (for k = 6 algorithms).

    The symmetric p-value matrix returned by scikit-posthocs is converted to a
    long-form table that mirrors the Tukey HSD layout so that downstream
    consumers (plots, reports) can handle both formats uniformly.  The
    ``meandiff`` column reports the difference of raw group means (not rank
    means) to keep the output directly interpretable.
    """
    p_matrix = sp.posthoc_dunn(
        df, val_col=value_col, group_col="Algorithm", p_adjust="bonferroni",
    )

    # Convert the symmetric p-value matrix to a long-form table that
    # mirrors the Tukey HSD output.
    groups = sorted(p_matrix.columns)
    rows = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1:]:
            p_adj = p_matrix.loc[g1, g2]
            mean1 = df.loc[df["Algorithm"] == g1, value_col].mean()
            mean2 = df.loc[df["Algorithm"] == g2, value_col].mean()
            rows.append({
                "test": "Dunn",
                "group1": g1,
                "group2": g2,
                "meandiff": round(mean2 - mean1, 4),
                "p-adj": p_adj,
                "reject": p_adj < ALPHA,
            })
    return pd.DataFrame(rows)


def _compute_assumptions(df: pd.DataFrame, value_col: str) -> tuple[bool, bool]:
    groups = {name: g[value_col].to_numpy() for name, g in df.groupby("Algorithm") if len(g)}
    if len(groups) < 2:
        return True, True
    return _check_normality(groups), _check_equal_variance(groups)


def save_tables(df: pd.DataFrame, value_col: str, output_prefix: Path) -> None:
    normality, equal_var = _compute_assumptions(df, value_col)

    mean_table(df, value_col).to_csv(output_prefix.with_name(output_prefix.name + "-means.csv"), index=False)
    anova_table(df, value_col, normality, equal_var).to_csv(output_prefix.with_name(output_prefix.name + "-anova.csv"), index=False)
    posthoc_table(df, value_col, normality, equal_var).to_csv(output_prefix.with_name(output_prefix.name + "-posthoc.csv"), index=False)


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
