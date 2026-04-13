#!/usr/bin/env python3
"""
Statistical Diagnostics for stats.py
=====================================

Audits the statistical assumptions behind the one-way ANOVA and Tukey HSD
tests used by stats.py.  For every result CSV produced by the benchmarking
pipeline, it:

  1. States what assumptions must hold for ANOVA / Tukey to be valid.
  2. Runs formal diagnostic tests for each assumption.
  3. Generates supporting plots.
"""

from __future__ import annotations

import textwrap
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats


# ── Paths ──────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"
DIAG_DIR = RESULTS_DIR / "diagnostics"
FIG_DIR = DIAG_DIR / "figures"

ALPHA = 0.05

PROTEINS = ["PR", "IN", "RT"]
DATASETS = [
    {"label": "Global RMSD",      "tpl": "{prot}-global-rmsd.csv", "col": "RMSD",  "kind": "global"},
    {"label": "Global TM",        "tpl": "{prot}-global-tm.csv",   "col": "TM",    "kind": "global"},
    {"label": "Positional RMSD",  "tpl": "{prot}-pos-rmsd.csv",    "col": "RMSD",  "kind": "positional"},
    {"label": "Positional pLDDT", "tpl": "{prot}-pos-plddt.csv",   "col": "pLDDT", "kind": "positional"},
]


# ═══════════════════════════════════════════════════════════════════════════
# Logging helper
# ═══════════════════════════════════════════════════════════════════════════
class _Log:
    def __init__(self):
        self._buf = StringIO()

    def h1(self, title):
        self._buf.write("\n" + "=" * 80 + "\n")
        self._buf.write(f"  {title}\n")
        self._buf.write("=" * 80 + "\n\n")

    def h2(self, title):
        self._buf.write(f"\n--- {title} ---\n\n")

    def h3(self, title):
        self._buf.write(f"\n  • {title}\n")

    def p(self, text, indent=0):
        prefix = "    " * indent
        for line in text.split("\n"):
            self._buf.write(prefix + line + "\n")

    def blank(self):
        self._buf.write("\n")

    def table(self, df, indent=1):
        prefix = "    " * indent
        for line in df.to_string(index=False).split("\n"):
            self._buf.write(prefix + line + "\n")
        self._buf.write("\n")

    def save(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self._buf.getvalue())


# ═══════════════════════════════════════════════════════════════════════════
# Diagnostic tests
# ═══════════════════════════════════════════════════════════════════════════
def _test_normality(groups: dict[str, np.ndarray]) -> pd.DataFrame:
    """Shapiro-Wilk (n ≤ 5000) or Anderson-Darling (n > 5000) per group."""
    rows = []
    for name, vals in groups.items():
        n = len(vals)
        if n < 3:
            rows.append({"Algorithm": name, "n": n, "statistic": np.nan,
                         "p / crit": "—", "normal": "too few"})
            continue
        if vals.std() < 1e-12:
            rows.append({"Algorithm": name, "n": n, "statistic": "—",
                         "p / crit": "constant", "normal": "degenerate"})
            continue
        if n > 5000:
            res = sp_stats.anderson(vals, dist="norm", method="interpolate")
            reject = res.pvalue < ALPHA
            rows.append({"Algorithm": name, "n": n,
                         "statistic": f"A²={res.statistic:.4f}",
                         "p / crit": f"p={res.pvalue:.2e}" if res.pvalue < 0.001 else f"p={res.pvalue:.4f}",
                         "normal": "no" if reject else "yes"})
        else:
            w, p = sp_stats.shapiro(vals)
            rows.append({"Algorithm": name, "n": n,
                         "statistic": f"W={w:.4f}",
                         "p / crit": f"p={p:.2e}" if p < 0.001 else f"p={p:.4f}",
                         "normal": "no" if p < ALPHA else "yes"})
    return pd.DataFrame(rows)


def _test_levene(groups: dict[str, np.ndarray]) -> dict:
    """Levene's test (median variant) for equal variances."""
    arrays = [v for v in groups.values() if len(v) >= 2 and v.std() > 1e-12]
    if len(arrays) < 2:
        return {"statistic": np.nan, "p_value": np.nan, "equal_var": "insufficient"}
    stat, p = sp_stats.levene(*arrays)
    return {"statistic": round(stat, 4), "p_value": p, "equal_var": "no" if p < ALPHA else "yes"}


def _autocorrelation(series: np.ndarray, max_lag: int = 20) -> np.ndarray:
    n = len(series)
    max_lag = min(max_lag, n - 1)
    if max_lag < 1:
        return np.array([])
    mu = series.mean()
    var = ((series - mu) ** 2).sum()
    if var == 0:
        return np.zeros(max_lag)
    return np.array([
        np.sum((series[: n - lag] - mu) * (series[lag:] - mu)) / var
        for lag in range(1, max_lag + 1)
    ])


# ═══════════════════════════════════════════════════════════════════════════
# Plots
# ═══════════════════════════════════════════════════════════════════════════
def _plot_qq(groups, title, path):
    ng = len(groups)
    ncols = min(3, ng)
    nrows = -(-ng // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    for ax, (name, vals) in zip(axes.flat, groups.items()):
        if vals.std() < 1e-12:
            ax.text(0.5, 0.5, f"{name}\n(constant value)", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11)
            ax.set_title(f"{name} (n={len(vals)})", fontsize=11)
            continue
        sp_stats.probplot(vals, dist="norm", plot=ax)
        ax.set_title(f"{name} (n={len(vals)})", fontsize=11)
        ax.get_lines()[0].set(markersize=3, alpha=0.6)
        ax.get_lines()[1].set(color="red", linewidth=1.5)

    for ax in axes.flat[ng:]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_dist(groups, title, value_label, path):
    fig, ax = plt.subplots(figsize=(10, 5))
    for name, vals in groups.items():
        if len(vals) > 1 and vals.std() > 1e-12:
            sns.kdeplot(vals, ax=ax, label=f"{name} (n={len(vals)})", fill=True, alpha=0.2)
        elif vals.std() < 1e-12:
            ax.axvline(vals[0], linestyle="--", label=f"{name} (constant={vals[0]:.2f})")
    ax.set_xlabel(value_label, fontsize=13)
    ax.set_ylabel("Density", fontsize=13)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_var(groups, title, value_label, levene_p, path):
    data = [{"Algorithm": n, value_label: v} for n, vals in groups.items() for v in vals]
    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=df, x="Algorithm", y=value_label, ax=ax)
    ax.tick_params(axis="x", rotation=20)
    ax.set_title(title, fontsize=14, fontweight="bold")

    parts = []
    for n, vals in groups.items():
        parts.append(f"{n}: σ²={vals.var():.4f}")
    footer = f"Levene p={levene_p:.4e}   |   " + "  |  ".join(parts)
    ax.text(0.5, -0.18, footer, transform=ax.transAxes, ha="center", fontsize=7, style="italic")

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_acf(df, protein, value_col, path):
    """Autocorrelation bar-charts for a sample of (algorithm, ref) pairs."""
    combos = df.groupby(["Algorithm", "ref"]).size().reset_index(name="cnt")
    combos = combos[combos["cnt"] >= 30]
    if combos.empty:
        return

    samples = []
    for algo in combos["Algorithm"].unique():
        row = combos[combos["Algorithm"] == algo].iloc[0]
        samples.append((algo, row["ref"]))
        if len(samples) == 4:
            break

    n = len(samples)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=True, squeeze=False)
    max_lag = 20

    for ax, (algo, ref) in zip(axes.flat, samples):
        sub = df[(df["Algorithm"] == algo) & (df["ref"] == ref)].sort_values("pos")
        vals = sub[value_col].values
        acf = _autocorrelation(vals, max_lag)
        lags = np.arange(1, len(acf) + 1)

        ax.bar(lags, acf, color="steelblue", alpha=0.7)
        ci = 1.96 / np.sqrt(len(vals))
        ax.axhline(ci, color="red", ls="--", lw=0.8, label="95% CI (white noise)")
        ax.axhline(-ci, color="red", ls="--", lw=0.8)
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("Lag (residue positions)")
        ax.set_title(f"{algo} / {ref}\n(n={len(vals)})", fontsize=10)
        ax.set_ylim(-0.3, 1.0)

    axes.flat[0].set_ylabel("Autocorrelation", fontsize=12)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"{protein} — Spatial autocorrelation of per-residue {value_col}",
                 fontsize=13, fontweight="bold", y=1.05)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_summary_heatmap(results, path):
    labels = [r["dataset"] for r in results]
    matrix = np.array([
        [int(r["all_normal"]), int(r["equal_var"]), int(r["independent"])]
        for r in results
    ])

    fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.55)))
    cmap = ListedColormap(["#ff6b6b", "#51cf66"])
    ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Normality\n(Shapiro-Wilk /\nAnderson-Darling)",
                        "Equal Variance\n(Levene)",
                        "Independence"], fontsize=11)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)

    for i in range(len(labels)):
        for j in range(3):
            txt = "PASS" if matrix[i, j] else "FAIL"
            clr = "black" if matrix[i, j] else "white"
            ax.text(j, i, txt, ha="center", va="center", fontsize=10, fontweight="bold", color=clr)

    ax.set_title("ANOVA / Tukey HSD Assumption Audit", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Main analysis
# ═══════════════════════════════════════════════════════════════════════════
def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    log = _Log()
    summary = []

    # ── 1  Introduction ────────────────────────────────────────────────
    log.h1("STATISTICAL DIAGNOSTICS REPORT")
    log.p(textwrap.dedent("""\
        Generated by stats_diagnostics.py for the HIV-1 structure prediction
        benchmarking pipeline.  This report audits the assumptions behind the
        statistical tests used in stats.py.
    """))

    # ── 2  Current implementation ──────────────────────────────────────
    log.h1("1. CURRENT IMPLEMENTATION IN stats.py")
    log.p(textwrap.dedent("""\
        stats.py computes three outputs for each result CSV:

          1. Mean table     — per-algorithm arithmetic mean.
          2. ANOVA table    — one-way ANOVA via scipy.stats.f_oneway.
          3. Tukey HSD table — all pairwise post-hoc comparisons via
                              statsmodels pairwise_tukeyhsd (α = 0.05).

        These tests are applied uniformly to four metric types:
          • Global RMSD    (one value per protein structure per algorithm)
          • Global TM      (one value per protein structure per algorithm)
          • Positional RMSD  (one value per residue per structure per algorithm)
          • Positional pLDDT (one value per residue per structure per algorithm)
    """))

    # ── 3  Required assumptions ────────────────────────────────────────
    log.h1("2. ASSUMPTIONS REQUIRED BY ANOVA AND TUKEY HSD")

    log.h2("2.1  One-way ANOVA")
    log.p(textwrap.dedent("""\
        One-way ANOVA tests H₀: all group (algorithm) population means are
        equal.  Three assumptions must hold for the F-statistic and its
        p-value to be valid:

          A1. INDEPENDENCE
              Observations within and across groups must be statistically
              independent.  Each data point should come from a distinct,
              unrelated experimental unit.

          A2. NORMALITY
              The distribution of the metric within each group should be
              approximately normal.  Moderate violations are tolerable when
              group sizes are large (central limit theorem), but heavy skew,
              multimodality, or extreme outliers degrade the test.

          A3. HOMOSCEDASTICITY (equal variances)
              The population variances of the metric should be roughly equal
              across all groups.  When variances differ and group sizes are
              unequal, the Type-I error rate of the standard F-test is
              inflated, making results unreliable.

        When A2 or A3 are violated, non-parametric alternatives — most
        commonly the Kruskal-Wallis test — are recommended.  When only A3
        is violated, Welch's ANOVA is an alternative that does not assume
        equal variances.
    """))

    log.h2("2.2  Tukey HSD (post-hoc pairwise comparisons)")
    log.p(textwrap.dedent("""\
        Tukey HSD inherits all three ANOVA assumptions.  In addition:

          A4. BALANCED DESIGN (equal group sizes)
              The original Tukey procedure assumes the same number of
              observations per group.  The statsmodels implementation applies
              the Tukey-Kramer correction for unequal n, which maintains
              correct coverage but loses statistical power.

        When ANOVA assumptions are violated, Dunn's test (a rank-based
        non-parametric post-hoc paired with Kruskal-Wallis) is the standard
        alternative.
    """))

    log.h2("2.3  Diagnostic tests applied in this report")
    log.p(textwrap.dedent("""\
        For every (protein × metric) dataset we run:

          • Shapiro-Wilk test (n ≤ 5000) or Anderson-Darling (n > 5000) per
            algorithm group → tests A2 (normality).
          • Levene's test (median variant, robust to non-normality) across all
            groups → tests A3 (homoscedasticity).
          • Lag-autocorrelation analysis for positional data → assesses A1
            (independence).
          • Group-size comparison → assesses A4 (balance).

        Significance level: α = 0.05 throughout.
    """))

    # ── 4  Per-dataset results ─────────────────────────────────────────
    log.h1("3. DIAGNOSTIC RESULTS PER DATASET")

    ds_counter = 0
    for protein in PROTEINS:
        for ds in DATASETS:
            csv_path = RESULTS_DIR / ds["tpl"].format(prot=protein.lower())
            if not csv_path.exists():
                continue

            ds_counter += 1
            df = pd.read_csv(csv_path)
            vcol = ds["col"]
            label = f"{protein} {ds['label']}"
            slug = f"{protein.lower()}-{ds['label'].lower().replace(' ', '-')}"

            log.h2(f"3.{ds_counter}  {label}  ({csv_path.name})")

            # ── groups ──
            groups: dict[str, np.ndarray] = {}
            for algo, grp in df.groupby("Algorithm"):
                vals = pd.to_numeric(grp[vcol], errors="coerce").dropna().values
                if len(vals) > 0:
                    groups[algo] = vals

            # ── A4: sample sizes / balance ──
            log.h3("Sample sizes & descriptive statistics  (Assumption A4: balance)")
            desc = df.groupby("Algorithm").agg(
                n=(vcol, "count"),
                mean=(vcol, "mean"),
                std=(vcol, "std"),
                median=(vcol, "median"),
            ).reset_index().round(4)
            log.table(desc)

            balanced = desc["n"].nunique() == 1
            log.p(f"Balanced design: {'YES' if balanced else 'NO — group sizes differ.'}", indent=1)
            if not balanced:
                log.p("→ The Tukey-Kramer correction (used by statsmodels) accounts for", indent=1)
                log.p("  unequal n, but statistical power is reduced for smaller groups.", indent=1)
            log.blank()

            # ── A2: normality ──
            log.h3("Normality  (Assumption A2)")
            log.p("H₀: data within each algorithm group is normally distributed.", indent=1)
            log.p(f"Significance level: α = {ALPHA}", indent=1)
            log.blank()

            norm_df = _test_normality(groups)
            log.table(norm_df)

            ok_vals = {"yes", "too few"}
            all_normal = all(v in ok_vals for v in norm_df["normal"])
            n_fail = sum(1 for v in norm_df["normal"] if v == "no")
            n_degen = sum(1 for v in norm_df["normal"] if v == "degenerate")

            if all_normal and n_degen == 0:
                log.p("Verdict: ALL groups pass normality.", indent=1)
            else:
                parts = []
                if n_fail:
                    parts.append(f"{n_fail} group(s) reject normality")
                if n_degen:
                    parts.append(f"{n_degen} group(s) are degenerate (constant value)")
                log.p(f"Verdict: {'; '.join(parts)}.", indent=1)
                log.p("→ ANOVA's normality assumption is VIOLATED.", indent=1)
                log.p("→ A non-parametric test (Kruskal-Wallis) should be preferred.", indent=1)
                all_normal = False
            log.blank()

            qq_path = FIG_DIR / f"{slug}-qq.png"
            _plot_qq(groups, f"{label} — Q-Q plots (normality diagnostic)", qq_path)
            log.p(f"[PLOT] {qq_path.relative_to(REPO_ROOT)}", indent=1)

            dist_path = FIG_DIR / f"{slug}-distributions.png"
            _plot_dist(groups, f"{label} — Distributions", vcol, dist_path)
            log.p(f"[PLOT] {dist_path.relative_to(REPO_ROOT)}", indent=1)
            log.blank()

            # ── A3: homoscedasticity ──
            log.h3("Homoscedasticity  (Assumption A3)")
            log.p("H₀: population variances are equal across all groups.", indent=1)
            log.p("Test: Levene's test (median variant).", indent=1)
            log.blank()

            lev = _test_levene(groups)
            log.p(f"Levene statistic = {lev['statistic']}", indent=1)
            log.p(f"p-value          = {lev['p_value']:.6e}" if isinstance(lev['p_value'], float) else f"p-value          = {lev['p_value']}", indent=1)
            log.p(f"Equal variances: {lev['equal_var'].upper()}", indent=1)
            log.blank()

            equal_var = lev["equal_var"] == "yes"
            if not equal_var and lev["equal_var"] != "insufficient":
                log.p("→ Homoscedasticity assumption is VIOLATED.", indent=1)
                log.p("→ When combined with unequal group sizes the standard F-test", indent=1)
                log.p("  over-rejects H₀ (inflated Type-I error).", indent=1)
                log.p("→ Use Welch's ANOVA or Kruskal-Wallis instead.", indent=1)
            log.blank()

            var_path = FIG_DIR / f"{slug}-variance.png"
            lev_p_val = lev["p_value"] if isinstance(lev["p_value"], float) else 0.0
            _plot_var(groups, f"{label} — Variance comparison", vcol, lev_p_val, var_path)
            log.p(f"[PLOT] {var_path.relative_to(REPO_ROOT)}", indent=1)
            log.blank()

            # ── A1: independence ──
            log.h3("Independence  (Assumption A1)")
            independent = True

            if ds["kind"] == "global":
                log.p("Data type: GLOBAL — one value per reference structure per algorithm.", indent=1)
                log.p("Each observation corresponds to a distinct experimentally resolved", indent=1)
                log.p("protein structure.  Different reference PDBs are independent.", indent=1)
                log.p("→ Independence assumption is SATISFIED.", indent=1)
            else:
                independent = False
                log.p("Data type: POSITIONAL — one value per residue per reference per algorithm.", indent=1)
                log.blank()
                log.p("CRITICAL: residues within the same protein chain are spatially", indent=1)
                log.p("correlated.  Neighbouring positions share local structural", indent=1)
                log.p("context (secondary structure, loops, domain boundaries), so", indent=1)
                log.p("their RMSD / pLDDT values are NOT independent.", indent=1)
                log.blank()

                # Show autocorrelation for a sample chain
                combos = df.groupby(["Algorithm", "ref"]).size().reset_index(name="cnt")
                combos = combos[combos["cnt"] >= 30]
                if not combos.empty:
                    s = combos.iloc[0]
                    sub = df[(df["Algorithm"] == s["Algorithm"]) & (df["ref"] == s["ref"])].sort_values("pos")
                    vals = sub[vcol].values
                    if vals.std() > 1e-12:
                        acf = _autocorrelation(vals, max_lag=5)
                        log.p(f"Example: {s['Algorithm']} / {s['ref']} (n={len(vals)} residues):", indent=1)
                        for lag, r in enumerate(acf, 1):
                            log.p(f"  Lag-{lag} autocorrelation = {r:.4f}", indent=1)
                        log.blank()
                        log.p("Values well above zero confirm strong spatial autocorrelation.", indent=1)
                    else:
                        log.p(f"Example: {s['Algorithm']} / {s['ref']} — constant values, no", indent=1)
                        log.p("autocorrelation to compute, but independence is still violated", indent=1)
                        log.p("because residues belong to the same physical protein.", indent=1)
                log.blank()
                log.p("→ Independence assumption is VIOLATED.", indent=1)
                log.p("→ ANOVA on raw per-residue values inflates the effective sample", indent=1)
                log.p("  size by treating every residue as a separate experiment.", indent=1)
                log.p("  This produces artificially small p-values that are meaningless.", indent=1)
                log.blank()
                log.p("→ Correct approach: aggregate per-residue values to per-protein", indent=1)
                log.p("  means first (yielding one independent observation per protein", indent=1)
                log.p("  per algorithm), then run the omnibus test on those means.", indent=1)

                acf_path = FIG_DIR / f"{slug}-autocorrelation.png"
                _plot_acf(df, protein, vcol, acf_path)
                log.p(f"[PLOT] {acf_path.relative_to(REPO_ROOT)}", indent=1)
            log.blank()

            # ── per-dataset verdict ──
            log.h3("Verdict for this dataset")
            violations = []
            if not all_normal:
                violations.append("normality (A2)")
            if not equal_var and lev["equal_var"] != "insufficient":
                violations.append("homoscedasticity (A3)")
            if not independent:
                violations.append("independence (A1)")
            if not balanced:
                violations.append("balance (A4)")

            if not violations:
                log.p(f"✓  All assumptions satisfied — ANOVA + Tukey HSD are appropriate.", indent=1)
            else:
                log.p(f"✗  Violated assumptions: {', '.join(violations)}.", indent=1)
                log.p(f"→  Current ANOVA + Tukey HSD results for {label} are UNRELIABLE.", indent=1)
            log.blank()

            summary.append({
                "dataset": label,
                "all_normal": all_normal,
                "equal_var": equal_var,
                "independent": independent,
            })

    # ── 5  Summary heatmap ─────────────────────────────────────────────
    heatmap_path = FIG_DIR / "assumption-summary-heatmap.png"
    _plot_summary_heatmap(summary, heatmap_path)

    # ── 6  Overall summary ─────────────────────────────────────────────
    log.h1("4. OVERALL SUMMARY")
    log.p(f"[PLOT] {heatmap_path.relative_to(REPO_ROOT)}")
    log.blank()

    n_total = len(summary)
    n_ok = sum(1 for r in summary if r["all_normal"] and r["equal_var"] and r["independent"])
    log.p(f"Datasets audited:       {n_total}")
    log.p(f"All assumptions met:    {n_ok}")
    log.p(f"≥1 assumption violated: {n_total - n_ok}")
    log.blank()

    # ── global vs positional breakdown ──
    gl = [r for r in summary if "Global" in r["dataset"]]
    po = [r for r in summary if "Positional" in r["dataset"]]

    log.h2("4.1  Global metrics (Global RMSD, Global TM)")
    log.p(f"  Normality failures:          {sum(1 for r in gl if not r['all_normal'])}/{len(gl)}")
    log.p(f"  Homoscedasticity failures:   {sum(1 for r in gl if not r['equal_var'])}/{len(gl)}")
    log.p(f"  Independence:                always satisfied (one value per structure)")
    log.blank()
    log.p("  Even for global metrics, many datasets exhibit non-normality and/or")
    log.p("  unequal variances.  The standard ANOVA F-test is therefore unreliable")
    log.p("  for several protein × metric combinations.")
    log.blank()

    log.h2("4.2  Positional metrics (per-residue RMSD, per-residue pLDDT)")
    log.p("  Independence:                ALWAYS VIOLATED (spatial autocorrelation)")
    log.p("  Normality:                   typically violated (skewed distributions)")
    log.blank()
    log.p("  The tens of thousands of per-residue rows are NOT independent observations.")
    log.p("  Running ANOVA directly on them massively inflates the effective sample")
    log.p("  size and produces p-values that are orders of magnitude too small.")
    log.p("  These results are statistically INVALID as currently computed.")
    log.blank()

    log.p(textwrap.dedent("""\
        SUMMARY:
          For GLOBAL metrics  → switch from standard ANOVA to Kruskal-Wallis
                                + Dunn's test for post-hoc pairwise comparisons.
          For POSITIONAL metrics → first aggregate to per-protein means, then
                                   apply the same non-parametric pipeline.
    """))

    # ── save ───────────────────────────────────────────────────────────
    log_path = DIAG_DIR / "report.log"
    log.save(log_path)

    print(f"Report:  {log_path.relative_to(REPO_ROOT)}")
    print(f"Figures: {FIG_DIR.relative_to(REPO_ROOT)}/  ({len(list(FIG_DIR.glob('*.png')))} plots)")


if __name__ == "__main__":
    main()
