from typing import Optional, Tuple

import numpy as np


def linear_regression_fit(
    x: np.ndarray, y: np.ndarray, n_points: int = 200
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Compute a simple linear regression fit over finite values.

    Parameters
    ----------
    x, y : array-like
        Observed data.  Non-finite entries (NaN / Inf) are silently dropped.
    n_points : int
        Number of evenly spaced points for the fitted line.

    Returns
    -------
    (x_line, y_line, coeffs) or ``None``
        *x_line* and *y_line* are the coordinates of the fit line, *coeffs*
        is ``[slope, intercept]``.  Returns ``None`` when fewer than two
        finite observations remain.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return None
    coeffs = np.polyfit(x[mask], y[mask], 1)
    x_line = np.linspace(x[mask].min(), x[mask].max(), n_points)
    y_line = np.polyval(coeffs, x_line)
    return x_line, y_line, coeffs


def bin_means(
    x, y, n_bins: int = 4
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Compute per-bin means over equal-width bins along *x*.

    Finite (x, y) pairs are divided into *n_bins* equal-width bins spanning
    ``[x.min(), x.max()]``.  For each non-empty bin the mean x, mean y, and
    observation count are returned.

    Parameters
    ----------
    x, y : array-like
        Observed data.  Non-finite entries are silently dropped.
    n_bins : int
        Number of equal-width bins along *x*.

    Returns
    -------
    (mean_x, mean_y, weights) or ``None``
        Arrays of length equal to the number of non-empty bins.
        *weights* contains the count of observations per bin.
        Returns ``None`` when fewer than two non-empty bins exist.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return None

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
        return None

    return np.asarray(mean_x), np.asarray(mean_y), np.asarray(weights, dtype=float)


def binned_mean_regression_fit(
    x: np.ndarray, y: np.ndarray, n_bins: int = 4, n_points: int = 200
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Linear regression through per-bin mean points.

    The finite (x, y) observations are divided into *n_bins* equal-width bins
    along *x*.  For each bin the mean x and mean y are computed, yielding
    *n_bins* representative points.  A standard linear regression
    (via :func:`linear_regression_fit`) is then fitted through those means.

    Parameters
    ----------
    x, y : array-like
        Observed data.  Non-finite entries are silently dropped.
    n_bins : int
        Number of equal-width bins along *x*.
    n_points : int
        Number of evenly spaced points for the returned fit line.

    Returns
    -------
    (x_line, y_line, coeffs, mean_x, mean_y) or ``None``
        First three elements follow the same contract as
        :func:`linear_regression_fit`.  *mean_x* and *mean_y* are the
        per-bin mean coordinates used for the fit.  Returns ``None``
        when fewer than two non-empty bins remain.
    """
    result = bin_means(x, y, n_bins)
    if result is None:
        return None
    mean_x, mean_y, _ = result
    fit = linear_regression_fit(mean_x, mean_y, n_points)
    if fit is None:
        return None
    x_line, y_line, coeffs = fit
    return x_line, y_line, coeffs, mean_x, mean_y


def weighted_binned_correlation(x, y, n_bins: int = 4):
    """Weighted Pearson correlation over per-bin means.

    Observations are binned into *n_bins* equal-width bins along *x*
    (via :func:`bin_means`).  A weighted Pearson *r* is then computed
    over the bin means, using the bin counts as weights.

    Returns ``NaN`` when fewer than two non-empty bins exist or when
    either dimension has zero variance across bins.
    """
    result = bin_means(x, y, n_bins)
    if result is None:
        return np.nan
    mean_x, mean_y, weights = result

    mx = np.average(mean_x, weights=weights)
    my = np.average(mean_y, weights=weights)
    cov = np.average((mean_x - mx) * (mean_y - my), weights=weights)
    vx = np.average((mean_x - mx) ** 2, weights=weights)
    vy = np.average((mean_y - my) ** 2, weights=weights)
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)


def count_overlapping(seq, motif):
    m = len(motif)
    count = 0
    for i in range(len(seq) - m + 1):
        if seq[i:i + m] == motif:
            count += 1
    return count


def _to_alignment_string(alignment) -> str:
    """Normalize supported alignment inputs to a plain sequence string."""
    if hasattr(alignment, "seq"):
        return str(alignment.seq)
    return str(alignment)


def pretty_print_alignments(ref_align, pred_align, line_width: int = 80) -> None:
    """Print two alignments in visually aligned blocks for core-logic debugging.

    This is a development helper intended for inspecting how reference and
    predicted alignment strings compare after each transformation step in
    `core.py`. Inputs can be plain strings or objects with a `.seq` attribute
    (e.g., BioPython `SeqRecord`).
    """
    if line_width <= 0:
        raise ValueError("line_width must be a positive integer")

    ref = _to_alignment_string(ref_align)
    pred = _to_alignment_string(pred_align)

    # Pad the shorter sequence so all chunks keep their columns aligned.
    width = max(len(ref), len(pred))
    ref = ref.ljust(width)
    pred = pred.ljust(width)

    for start in range(0, width, line_width):
        end = min(start + line_width, width)
        ref_chunk = ref[start:end]
        pred_chunk = pred[start:end]
        match_chunk = "".join("|" if r == p else " " for r, p in zip(ref_chunk, pred_chunk))

        print(f"[{start:>4}:{end:<4}]")
        print(f"REF  {ref_chunk}")
        print(f"     {match_chunk}")
        print(f"PRED {pred_chunk}")
        print()
