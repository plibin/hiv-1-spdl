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
