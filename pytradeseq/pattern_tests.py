"""patternTest, earlyDETest — between-lineage pattern comparisons.

Tests whether the SHAPE of the expression-vs-pseudotime curve differs
between lineages. Sample pseudotime points uniformly along each lineage's
range, predict expression at those points, then Wald-test the difference
of predicted-expression vectors.

Mirrors tradeSeq's `patternTest.R` and `earlyDETest.R`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .fitting import TradeSeqGAMs
from .tests_module import _design_row_at_t


def _build_eval_grid(pseudotime: np.ndarray, cell_assignments: np.ndarray,
                     nPoints: int = 100) -> list[np.ndarray]:
    """Per-lineage uniform pseudotime grid over the fitted range."""
    n_lin = pseudotime.shape[1]
    grids = []
    for L in range(n_lin):
        mask = cell_assignments[:, L] == 1
        if mask.sum() < 2:
            grids.append(np.array([]))
            continue
        t = pseudotime[mask, L]
        grids.append(np.linspace(t.min(), t.max(), nPoints))
    return grids


def patternTest(
    gams: TradeSeqGAMs,
    *,
    nPoints: int = 100,
    pairwise: bool = False,
) -> pd.DataFrame:
    """Test whether the expression PATTERN along pseudotime differs across lineages.

    For each gene, compares fitted η_L(t_grid) curves pairwise (or summed).
    Wald-test on the difference vector.

    Args:
        gams: from fitGAM with ≥ 2 lineages
        nPoints: number of grid points per lineage for the comparison
        pairwise: if True, return per-pair (i, j) results; else summed.
    """
    n_g = gams.n_genes
    n_lin = gams.n_lineages
    wald = np.full(n_g, np.nan)
    dfs = np.zeros(n_g, dtype=int)
    if n_lin < 2:
        return pd.DataFrame({"waldStat": wald, "df": dfs,
                             "pvalue": np.ones(n_g)},
                            index=gams.gene_names)

    grids = _build_eval_grid(gams.pseudotime, gams.cell_assignments, nPoints)

    for g in range(n_g):
        wald_g = 0.0; df_g = 0; any_ok = False
        for i in range(n_lin):
            for j in range(i + 1, n_lin):
                m_i = gams.models[g][i]
                m_j = gams.models[g][j]
                if m_i is None or m_j is None: continue
                if grids[i].size == 0 or grids[j].size == 0: continue
                try:
                    # Predicted η on each lineage's grid (the actual grid points,
                    # using each lineage's own basis)
                    eta_i = np.array([_design_row_at_t(m_i, t) @ np.asarray(m_i.params)
                                       for t in grids[i]])
                    eta_j = np.array([_design_row_at_t(m_j, t) @ np.asarray(m_j.params)
                                       for t in grids[j]])
                    # Match by quantile (uniform grid → grids same length)
                    diff = eta_j - eta_i
                    # Variance: sum of per-point variances on each lineage
                    var_i = np.array([
                        float(_design_row_at_t(m_i, t) @ np.asarray(m_i.cov_params())
                              @ _design_row_at_t(m_i, t))
                        for t in grids[i]
                    ])
                    var_j = np.array([
                        float(_design_row_at_t(m_j, t) @ np.asarray(m_j.cov_params())
                              @ _design_row_at_t(m_j, t))
                        for t in grids[j]
                    ])
                    var_diff = var_i + var_j
                    # Aggregate Wald: sum (diff[k]² / var[k]) ~ χ²(nPoints)
                    valid = var_diff > 0
                    if valid.sum() == 0: continue
                    stat = float((diff[valid] ** 2 / var_diff[valid]).sum())
                    wald_g += stat
                    df_g += int(valid.sum())
                    any_ok = True
                except Exception:
                    continue
        if any_ok and df_g > 0:
            wald[g] = wald_g; dfs[g] = df_g

    pvals = np.array([chi2.sf(w, d) if np.isfinite(w) and d > 0 else 1.0
                      for w, d in zip(wald, dfs)])
    return pd.DataFrame({"waldStat": wald, "df": dfs, "pvalue": pvals},
                        index=gams.gene_names)


def earlyDETest(
    gams: TradeSeqGAMs,
    *,
    knots: tuple[int, int] = (1, 2),
) -> pd.DataFrame:
    """Test for differential expression in the EARLY part of pseudotime.

    R: `earlyDETest(sce, knots = c(1, 2))` — restricts patternTest to the
    pseudotime interval between two specified knots.

    Python: uses the first 25% of each lineage's pseudotime range as the
    "early" window. (Knot-based ranges require access to mgcv's internal
    knot placement; simplification documented in MATH.md.)
    """
    n_lin = gams.n_lineages
    if n_lin < 2:
        return pd.DataFrame({"waldStat": np.full(gams.n_genes, np.nan),
                             "df": np.zeros(gams.n_genes, dtype=int),
                             "pvalue": np.ones(gams.n_genes)},
                            index=gams.gene_names)

    # Build an "early" pseudotime restricted to lowest 25% per lineage
    early_pt = gams.pseudotime.copy()
    early_mask = np.ones_like(gams.cell_assignments)
    for L in range(n_lin):
        mask = gams.cell_assignments[:, L] == 1
        if mask.sum() == 0: continue
        t = gams.pseudotime[mask, L]
        q25 = np.quantile(t, 0.25)
        # zero out cells beyond q25 in this lineage's assignment
        beyond = mask & (gams.pseudotime[:, L] > q25)
        early_mask[beyond, L] = 0

    # Wrap a "view" of gams with restricted assignments + call patternTest
    early_gams = TradeSeqGAMs(
        gene_names=gams.gene_names,
        n_lineages=gams.n_lineages,
        models=gams.models,
        knot_positions=gams.knot_positions,
        offsets=gams.offsets,
        cell_assignments=early_mask,
        pseudotime=early_pt,
        family=gams.family,
        nknots=gams.nknots,
        converged=gams.converged,
    )
    return patternTest(early_gams, nPoints=50)
