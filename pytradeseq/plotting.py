"""Visualisation — 1:1 port of tradeSeq::plotSmoothers using ggplot2-python.

R signature::

    plotSmoothers(models, counts, gene, nPoints = 100, lwd = 2,
                  size = 2/3, xlab = "Pseudotime",
                  ylab = "Log(expression + 1)",
                  border = FALSE, alpha = 1, sample = 1)

Returns a ggplot2-python plot. Render with::

    from ggplot2_py import ggsave
    ggsave("smoother.png", plot=p, width=6, height=4, dpi=120)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    geom_line,
    geom_point,
    ggplot,
    labs,
    scale_color_viridis_d,
    theme_classic,
)


def plotSmoothers(
    gam,
    counts: np.ndarray | pd.DataFrame,
    gene,
    *,
    nPoints: int = 100,
    lwd: float = 2.0,
    size: float = 2.0 / 3.0,
    xlab: str = "Pseudotime",
    ylab: str = "Log(expression + 1)",
    border: bool = False,
    alpha: float = 1.0,
    sample: float = 1.0,
):
    """1:1 port of tradeSeq::plotSmoothers.

    Args:
        gam: a fitted ``TradeSeqGAMs`` object.
        counts: full (n_cells × n_genes) count matrix.
        gene: index or name of the gene.
        nPoints: number of grid points for smoother lines.
        lwd: smoother line width.
        size: scatter point size.
        border: draw a white halo around the smoother (R parity).
        alpha: viridis discrete-palette alpha.
        sample: fraction of cells to subsample for the scatter (1.0 = all).
    """
    # Resolve gene → integer index in the gam.models list
    if isinstance(gene, str):
        gene_idx = gam.gene_names.index(gene)
    else:
        gene_idx = int(gene)

    if isinstance(counts, pd.DataFrame):
        # counts may be cells × genes or genes × cells; pick row/col matching gene name
        if gene in counts.columns:
            y = counts[gene].to_numpy(dtype=np.float64)
        elif gene in counts.index:
            y = counts.loc[gene].to_numpy(dtype=np.float64)
        else:
            y = counts.iloc[:, gene_idx].to_numpy(dtype=np.float64)
    else:
        y = np.asarray(counts, dtype=np.float64)[:, gene_idx]

    # Pseudotime + lineage assignment per cell
    pseudotime = np.asarray(gam.pseudotime, dtype=np.float64)
    cell_w = np.asarray(gam.cell_assignments, dtype=np.float64)
    n_cells, n_lineages = cell_w.shape

    time_all = np.zeros(n_cells)
    lineage_all = np.zeros(n_cells, dtype=int)
    for jj in range(n_lineages):
        mask = cell_w[:, jj] == 1
        time_all[mask] = pseudotime[mask, jj] if pseudotime.ndim == 2 else pseudotime[mask]
        lineage_all[mask] = jj + 1

    df = pd.DataFrame(
        {
            "time": time_all,
            "gene_count": y,
            "lineage": [str(int(c)) for c in lineage_all],
        }
    )
    if sample < 1.0:
        df = df.sample(frac=sample, random_state=0)

    df = df.assign(log_expr=np.log1p(df["gene_count"].to_numpy()))

    p = (
        ggplot(df, aes(x="time", y="log_expr", colour="lineage"))
        + geom_point(size=size)
        + labs(x=xlab, y=ylab)
        + theme_classic()
        + scale_color_viridis_d(alpha=alpha)
    )

    # Smoother curves per lineage
    from .tests_module import _design_row_at_t
    for jj in range(n_lineages):
        m = gam.models[gene_idx][jj]
        if m is None or not gam.converged[gene_idx, jj]:
            continue
        lin_pt = pseudotime[:, jj] if pseudotime.ndim == 2 else pseudotime
        finite = np.isfinite(lin_pt)
        if finite.sum() == 0:
            continue
        # Clip grid to the knot range so the spline basis is defined
        knots = gam.knot_positions[jj]
        t_lo = max(float(knots.min()), float(lin_pt[finite].min()))
        t_hi = min(float(knots.max()), float(lin_pt[finite].max()))
        t_grid = np.linspace(t_lo, t_hi, nPoints)
        eta = np.array([_design_row_at_t(m, t) @ np.asarray(m.params) for t in t_grid])
        # Add the mean offset (log library size) for cells assigned to this lineage,
        # matching tradeSeq::predict(..., type="response") on an "average" cell.
        cell_mask = (cell_w[:, jj] == 1)
        mean_off = float(np.mean(gam.offsets[cell_mask])) if cell_mask.any() else 0.0
        yhat = np.exp(eta + mean_off)
        smooth_df = pd.DataFrame(
            {
                "time": t_grid,
                "log_expr": np.log1p(yhat),
                "lineage": str(jj + 1),
            }
        )
        if border:
            p = p + geom_line(
                aes(x="time", y="log_expr"),
                data=smooth_df.assign(lineage=str(jj + 1)),
                size=lwd + 1,
                colour="white",
            )
        p = p + geom_line(
            aes(x="time", y="log_expr", colour="lineage"),
            data=smooth_df,
            size=lwd,
        )

    return p
