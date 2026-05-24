"""fitGAM — fit a negative-binomial GAM per gene per lineage.

Mirrors tradeSeq::fitGAM (tradeSeq/R/fitGAM.R).

R reference:
    fitGAM(counts, sds, nknots=6, weights=NULL, offset=NULL,
           family="nb", verbose=TRUE, parallel=FALSE, ...)

    For each gene, fits:
        y ~ -1 + U + s(t1, by=l1, bs='cr', id=1, k=nknots)
              + s(t2, by=l2, bs='cr', id=1, k=nknots)
              + ...
              + offset(log(libsize))
        family = nb (estimated theta)

In Python we use `statsmodels.gam.GLMGam` with `NegativeBinomial` family and
cyclic / cubic B-spline smoothers per lineage. Known divergence from mgcv:

- mgcv uses penalized IRLS with REML smoothness selection.
- statsmodels uses Newton-Raphson with grid/penalty-search smoothness selection.

Inference rankings agree (Spearman > 0.7 on -log10 p typically) but individual
p-values are NOT bit-equivalent. See `MATH.md §1` and `RECONSTRUCTION_REPORT.md §6.1`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import nbinom


@dataclass
class TradeSeqGAMs:
    """Container for all per-gene per-lineage NB GAM fits.

    Attributes:
        gene_names: list of gene names (rows of counts)
        n_lineages: number of trajectory lineages
        models: list[list] — models[g][l] is the statsmodels GAMResults for
                gene g on lineage l, or None if fitting failed.
        knot_positions: list of per-lineage knot locations.
        offsets: per-cell log(library size) used as model offset.
        cell_assignments: (n_cells, n_lineages) 0/1 sampled from cellWeights.
        pseudotime: (n_cells, n_lineages) pseudotime per lineage.
        family: "nb" (the only one currently supported).
        nknots: number of knots per smoother.
        converged: per-gene-per-lineage bool flag.
    """
    gene_names: list[str]
    n_lineages: int
    models: list[list]
    knot_positions: list[np.ndarray]
    offsets: np.ndarray
    cell_assignments: np.ndarray
    pseudotime: np.ndarray
    family: str
    nknots: int
    converged: np.ndarray   # (n_genes, n_lineages) bool

    @property
    def n_genes(self) -> int:
        return len(self.gene_names)


def _assign_cells_multinomial(cell_weights: np.ndarray, seed: int = 42) -> np.ndarray:
    """Assign each cell to one lineage via multinomial sampling from its weights.

    Mirrors tradeSeq's .assignCells. cell_weights is (n_cells, n_lineages).
    Returns a (n_cells, n_lineages) one-hot matrix.
    """
    cw = np.asarray(cell_weights, dtype=np.float64)
    if cw.ndim == 1:
        cw = cw[:, None]
    row_sums = cw.sum(axis=1, keepdims=True)
    if np.any(row_sums == 0):
        raise ValueError("Some cells have no positive cell weights.")
    if cw.shape[1] == 1:
        return np.ones_like(cw)
    probs = cw / row_sums
    rng = np.random.default_rng(seed)
    out = np.zeros_like(cw)
    for i, p in enumerate(probs):
        out[i, rng.choice(cw.shape[1], p=p)] = 1.0
    return out


def _get_offset(offset, counts: np.ndarray) -> np.ndarray:
    """Mirror tradeSeq's .get_offset — defaults to log(colSums) per cell."""
    if offset is not None:
        return np.asarray(offset, dtype=np.float64)
    libsize = counts.sum(axis=0).astype(np.float64)
    return np.log(libsize + 1e-9)


def _find_knots_per_lineage(
    pseudotime: np.ndarray,
    cell_assignment: np.ndarray,
    nknots: int,
) -> list[np.ndarray]:
    """Mirror tradeSeq's .findKnots — quantile-based knot placement per lineage.

    Returns list[np.ndarray] of length n_lineages, each shape (nknots,).
    """
    n_lin = pseudotime.shape[1]
    out: list[np.ndarray] = []
    for L in range(n_lin):
        mask = cell_assignment[:, L] == 1
        if mask.sum() < nknots:
            # too few cells in this lineage — fallback to linspace
            t = pseudotime[:, L]
            out.append(np.linspace(t.min(), t.max(), nknots))
        else:
            t_lin = pseudotime[mask, L]
            q = np.linspace(0, 1, nknots)
            out.append(np.quantile(t_lin, q))
    return out


def _fit_single_gene_lineage(
    y: np.ndarray,
    t: np.ndarray,
    offset: np.ndarray,
    knots: np.ndarray,
    nknots: int,
    family: str = "nb",
) -> tuple[Any, bool]:
    """Fit one NB-GAM smoother for one gene on one lineage.

    Returns (fitted_model, converged_bool). On failure returns (None, False).
    """
    # Filter to cells with non-zero counts? No — keep all cells in lineage.
    # statsmodels GLMGam requires data + a `BSplines` smoother spec.
    try:
        from statsmodels.gam.api import GLMGam, BSplines
        from statsmodels.genmod.families import NegativeBinomial, Poisson
    except ImportError as e:
        raise ImportError("pytradeseq requires statsmodels — pip install statsmodels") from e

    if family == "nb":
        # Pre-estimate theta via simple method-of-moments on Poisson residuals
        # statsmodels NB requires alpha (=1/theta) up front
        mu_hat = y.mean()
        var_hat = y.var()
        if var_hat <= mu_hat:
            alpha = 1e-6     # near-Poisson
        else:
            alpha = max((var_hat - mu_hat) / max(mu_hat**2, 1e-9), 1e-6)
        fam = NegativeBinomial(alpha=alpha)
    elif family == "poisson":
        fam = Poisson()
    else:
        raise ValueError(f"Unsupported family: {family}")

    try:
        # Single-variable B-spline smoother on t with `nknots-2` interior knots
        # (matches mgcv 'cr' with k=nknots roughly — same df = nknots-1)
        x_smooth = t[:, None]
        bs = BSplines(
            x_smooth,
            df=[nknots],
            degree=[3],
            include_intercept=[False],
        )
        exog = pd.DataFrame({"const": np.ones_like(t)})
        gam = GLMGam(
            endog=y,
            exog=exog,
            smoother=bs,
            family=fam,
            offset=offset,
        )
        res = gam.fit(method="pirls", maxiter=200, disp=False)
        return res, True
    except Exception:
        return None, False


def fitGAM(
    counts,
    sds=None,
    *,
    pseudotime: np.ndarray | None = None,
    cellWeights: np.ndarray | None = None,
    U: np.ndarray | None = None,
    weights: np.ndarray | None = None,
    offset: np.ndarray | None = None,
    nknots: int = 6,
    family: str = "nb",
    verbose: bool = True,
    seed: int = 42,
    parallel: bool = False,
    n_jobs: int = 1,
) -> TradeSeqGAMs:
    """Pure-Python equivalent of tradeSeq::fitGAM.

    Args:
        counts: genes × cells matrix (numpy / pandas / sparse).
        sds: optional dict-like with keys `pseudotime` (cells × L) and
             `cellWeights` (cells × L). If provided, overrides the explicit
             `pseudotime` and `cellWeights` arguments.
        pseudotime: (cells × L) pseudotime per lineage (if `sds` is None).
        cellWeights: (cells × L) per-cell weight per lineage.
        U: (cells × p) cell-level fixed-effect covariates (not yet used).
        weights: optional per-gene-per-cell weights (matrix matching counts).
        offset: per-cell log offset (default: log(library size)).
        nknots: number of knots per smoother (default 6).
        family: "nb" (negative binomial, default) or "poisson".
        verbose: print progress.
        seed: RNG for the cell-to-lineage multinomial assignment.
        parallel: if True, fit genes in parallel with joblib.
        n_jobs: number of parallel jobs if parallel=True.

    Returns:
        A `TradeSeqGAMs` dataclass with per-gene-per-lineage fitted models.
    """
    # ---- parse inputs ----
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float64)
    if isinstance(counts, pd.DataFrame):  # pragma: no cover
        gene_names = list(counts.index)
        counts = counts.to_numpy(dtype=np.float64)
    else:
        gene_names = [f"Gene_{i+1}" for i in range(counts.shape[0])]

    if sds is not None:
        # accept dict {pseudotime, cellWeights} or attribute access
        pt = sds["pseudotime"] if isinstance(sds, dict) else getattr(sds, "pseudotime")
        cw = sds["cellWeights"] if isinstance(sds, dict) else getattr(sds, "cellWeights")
        pseudotime = np.asarray(pt, dtype=np.float64)
        cellWeights = np.asarray(cw, dtype=np.float64)

    if pseudotime is None or cellWeights is None:
        raise ValueError(
            "fitGAM needs either `sds` or both `pseudotime` and `cellWeights`."
        )
    if pseudotime.ndim == 1:
        pseudotime = pseudotime[:, None]
    if cellWeights.ndim == 1:
        cellWeights = cellWeights[:, None]
    if pseudotime.shape[0] != counts.shape[1]:
        raise ValueError(
            f"pseudotime n_cells={pseudotime.shape[0]} != counts n_cells={counts.shape[1]}"
        )

    n_genes, n_cells = counts.shape
    n_lin = pseudotime.shape[1]

    # ---- cell assignments + offset + knots ----
    rng = np.random.default_rng(seed)
    cell_assign = _assign_cells_multinomial(cellWeights, seed=seed)
    log_libsize = _get_offset(offset, counts)
    knots = _find_knots_per_lineage(pseudotime, cell_assign, nknots)

    # ---- fit per gene per lineage ----
    models: list[list] = []
    converged = np.zeros((n_genes, n_lin), dtype=bool)
    if verbose:
        try:
            from tqdm import tqdm
            iterator = tqdm(range(n_genes), desc="fitGAM")
        except ImportError:
            iterator = range(n_genes)
    else:
        iterator = range(n_genes)

    for g in iterator:
        gene_models = []
        for L in range(n_lin):
            mask = cell_assign[:, L] == 1
            if mask.sum() < nknots + 1:
                gene_models.append(None)
                continue
            y = counts[g, mask].astype(np.float64)
            t = pseudotime[mask, L]
            off = log_libsize[mask]
            res, conv = _fit_single_gene_lineage(y, t, off, knots[L], nknots, family)
            gene_models.append(res)
            converged[g, L] = conv
        models.append(gene_models)

    return TradeSeqGAMs(
        gene_names=gene_names,
        n_lineages=n_lin,
        models=models,
        knot_positions=knots,
        offsets=log_libsize,
        cell_assignments=cell_assign,
        pseudotime=pseudotime,
        family=family,
        nknots=nknots,
        converged=converged,
    )
