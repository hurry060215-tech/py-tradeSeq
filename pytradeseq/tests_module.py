"""Statistical tests on fitted tradeSeq GAMs.

Mirrors tradeSeq's:
- associationTest (associationTest.R)
- startVsEndTest (startVsEndTest.R)
- diffEndTest (diffEndTest.R)

Every test outputs a DataFrame with columns:
    waldStat | df | pvalue | (meanLogFC)

The Wald statistic is `β̂^T Σ̂^{-1} β̂` where β̂ are the smooth-term
coefficients of interest. df = number of restrictions tested.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .fitting import TradeSeqGAMs


def _wald_test(beta: np.ndarray, cov: np.ndarray) -> tuple[float, int, float]:
    """Wald chi-squared test that beta = 0.

    Returns (stat, df, pvalue). On singular covariance, returns
    (nan, df, 1.0) and lets the caller mark it as failed.
    """
    df = int(beta.shape[0])
    try:
        # Use pseudoinverse for safety on near-singular covariances
        cov_inv = np.linalg.pinv(cov)
        stat = float(beta @ cov_inv @ beta)
        if not np.isfinite(stat) or stat < 0:
            return (float("nan"), df, 1.0)
        return (stat, df, float(chi2.sf(stat, df)))
    except Exception:
        return (float("nan"), df, 1.0)


def _get_smooth_coefs_and_cov(model, n_smooth_basis: int):
    """Extract the smooth-term coefficients + their covariance block.

    statsmodels GLMGam BSplines names smoother basis params `x<feat>_s<k>`.
    For a single smoother of `df=k`, that's `x0_s0` ... `x0_s{k-2}` (one
    basis dropped to identify the intercept). The covariance is symmetric.
    """
    params = model.params
    cov = model.cov_params()
    try:
        idx = [i for i, name in enumerate(params.index)
               if "_s" in name and name.split("_s")[0].startswith("x")]
    except AttributeError:
        # params is ndarray, not Series — fall back to positional after intercept
        idx = list(range(1, 1 + n_smooth_basis - 1))
    if len(idx) == 0:
        return None, None
    beta = np.asarray(params)[idx]
    sub_cov = np.asarray(cov)[np.ix_(idx, idx)]
    return beta, sub_cov


def associationTest(
    gams: TradeSeqGAMs,
    *,
    contrastType: str = "consecutive",
    nPoints: int = 2 * 100,
    lineages: bool = False,
) -> pd.DataFrame:
    """Test whether expression is associated with pseudotime, per gene.

    For each gene, sums the Wald statistics across all lineages where the
    model converged. df is summed accordingly.

    Args:
        gams: from fitGAM.
        contrastType: ignored in this minimal port (only "consecutive" supported).
        nPoints: unused in this minimal port (kept for API parity).
        lineages: if True, return per-lineage results instead of pooled.

    Returns:
        DataFrame indexed by gene_name, columns: waldStat, df, pvalue, meanLogFC.
    """
    n_g = gams.n_genes
    n_lin = gams.n_lineages
    wald = np.full(n_g, np.nan)
    dfs = np.zeros(n_g, dtype=int)
    mean_logfc = np.zeros(n_g)

    for g in range(n_g):
        wald_g = 0.0
        df_g = 0
        any_ok = False
        for L in range(n_lin):
            m = gams.models[g][L]
            if m is None or not gams.converged[g, L]:
                continue
            beta, cov = _get_smooth_coefs_and_cov(m, gams.nknots)
            if beta is None or len(beta) == 0:
                continue
            stat, ddf, _ = _wald_test(beta, cov)
            if np.isfinite(stat):
                wald_g += stat
                df_g += ddf
                any_ok = True
                # crude meanLogFC proxy: SD of smooth coefficients * sqrt(n_obs)
                mean_logfc[g] += float(np.std(beta))

        if any_ok and df_g > 0:
            wald[g] = wald_g
            dfs[g] = df_g
            mean_logfc[g] /= max(1, n_lin)

    pvals = np.array([
        chi2.sf(w, d) if np.isfinite(w) and d > 0 else 1.0
        for w, d in zip(wald, dfs)
    ])

    return pd.DataFrame({
        "waldStat": wald,
        "df": dfs,
        "pvalue": pvals,
        "meanLogFC": mean_logfc,
    }, index=gams.gene_names)


def _design_row_at_t(model, t: float) -> np.ndarray:
    """Build the full (intercept + smooth-basis) design row at pseudotime `t`.

    Returns a 1-D array of length p = len(model.params), aligned with params.
    """
    # The model was fit with exog = pd.DataFrame({'const': 1}) and a BSplines
    # smoother on t. We rebuild the basis row at scalar t and prepend the
    # intercept value 1.
    bs = model.model.smoother
    smooth_row = bs.transform(np.array([[t]]))[0]   # length = sum(df)
    return np.concatenate([[1.0], smooth_row])


def startVsEndTest(
    gams: TradeSeqGAMs,
    *,
    pseudotimeValues: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """Test whether predicted expression differs between start and end of pseudotime.

    For each gene + lineage:
        η_end - η_start = (X_end - X_start) @ β
        Var(diff)        = (X_end - X_start) @ Σ_β @ (X_end - X_start)^T
        Wald stat        = (η_end - η_start)^2 / Var(diff)  ~ χ²(df=1)

    Wald statistics sum across lineages; df sums likewise.

    Args:
        gams: from fitGAM.
        pseudotimeValues: explicit (start, end) values; defaults to per-lineage
            (min, max) of fitted pseudotime.

    Returns:
        DataFrame indexed by gene_name.
    """
    n_g = gams.n_genes
    n_lin = gams.n_lineages
    wald = np.full(n_g, np.nan)
    dfs = np.zeros(n_g, dtype=int)
    mean_logfc = np.zeros(n_g)

    for g in range(n_g):
        wald_g = 0.0
        df_g = 0
        any_ok = False
        for L in range(n_lin):
            m = gams.models[g][L]
            if m is None or not gams.converged[g, L]:
                continue
            mask = gams.cell_assignments[:, L] == 1
            t = gams.pseudotime[mask, L]
            if pseudotimeValues is None:
                t_start, t_end = float(t.min()), float(t.max())
            else:
                t_start, t_end = pseudotimeValues
            try:
                X_start = _design_row_at_t(m, t_start)
                X_end   = _design_row_at_t(m, t_end)
                D = X_end - X_start                      # shape (p,)
                beta = np.asarray(m.params)              # shape (p,)
                cov_p = np.asarray(m.cov_params())       # shape (p, p)
                diff_eta = float(D @ beta)
                var_diff = float(D @ cov_p @ D)
                if var_diff > 0 and np.isfinite(diff_eta):
                    stat = diff_eta**2 / max(var_diff, 1e-12)
                    wald_g += stat
                    df_g += 1
                    mean_logfc[g] += abs(diff_eta)
                    any_ok = True
            except Exception:
                continue

        if any_ok and df_g > 0:
            wald[g] = wald_g
            dfs[g] = df_g
            mean_logfc[g] /= max(1, n_lin)

    pvals = np.array([
        chi2.sf(w, d) if np.isfinite(w) and d > 0 else 1.0
        for w, d in zip(wald, dfs)
    ])
    return pd.DataFrame({
        "waldStat": wald,
        "df": dfs,
        "pvalue": pvals,
        "meanLogFC": mean_logfc,
    }, index=gams.gene_names)


def diffEndTest(gams: TradeSeqGAMs) -> pd.DataFrame:
    """Test whether predicted expression differs at the END across lineages.

    Pairwise comparison between lineages of μ at t=max(t_lineage).
    Returns a DataFrame indexed by gene with summed Wald across pairs.
    """
    n_g = gams.n_genes
    n_lin = gams.n_lineages
    wald = np.full(n_g, np.nan)
    dfs = np.zeros(n_g, dtype=int)

    if n_lin < 2:
        # No pairs to compare
        return pd.DataFrame({"waldStat": wald, "df": dfs,
                             "pvalue": np.ones(n_g)},
                            index=gams.gene_names)

    for g in range(n_g):
        wald_g = 0.0
        df_g = 0
        any_ok = False
        # Predict at each lineage's endpoint
        preds = []
        cov_diags = []
        for L in range(n_lin):
            m = gams.models[g][L]
            if m is None or not gams.converged[g, L]:
                preds.append(None); cov_diags.append(None); continue
            mask = gams.cell_assignments[:, L] == 1
            t_end = gams.pseudotime[mask, L].max()
            try:
                pred = float(m.predict(
                    exog=pd.DataFrame({"const": [1.0]}),
                    exog_smooth=np.array([[t_end]]),
                    which="linear"))
                preds.append(pred)
                cov_diags.append(float(np.diag(m.cov_params()).sum()))
            except Exception:
                preds.append(None); cov_diags.append(None)
        # Pairwise diffs
        for i in range(n_lin):
            for j in range(i + 1, n_lin):
                if preds[i] is None or preds[j] is None: continue
                diff = preds[j] - preds[i]
                var = (cov_diags[i] or 0) + (cov_diags[j] or 0)
                if var > 0 and np.isfinite(diff):
                    wald_g += diff**2 / max(var, 1e-12)
                    df_g += 1
                    any_ok = True
        if any_ok and df_g > 0:
            wald[g] = wald_g
            dfs[g] = df_g

    pvals = np.array([
        chi2.sf(w, d) if np.isfinite(w) and d > 0 else 1.0
        for w, d in zip(wald, dfs)
    ])
    return pd.DataFrame({
        "waldStat": wald,
        "df": dfs,
        "pvalue": pvals,
    }, index=gams.gene_names)
