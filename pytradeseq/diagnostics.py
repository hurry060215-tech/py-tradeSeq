"""Diagnostic helpers: nknots, evaluateK.

Mirrors tradeSeq::nknots and tradeSeq::evaluateK.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .fitting import TradeSeqGAMs, fitGAM


def nknots(gams: TradeSeqGAMs) -> int:
    """Return the number of knots used to fit a TradeSeqGAMs object.

    Mirror R: `nknots(sce)` returns the integer nknots passed to fitGAM.
    """
    return int(gams.nknots)


def evaluateK(
    counts,
    *,
    sds=None,
    pseudotime=None,
    cellWeights=None,
    k_range: range = range(3, 12),
    nGenes: int = 50,
    seed: int = 42,
    verbose: bool = True,
) -> pd.DataFrame:
    """Scan candidate knot values; return per-k mean AIC for nGenes random genes.

    Mirror R: `evaluateK(counts, sds, k = 3:10)` — for each candidate k, fits
    GAMs to a random subsample of genes and reports goodness-of-fit metrics
    (AIC, deviance, R²). The user picks the smallest k where AIC plateaus.

    Args:
        counts: genes × cells matrix
        sds: SlingshotDataSet-like (dict with pseudotime, cellWeights) — OR
             provide both `pseudotime` and `cellWeights` arrays directly.
        k_range: candidate nknots values to try
        nGenes: number of genes to sample for the scan (default 50)
        seed: RNG

    Returns:
        DataFrame indexed by k, columns: mean_aic, mean_deviance, mean_loglik
    """
    rng = np.random.default_rng(seed)
    counts = np.asarray(counts, dtype=np.float64)
    n_genes_total = counts.shape[0]
    sample = rng.choice(n_genes_total, size=min(nGenes, n_genes_total), replace=False)
    counts_sub = counts[sample]

    results = []
    for k in k_range:
        gams = fitGAM(
            counts_sub, sds=sds, pseudotime=pseudotime, cellWeights=cellWeights,
            nknots=k, verbose=False, seed=seed,
        )
        aics, devs, lls = [], [], []
        for g_models in gams.models:
            for m in g_models:
                if m is None: continue
                try:
                    aics.append(float(m.aic))
                    devs.append(float(m.deviance))
                    lls.append(float(m.llf))
                except Exception:
                    continue
        results.append({
            "k": k,
            "mean_aic": float(np.mean(aics)) if aics else np.nan,
            "mean_deviance": float(np.mean(devs)) if devs else np.nan,
            "mean_loglik": float(np.mean(lls)) if lls else np.nan,
            "n_fits": len(aics),
        })
        if verbose:
            print(f"[evaluateK] k={k}: mean_aic={results[-1]['mean_aic']:.2f}")

    return pd.DataFrame(results).set_index("k")
