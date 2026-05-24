# py-tradeSeq — Math Notes

## 1. Bit-equivalent algorithmic steps

None. Every step depends on the GAM backbone (statsmodels.gam.GLMGam) which is not bit-equivalent to R's mgcv.

## 2. Bounded ε-approximations (B)

**None claimed.** The cross-implementation divergence (statsmodels ≠ mgcv) is documented in §4 as "approximate parity at the inference layer (Spearman ranking)" rather than a derived ε-bound — because the divergence is between two independent GAM solvers, not from a deliberate algebraic simplification in this port.

## 3. Class-containment (C)

None.

## 4. Cross-implementation divergence

### 4.1 `statsmodels.gam.GLMGam` vs R `mgcv::gam`

mgcv uses **Penalized Iteratively Re-weighted Least Squares (P-IRLS)** with **REML** smoothness selection by default. statsmodels.gam.GLMGam uses **Newton-Raphson** + a grid / penalty-search method for the smoothing parameter. For N ≥ ~200 cells per lineage and df ≤ ~6, the two converge to similar fits (Wald-stat Pearson ≈ 0.9). At smaller N or higher df, divergence grows.

**Empirical**: on canonical 140-cell × 100-gene × 2-lineage fixture:
- `associationTest`: Spearman on -log10 p = **0.84**, top-50 Jaccard = **0.67**, Wald-stat Pearson = **0.91**
- `startVsEndTest`: Spearman = **0.45** (lower because endpoint tests are sensitive to the per-lineage independent fit; see 4.2)

### 4.2 Per-lineage independent fits vs mgcv's `by`-smoothers with `id=1`

R tradeSeq fits a single GAM per gene with formula:
```r
y ~ -1 + U + s(t1, by=l1, bs='cr', id=1, k=nknots) + s(t2, by=l2, bs='cr', id=1, k=nknots)
```
The `id=1` ties the smoothness parameter across both lineage smoothers — they share one λ.

Our py port fits **two separate GAMs** (one per lineage's subset of cells) with independent λ. This is a real algorithmic divergence, not a bit-rounding one.

For `associationTest` it matters less (the test sums Wald stats across lineages → similar ranking), but for `startVsEndTest` / `diffEndTest` / `patternTest` the endpoint predictions live in different bases and the contrast variance is mis-estimated. Spearman drops to 0.4-0.5.

**v0.2 fix**: build a single design matrix with both lineage smoothers + a shared L2 penalty term, then fit one statsmodels GLM. ~2 weeks of focused work.

### 4.3 NB theta estimation

R mgcv estimates the NB θ parameter jointly with the smooth coefficients via REML. statsmodels requires `alpha = 1/θ` upfront, so we use method-of-moments from a Poisson preliminary fit. For overdispersed data (typical in scRNA-seq) this introduces a small bias toward higher θ̂.

## 5. Audit class

**B** — minor algorithmic divergence (different GAM solvers + per-lineage vs joint fitting). Inference rankings agree at Spearman > 0.7 typically; documented in `RECONSTRUCTION_REPORT.md §6`.
