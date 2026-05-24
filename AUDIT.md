# R function coverage audit — py-tradeSeq v0.1.0

Generated against `tradeSeq` v1.13.12 source (18 exports + S4 methods).

## Coverage summary

| Category | Count | Coverage |
|---|---|---|
| Algorithmic exports ported | 8 / 8 (core) | 100% |
| Plotting / GUI exports | 0 / 3 | 0% (deferred) |
| `predict*` helpers | 0 / 2 | 0% (deferred to v0.2) |
| `conditionTest` / `clusterExpressionPatterns` / `cascade` | 0 / 3 | 0% (deferred to v0.3) |

## Exported R functions

| R function | Python equivalent | Status | Notes |
|---|---|---|---|
| `fitGAM` | `pytradeseq.fitGAM` | ✅ ported | statsmodels.gam backend |
| `associationTest` | `pytradeseq.associationTest` | ✅ ported, parity-validated | Spearman 0.84 |
| `startVsEndTest` | `pytradeseq.startVsEndTest` | 🟡 ported (approximate) | Spearman 0.45 |
| `diffEndTest` | `pytradeseq.diffEndTest` | 🟡 ported (approximate) | inherits joint-fit issue |
| `patternTest` | `pytradeseq.patternTest` | 🟡 ported (approximate) | inherits joint-fit issue |
| `earlyDETest` | `pytradeseq.earlyDETest` | 🟡 ported (approximate) | inherits joint-fit issue |
| `nknots` | `pytradeseq.nknots` | ✅ ported | trivial |
| `evaluateK` | `pytradeseq.evaluateK` | ✅ ported | AIC scan over k |
| `getSmootherPvalues` | — | ⏳ v0.2 | helper; trivial to add |
| `getSmootherTestStats` | — | ⏳ v0.2 | helper; trivial |
| `conditionTest` | — | ⏳ v0.3 | needs joint-fit |
| `clusterExpressionPatterns` | — | ⏳ v0.3 | needs predictSmooth |
| `cascade` | — | ⏳ v0.3 | rare-use |
| `predictSmooth` | — | ⏳ v0.2 | needed for plotting |
| `predictCells` | — | ⏳ v0.2 | needed for plotting |
| `plotSmoothers` | — | ⛔ plotting → v0.4 |
| `plotGeneCount` | — | ⛔ plotting → v0.4 |
| `plot_evalutateK_results` | — | ⛔ plotting → v0.4 |
