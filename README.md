# py-tradeSeq

A **Python port of [tradeSeq](https://github.com/statOmics/tradeSeq)** (Van den Berge et al., *Nature Communications* 2020) — trajectory-based differential expression analysis for single-cell RNA-seq.

- AnnData-compatible
- Fits NB GAM per gene per lineage via `statsmodels.gam.GLMGam`
- Inference-class parity: Spearman = **0.84** on -log10 p, top-50 Jaccard = **0.67** vs R `associationTest`
- 5× faster than R `fitGAM` on canonical fixture

## Install

```bash
pip install pytradeseq
```

## Quick-start

```python
import numpy as np
from pytradeseq import fitGAM, associationTest, startVsEndTest

# pseudotime: cells × lineages; cellWeights: cells × lineages
gams = fitGAM(counts, pseudotime=pt, cellWeights=cw, nknots=6)

# Which genes vary along pseudotime?
at = associationTest(gams)
print(at.sort_values('pvalue').head(20))

# Which genes differ between trajectory endpoints?
svet = startVsEndTest(gams)
```

## Function map

| Python | R counterpart | Purpose |
|---|---|---|
| `fitGAM` | `fitGAM` | fit NB-GAM per gene per lineage |
| `associationTest` | `associationTest` | Wald test against null (pseudotime-independent) — **parity-validated** |
| `startVsEndTest` | `startVsEndTest` | endpoint comparison — approximate |
| `diffEndTest` | `diffEndTest` | between-lineage endpoint comparison |
| `patternTest` | `patternTest` | between-lineage curve-shape comparison |
| `earlyDETest` | `earlyDETest` | restricted patternTest on early pt |
| `nknots` | `nknots` | return knot count of a fit |
| `evaluateK` | `evaluateK` | AIC scan over candidate knots |

## Parity status

| Test | Metric | Threshold | Measured | Pass |
|---|---|---|---|---|
| `associationTest` p-value | Spearman on -log10 p | ≥ 0.70 | **0.84** | ✅ |
| `associationTest` top-50 | Jaccard | ≥ 0.60 | **0.67** | ✅ |
| `startVsEndTest` | Spearman | — | 0.45 (approximate) | 🟡 |

The `startVsEndTest` lower number reflects a known divergence from R's joint `id=1` tied-smoothness fitting. v0.2 will switch to a joint-GLM formulation.

## Known limitations (v0.1)

1. **`statsmodels.gam` ≠ `mgcv`**: different penalty + smoothness-selection algorithms. Inference rankings agree (Spearman > 0.7 on associationTest) but individual p-values are NOT bit-equivalent.
2. **`startVsEndTest` / `diffEndTest`**: per-lineage independent smoother fitting differs from R's joint-fit-with-tied-smoothness. v0.2.
3. **9/18 functions deferred** to v0.2+: `conditionTest`, `clusterExpressionPatterns`, `cascade`, `predictSmooth`, `predictCells`, all plotting.
4. **No AnnData class API yet** — functional API only. v0.2.

## Citation

> Van den Berge, K. et al. **Trajectory-based differential expression analysis for single-cell sequencing data.** *Nature Communications* 11, 1201 (2020).

## License

MIT.
