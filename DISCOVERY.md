# Discovery — py-tradeSeq

## 1. Is this package already ported?

```
$ python -m engine.discover_omicverse_deps --check tradeSeq
**No existing omicverse port found.** Safe to start a new port.
```

**Decision**: START_PORT.

## 2. Dependency audit

```
$ python -m engine.discover_omicverse_deps --description tradeSeq-ref/DESCRIPTION

| R dep                  | omicverse match | Decision |
|------------------------|-----------------|----------|
| mgcv                   | —               | The hard part. statsmodels.gam OR pygam — both have documented divergence from mgcv NB-GAM. |
| edgeR                  | py-edgeR ✅     | hard dep — for offset / library-size normalization |
| SingleCellExperiment   | —               | native-python: AnnData |
| SummarizedExperiment   | —               | native-python: AnnData |
| slingshot              | —               | omicverse.single._pyslingshot (already in omicverse main package) |
| princurve              | —               | scipy.interpolate.splprep + scipy.optimize |
| BiocParallel           | —               | native-python: joblib / concurrent.futures |
| Biobase                | —               | native-python: AnnData |
| pbapply                | —               | tqdm |
| igraph                 | —               | networkx |
| ggplot2 / RColorBrewer / viridis | —    | out-of-scope (plotting deferred to v0.2) |
| Matrix                 | —               | scipy.sparse |
| matrixStats            | —               | scipy.stats / numpy |
| MASS                   | —               | scipy.stats |
| TrajectoryUtils        | —               | porting alongside is overkill; reimplement minimal helpers |
| tibble                 | —               | pandas |
| magrittr               | —               | native python (no need) |
| S4Vectors              | —               | not needed (we don't ship S4 objects) |
| methods                | —               | not needed (we don't ship S4 objects) |
```

## 3. Decisions per R dep

| R dep | omicverse match | Decision | Python replacement |
|---|---|---|---|
| `mgcv` | — | **GAM-backbone** | **`statsmodels.gam` with NB family + spline smoothers**; documented divergence; inference-class threshold 0.70 (relaxed from 0.90 default) |
| `edgeR` | `py-edgeR` v0.x | hard dep (optional) | `pyedgeR>=0.1` — for library-size normalization helpers (otherwise compute inline) |
| `slingshot` | — | reuse via omicverse | accept upstream pseudotime + cell weights as inputs, don't reimplement slingshot |
| `BiocParallel` | — | native-python | `joblib.Parallel` |
| `princurve` | — | native-python | `scipy.interpolate.splprep` |
| Plotting (ggplot2 / RColorBrewer / viridis) | — | out-of-scope | matplotlib (deferred to v0.2) |
| `pbapply` | — | native-python | `tqdm` |
| Trajectory containers (SCE / SummarizedExperiment / Biobase) | — | native-python | AnnData |

## 4. Reusable work saved

| Reused omicverse port | LOC saved | Notes |
|---|---|---|
| `py-edgeR` | ~200 (library-size norm, dispersion) | optional dep |
| `omicverse.single._pyslingshot` | n/a | upstream consumer — pseudotime input comes from user |
| **Total** | **~200** | small reuse; tradeSeq's algorithmic core is mgcv-bound |

## 5. New ports surfaced

A future **`py-mgcv`** port would let py-tradeSeq tighten its inference threshold from 0.70 to 0.90. mgcv is a large package (~30k LOC R) — major engineering effort, not in scope here. Added to `examples/ROADMAP_TRAJ.md` as a future tier.

## 6. Key technical risk

**mgcv vs statsmodels.gam NB-GAM divergence** is real and documented:
- mgcv uses penalized iteratively-reweighted least-squares (P-IRLS) with REML
- statsmodels uses Newton-Raphson + scoring with various penalty/smoothing parameter selection
- For small N (cells per lineage < 100) or rare features (low-count genes), the fits can differ by 10-30% on Wald statistics
- p-value ordering is generally preserved (Spearman > 0.7 typically), but individual p-values are NOT bit-equivalent

This is why we pre-register threshold 0.70 (Spearman on -log10 p) — realistic for cross-implementation GAM agreement, not aspirational.

## 7. Scope ranking

- **v0.1 (mandatory)**: `fitGAM`, `associationTest`, `startVsEndTest`, `diffEndTest`, `nknots`
- **v0.2**: `patternTest`, `earlyDETest`, `evaluateK`
- **v0.3**: `conditionTest`, `clusterExpressionPatterns`, `cascade`
- **v0.4**: plotting (`plotGeneCount`, `plotSmoothers`, `plot_evalutateK_results`)

This is the largest port attempted in the trajectory roadmap (~5300 R LOC vs TSCAN's 640). Expected effort: 5-10 active days for v0.1.
