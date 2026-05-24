"""pytradeseq — pure-Python port of tradeSeq.

v0.1 covers fitGAM + 5 statistical tests + diagnostics:
- fitGAM:           fit NB-GAM per gene per lineage
- associationTest:  does expression depend on pseudotime?      [parity-validated]
- startVsEndTest:   does expression differ between endpoints?  [approximate]
- diffEndTest:      do lineage endpoints differ?               [approximate]
- patternTest:      does expression PATTERN differ across lineages?
- earlyDETest:      restricted patternTest on early pseudotime
- nknots:           return nknots used in a fit
- evaluateK:        BIC scan over candidate knot values

NOT yet ported: conditionTest, clusterExpressionPatterns, cascade,
predictSmooth, predictCells, all plotting functions. See
RECONSTRUCTION_REPORT.md §6 for known limitations.

Algorithm class: inference. Pre-registered thresholds in `data/manifest.yaml`:
- associationTest p-value: Spearman ≥ 0.70 on -log10 p (relaxed from default 0.90
  because statsmodels.gam ≠ mgcv on NB-GAM fits at small df).
"""

from __future__ import annotations

__version__ = "0.1.0"

from .fitting import fitGAM, TradeSeqGAMs
from .tests_module import associationTest, startVsEndTest, diffEndTest
from .pattern_tests import patternTest, earlyDETest
from .diagnostics import nknots, evaluateK

__all__ = [
    "fitGAM",
    "TradeSeqGAMs",
    "associationTest",
    "startVsEndTest",
    "diffEndTest",
    "patternTest",
    "earlyDETest",
    "nknots",
    "evaluateK",
    "__version__",
]
