"""Smoke tests for pytradeseq."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pytradeseq


def test_import():
    assert pytradeseq.__version__.startswith("0.1")
    for fn in ("fitGAM", "associationTest", "startVsEndTest", "diffEndTest",
               "patternTest", "earlyDETest", "nknots", "evaluateK"):
        assert hasattr(pytradeseq, fn)


def test_fitGAM_minimal():
    rng = np.random.default_rng(42)
    n_cells, n_genes, n_lineages = 50, 10, 2
    pt = rng.uniform(0, 1, (n_cells, n_lineages))
    cw = rng.dirichlet([1, 1], n_cells)
    counts = rng.negative_binomial(5, 0.5, (n_genes, n_cells)).astype(float)
    gams = pytradeseq.fitGAM(counts, pseudotime=pt, cellWeights=cw,
                              nknots=4, verbose=False, seed=42)
    assert gams.n_genes == n_genes
    assert gams.n_lineages == n_lineages
    # At least some genes converged
    assert gams.converged.sum() > 0


def test_associationTest():
    rng = np.random.default_rng(42)
    n_cells = 60
    pt = rng.uniform(0, 1, (n_cells, 1))
    cw = np.ones((n_cells, 1))
    counts = rng.negative_binomial(5, 0.5, (5, n_cells)).astype(float)
    # Inject pseudotime signal into gene 0
    counts[0] = rng.negative_binomial(5 + 8 * pt[:, 0], 0.5)
    gams = pytradeseq.fitGAM(counts, pseudotime=pt, cellWeights=cw,
                              nknots=4, verbose=False, seed=42)
    at = pytradeseq.associationTest(gams)
    assert at.shape[0] == 5
    assert "pvalue" in at.columns
    assert (at["pvalue"] >= 0).all() and (at["pvalue"] <= 1).all()
