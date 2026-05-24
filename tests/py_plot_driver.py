"""Render py-tradeSeq plotSmoothers on the same fixture."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))

import pytradeseq
from pytradeseq.plotting import plotSmoothers
from ggplot2_py import ggsave


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(_PORT / "data/_sidecar_counts.csv", index_col=0)
    pseudotime = pd.read_csv(_PORT / "data/_sidecar_pseudotime.csv").to_numpy(dtype=np.float64)
    weights = pd.read_csv(_PORT / "data/_sidecar_cellWeights.csv").to_numpy(dtype=np.float64)
    print(f"counts: {counts.shape}, pseudotime: {pseudotime.shape}, weights: {weights.shape}")

    # counts is genes × cells (R orientation). pytradeseq.fitGAM expects same
    gam = pytradeseq.fitGAM(
        counts=counts,
        pseudotime=pseudotime,
        cellWeights=weights,
        nknots=5,
        verbose=False,
    )

    # Read R's top gene
    top_gene = pd.read_csv(out_dir / "top_gene.csv")["gene"].iloc[0]
    print(f"using gene: {top_gene}")

    counts_T = counts.T  # cells × genes for plotting helper
    p = plotSmoothers(gam, counts=counts_T, gene=top_gene)
    ggsave(str(out_dir / "Py_smoother.png"), plot=p, width=6, height=4, dpi=100)
    print("done")


if __name__ == "__main__":
    main()
