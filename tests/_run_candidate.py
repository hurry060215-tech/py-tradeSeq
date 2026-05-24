"""Candidate runner — fits NB GAMs on the same fixture as the R reference.

Reads the sidecar CSVs produced by r_reference_driver.R (pseudotime,
cellWeights, counts) since pyreadr cannot deserialize SlingshotDataSet.
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

_PORT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT_DIR))

from pytradeseq import fitGAM, associationTest, startVsEndTest, diffEndTest


def main():
    fixture_path, output_path = sys.argv[1], sys.argv[2]
    sidecar_dir = Path(output_path).parent

    counts = pd.read_csv(sidecar_dir / "_sidecar_counts.csv", index_col=0).to_numpy(dtype=np.float64)
    pt = pd.read_csv(sidecar_dir / "_sidecar_pseudotime.csv").to_numpy(dtype=np.float64)
    cw = pd.read_csv(sidecar_dir / "_sidecar_cellWeights.csv").to_numpy(dtype=np.float64)
    gene_names = list(pd.read_csv(sidecar_dir / "_sidecar_counts.csv", index_col=0).index)
    print(f"[cand] counts {counts.shape}, pseudotime {pt.shape}, cellWeights {cw.shape}")

    t0 = time.perf_counter()
    gams = fitGAM(counts, pseudotime=pt, cellWeights=cw, nknots=4, verbose=True, seed=42)
    t_fit = time.perf_counter() - t0
    n_conv = int(gams.converged.sum())
    n_total = gams.converged.size
    print(f"[cand] fitGAM: {t_fit:.2f}s; converged {n_conv}/{n_total}")

    t0 = time.perf_counter()
    at = associationTest(gams)
    t_at = time.perf_counter() - t0
    print(f"[cand] associationTest: {t_at:.2f}s; head:")
    print(at.head().to_string())

    t0 = time.perf_counter()
    svet = startVsEndTest(gams)
    t_svet = time.perf_counter() - t0
    print(f"[cand] startVsEndTest: {t_svet:.2f}s")

    t0 = time.perf_counter()
    det = diffEndTest(gams)
    t_det = time.perf_counter() - t0
    print(f"[cand] diffEndTest: {t_det:.2f}s")

    # Top-50 by p-value (lowest first)
    top50_at = at.sort_values("pvalue", na_position="last").index.tolist()[:50]
    top50_svet = svet.sort_values("pvalue", na_position="last").index.tolist()[:50]

    out = {
        "gene_names": gene_names,
        "n_genes": len(gene_names),
        "n_cells": counts.shape[1],
        "n_lineages": pt.shape[1],
        "associationTest": {
            "pvalue": at.loc[gene_names, "pvalue"].tolist(),
            "waldStat": at.loc[gene_names, "waldStat"].tolist(),
            "df": at.loc[gene_names, "df"].astype(int).tolist(),
            "top50": top50_at,
        },
        "startVsEndTest": {
            "pvalue": svet.loc[gene_names, "pvalue"].tolist(),
            "waldStat": svet.loc[gene_names, "waldStat"].tolist(),
            "df": svet.loc[gene_names, "df"].astype(int).tolist(),
            "top50": top50_svet,
        },
        "timings": {
            "fitGAM": t_fit,
            "associationTest": t_at,
            "startVsEndTest": t_svet,
            "diffEndTest": t_det,
        },
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote {output_path}")


if __name__ == "__main__":
    main()
