"""Parity gate against R tradeSeq on canonical fixture."""
import json, sys, subprocess
from pathlib import Path
import numpy as np
import pytest
import yaml

PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT))
sys.path.insert(0, str(PORT.parent / "omicverse-rebuildr" / "engine"))
from parity_metrics import compute_parity


@pytest.fixture(scope="session")
def manifest():
    return yaml.safe_load((PORT / "data" / "manifest.yaml").read_text())


@pytest.fixture(scope="session")
def outputs():
    ref = PORT / "data" / "reference_output.json"
    cand = PORT / "data" / "candidate_output.json"
    if not (ref.exists() and cand.exists()):
        pytest.skip("Run r_reference_driver.R + _run_candidate.py first")
    return json.loads(ref.read_text()), json.loads(cand.read_text())


def test_associationTest_parity(manifest, outputs):
    r, p = outputs
    r_pval = np.array(r["associationTest"]["pvalue"])
    p_pval = np.array(p["associationTest"]["pvalue"])
    m = compute_parity(r_pval, p_pval, "inference", top_k=50)
    spec = next(o for o in manifest["outputs"] if o["name"] == "associationTest_pval")
    assert m["spearman_neglog10p"] >= spec["threshold"], (
        f"associationTest Spearman {m['spearman_neglog10p']:.3f} < {spec['threshold']}"
    )


def test_associationTest_top50(manifest, outputs):
    r, p = outputs
    r_top = list(r["associationTest"]["top50"])
    p_top = list(p["associationTest"]["top50"])
    jacc = len(set(r_top) & set(p_top)) / max(1, len(set(r_top) | set(p_top)))
    spec = next(o for o in manifest["outputs"] if o["name"] == "associationTest_top50_genes")
    assert jacc >= spec["threshold"], f"top50 Jaccard {jacc:.3f} < {spec['threshold']}"
