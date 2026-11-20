"""Compare published MI tests with a fresh pooling by the verified cabe340 export."""
from pathlib import Path
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "pinned/src"))
from unclassified_meta.da import combine_imputations, read_results

TERM = "grouptreated"


def verify_tree():
    for line in (ROOT / "pinned/source_blobs.tsv").read_text().splitlines()[1:]:
        blob, name = line.split("\t")
        b = (ROOT / "pinned" / name).read_bytes()
        assert hashlib.sha1(f"blob {len(b)}\0".encode() + b).hexdigest() == blob, name


def compare(folder):
    paths = sorted((folder / "da_em/imputations").glob("ancombc2_results_*.tsv"))
    assert len(paths) == 5
    inputs = [read_results(p) for p in paths]
    expected, terms = combine_imputations(inputs, .05, "BH")
    assert terms == [TERM]
    stats = [c for c in expected if c.endswith("_" + TERM)]
    actual = read_results(folder / "da_em/mi_results.tsv").set_index("feature_id").sort_index()
    expected = expected.set_index("feature_id").sort_index()
    assert actual.index.equals(expected.index)
    diffs = {}
    for c in stats:
        x, y = expected[c].astype(float), actual[c].astype(float)
        if c.startswith(("diff_", "n_fits_")):
            np.testing.assert_array_equal(x, y)
        else:
            np.testing.assert_allclose(x, y, rtol=2e-8, atol=1e-10)
        diffs[c] = float((x - y).abs().replace(np.inf,np.nan).max())
    rare = sum(int(((r["p_"+TERM] >= 1) & r["W_"+TERM].ne(0)).sum()) for r in inputs)
    return {"path": str(folder.relative_to(ROOT)), "rows": len(expected), "p1_nonzero_w": rare, **diffs}


def verify_duplicate(folder):
    alt = folder / "pinned_fit"
    for rel in ["lca/counts.tsv","em/counts.tsv", *(str(p.relative_to(folder)) for p in sorted((folder/"em/imputations").glob("counts_*.tsv")))]:
        assert (folder / rel).read_bytes() == (alt / rel).read_bytes(), (folder,rel)
    for rel in ["da_lca/ancombc2_results.tsv","da_em/ancombc2_results.tsv","da_em/mi_results.tsv",
                *(str(p.relative_to(folder)) for p in sorted((folder/"da_em/imputations").glob("ancombc2_results_*.tsv")))]:
        x = read_results(folder/rel).sort_values("feature_id").reset_index(drop=True)
        y = read_results(alt/rel).sort_values("feature_id").reset_index(drop=True)
        pd.testing.assert_frame_equal(x, y, rtol=2e-8, atol=1e-10)


if __name__ == "__main__":
    verify_tree()
    rows = []
    for dataset in sorted((ROOT/"runs").iterdir()):
        for base in (dataset, dataset/"pinned_fit"):
            p = base / "run.json"
            if not p.exists() or json.loads(p.read_text())["status"] != "ok":
                continue
            check = base/"repool_check.json"
            if not check.exists():
                info = compare(base)
                check.write_text(json.dumps(info,indent=2)+"\n")
            rows.append(json.loads(check.read_text()))
        if ((dataset/"pinned_fit/run.json").exists() and json.loads((dataset/"pinned_fit/run.json").read_text())["status"] == "ok"
                and json.loads((dataset/"run.json").read_text())["status"] == "ok"):
            verify_duplicate(dataset)
    pd.DataFrame(rows).to_csv(ROOT/"repool_audit.tsv",sep="\t",index=False)
    print("Verified",len(rows),"pooled results from their five component fits.",flush=True)
