"""`da --rand-formula`: repeated measures of the same unit (needs the R environment; skipped without it)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import da

HERE = Path(__file__).parent
TRUE = {"taxon_003", "taxon_017"}


def rscript():
    try:
        r = da.find_rscript(None)
    except Exception:
        return None
    return r if da.has_ancombc(r) else None


def cli(*args):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("mixed")
    subprocess.run([sys.executable, str(HERE / "make_synthetic_counts.py"), str(d)], check=True)
    return d


def base(data, out, r):
    return ("da", "--unsealed-input", "--counts", data / "counts.tsv", "--features", data / "features.tsv",
            "--metadata", data / "metadata.tsv", "--group", "group", "--reference-level", "H", "-o", out,
            "--rscript", r, "-t", 4)


def test_random_unit_effect_finds_the_true_features_without_false_positives(data, tmp_path):
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package (create the environment-r.yml environment)")
    rc, log = cli(*base(data, tmp_path / "da", r), "--rand-formula", "(1|unit)")
    assert rc == 0, log
    text = (tmp_path / "da" / "ancombc2_results.tsv").read_text()
    assert "# rand_formula: (1|unit) [unit: 12 units, 2-2 samples per unit, units per group H=6, M=6]" in text
    res = pd.read_csv(tmp_path / "da" / "ancombc2_results.tsv", sep="\t", comment="#").set_index("feature_id")
    hits = set(res.index[res["diff_robust_groupM"].astype(bool)])
    assert TRUE <= hits and len(hits - TRUE) == 0, hits
    assert (res.loc[sorted(TRUE), "lfc_groupM"] > 1.3).all()            # true value ln(8) = 2.08


def test_random_effect_column_is_checked_before_fitting(data, tmp_path):
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package")
    rc, log = cli(*base(data, tmp_path / "a", r), "--rand-formula", "(1|plant)")
    assert rc == 1 and "not in the metadata: plant" in log
    rc, log = cli(*base(data, tmp_path / "b", r), "--rand-formula", "(1|sample)")
    assert rc == 1 and "one sample per unit" in log
    rc, log = cli(*base(data, tmp_path / "c", r), "--rand-formula", "unit")
    assert rc == 1 and "must look like" in log
    rc, log = cli(*base(data, tmp_path / "d", r), "--lme-relax-checks")
    assert rc == 1 and "needs --rand-formula" in log
