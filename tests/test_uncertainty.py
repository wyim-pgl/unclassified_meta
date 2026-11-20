"""em-uncertainty: read bootstrap, posterior realizations, alternative mapping criteria and their seal."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import manifest, panel, uncertainty

FEATURES = ["ns|A|same", "ns|B|like", "ns|C|same", "ns|D|like"]
DECOY = "decoy_x"
# D is only ever seen in shared classes; A and B share a parent; C and D have none
CLASSES = {
    "s1": {("ns|A|same",): 200, ("ns|B|like",): 5, ("ns|A|same", "ns|B|like"): 100, ("ns|C|same",): 50,
           ("ns|B|like", "ns|D|like"): 40, ("ns|D|like", DECOY): 20, (DECOY,): 10, (): 25},
    "s2": {("ns|A|same",): 150, ("ns|B|like",): 12, ("ns|A|same", "ns|B|like"): 60, ("ns|C|same",): 80,
           ("ns|B|like", "ns|D|like"): 30, ("ns|D|like", DECOY): 15, (DECOY,): 8, (): 40},
}


def quant(tmp_path, name="quant"):
    """A sealed panel-quant folder with equivalence classes for samples s1 and s2."""
    q = tmp_path / name
    q.mkdir()
    samples = [SimpleNamespace(name=s) for s in CLASSES]
    pd.DataFrame({"feature_id": FEATURES, "label": ["A", "B", "C", "D-like"],
                  "identity_class": ["same", "like", "same", "like"],
                  "mixed_identity": [False, False, False, True],
                  "parent_feature": ["P", "P", "", ""]}).to_csv(q / "features.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": FEATURES, "s1": [200, 5, 50, 0], "s2": [150, 12, 80, 0]}).to_csv(
        q / "counts.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": list(CLASSES), "input_fragments": [sum(c.values()) for c in CLASSES.values()]}).to_csv(
        q / "ledger.tsv", sep="\t", index=False)
    panel.write_class_tables(q, samples, CLASSES, FEATURES)
    panel.write_quant_provenance(q, {"panel": {"panel_id": "p1"}, "ambiguous": "ledger"}, list(panel.CLASS_ARTIFACTS))
    return q


def run(tmp_path, out="unc", **kw):
    q = tmp_path / "quant" if (tmp_path / "quant").exists() else quant(tmp_path)
    kw = {"bootstraps": 4, "realizations": 3, "seed": 1, **kw}
    return uncertainty.run(q, tmp_path / out, **kw)


def table(path):
    return pd.read_csv(path, sep="\t").set_index("feature_id")


def test_bootstrap_replicate_keeps_the_input_fragment_total():
    rng = np.random.default_rng([1, 3, 0])
    for _ in range(5):
        rep = uncertainty.resample_classes(CLASSES["s1"], rng)
        assert sum(rep.values()) == sum(CLASSES["s1"].values())
        assert set(rep) <= set(CLASSES["s1"])


def test_bootstrap_is_reproducible_with_the_same_seed_and_differs_with_another(tmp_path):
    a, b, c = run(tmp_path, "a"), run(tmp_path, "b"), run(tmp_path, "c", seed=2)
    for n in range(1, 5):
        name = uncertainty.bootstrap_name(n)
        assert (a / name).read_bytes() == (b / name).read_bytes()
    assert any((a / uncertainty.bootstrap_name(n)).read_bytes() != (c / uncertainty.bootstrap_name(n)).read_bytes()
               for n in range(1, 5))


def test_results_do_not_depend_on_threads(tmp_path):
    one, two = run(tmp_path, "one", threads=1), run(tmp_path, "two", threads=2)
    for name in (uncertainty.bootstrap_name(2), uncertainty.realization_name(3), "em_bootstrap_abundance.tsv"):
        assert (one / name).read_bytes() == (two / name).read_bytes()


def test_unique_counts_are_the_quant_table_unaltered(tmp_path):
    out = run(tmp_path)
    assert (out / "unique_counts.tsv").read_bytes() == (tmp_path / "quant/unique_counts.tsv").read_bytes()
    em = table(out / "em_expected_counts.tsv")
    uniq = table(out / "unique_counts.tsv")
    assert uniq.loc["ns|A|same", "s1"] == 200 and em.loc["ns|A|same", "s1"] > 200   # EM adds, never replaces


def test_feature_supported_only_by_shared_fragments_is_em_driven_and_low(tmp_path):
    sup = uncertainty.read_support(run(tmp_path) / "viral_unique_read_support.tsv").set_index("feature_id")
    d = sup.loc["ns|D|like"]
    assert d.unique_total == 0 and d.em_total > 0 and bool(d.em_driven)
    assert d.assignment_confidence == "low" and d.samples_with_unique_support == 0
    assert sup.loc["ns|A|same", "assignment_confidence"] in ("high", "moderate")


def test_alternative_criteria_bracket_the_em_counts(tmp_path):
    a = pd.read_csv(run(tmp_path) / "em_assignment_uncertainty.tsv", sep="\t")
    tol = 1e-6
    assert (a["unique"] <= a["em"] + tol).all() and (a["em"] <= a["any_candidate"] + tol).all()
    assert (a["unique"] <= a["even_split"] + tol).all() and (a["even_split"] <= a["any_candidate"] + tol).all()
    row = a[(a.feature_id == "ns|B|like") & (a["sample"] == "s1")].iloc[0]
    assert row.shared == 140 and row.even_split == pytest.approx(5 + 50 + 20)
    assert row.em_from_shared == pytest.approx(row.em - row["unique"])


def test_realizations_are_integers_that_keep_every_mapped_fragment(tmp_path):
    mapped = {k: v for k, v in CLASSES["s1"].items() if k}
    res = uncertainty._sample_task(("s1", 0, CLASSES["s1"], 1, 2, 5))
    for r in res["real"]:
        assert all(isinstance(v, int) for v in r.values()) and sum(r.values()) == sum(mapped.values())
    out = run(tmp_path)
    for n in range(1, 4):
        t = table(out / uncertainty.realization_name(n))
        assert all(t[c].dtype.kind == "i" for c in t.columns)
        assert t["s1"].sum() <= sum(mapped.values())          # decoy draws are not feature rows


def test_normalized_abundance_is_per_million_input_fragments(tmp_path):
    out = run(tmp_path)
    em, norm = table(out / "em_expected_counts.tsv"), table(out / "normalized_abundance.tsv")
    total = sum(CLASSES["s2"].values())
    assert norm.loc["ns|C|same", "s2"] == pytest.approx(em.loc["ns|C|same", "s2"] / total * 1e6, abs=1e-5)


def test_manifest_seals_every_output_and_a_tampered_file_fails_verify(tmp_path):
    out = run(tmp_path)
    rec = uncertainty.verify_uncertainty(out)
    on_disk = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()} - {manifest.NAME, manifest.MARKER}
    assert set(rec["artifacts"]) == on_disk
    assert rec["step"] == "em-uncertainty" and rec["bootstraps"] == 4 and rec["realizations"] == 3
    assert rec["quant_provenance_sha256"] == manifest.sha256(tmp_path / "quant/provenance.json")
    with open(out / uncertainty.bootstrap_name(2), "a") as fh:
        fh.write("x\n")
    with pytest.raises(manifest.ManifestError, match="changed after"):
        uncertainty.verify_uncertainty(out)


def test_rerun_with_fewer_bootstraps_removes_the_stale_replicates(tmp_path):
    run(tmp_path)
    out = run(tmp_path, bootstraps=2)
    assert not (out / uncertainty.bootstrap_name(3)).exists()
    assert uncertainty.verify_uncertainty(out)["bootstraps"] == 2


def test_aggregation_sums_by_parent_feature(tmp_path):
    out = run(tmp_path)
    em, agg = table(out / "em_expected_counts.tsv"), table(out / "aggregated_counts.tsv")
    assert list(agg.index) == ["P", "ns|C|same", "ns|D|like"]
    assert agg.loc["P", "s1"] == pytest.approx(em.loc["ns|A|same", "s1"] + em.loc["ns|B|like", "s1"], abs=1e-5)
    au = table(out / "aggregated_unique_counts.tsv")
    assert au.loc["P", "s2"] == 162
    assert manifest.verify(out)["aggregate_by"] == "parent_feature"


def test_taxonomy_without_a_manifest_is_refused(tmp_path):
    tx = tmp_path / "tx" / "panel_taxonomy.tsv"
    tx.parent.mkdir()
    pd.DataFrame({"feature_id": FEATURES, "taxid": ["1", "1", "", ""]}).to_csv(tx, sep="\t", index=False)
    with pytest.raises(manifest.ManifestError, match="panel-taxonomy"):
        run(tmp_path, taxonomy=tx)


def test_taxonomy_aggregates_features_by_taxid(tmp_path):
    tx = tmp_path / "tx" / "panel_taxonomy.tsv"
    tx.parent.mkdir()
    pd.DataFrame({"feature_id": FEATURES, "taxid": ["", "9", "9", ""]}).to_csv(tx, sep="\t", index=False)
    manifest.seal(tx.parent, ["panel_taxonomy.tsv"], {"step": "panel-taxonomy", "panel_id": "p1"})
    out = run(tmp_path, taxonomy=tx)
    agg = table(out / "aggregated_unique_counts.tsv")
    assert list(agg.index) == ["ns|A|same", "taxid:9", "ns|D|like"] and agg.loc["taxid:9", "s1"] == 55


def test_a_tampered_quant_folder_is_refused(tmp_path):
    q = quant(tmp_path)
    with open(q / "equivalence_classes.tsv", "a") as fh:
        fh.write("s1\tns|A|same\t1\n")
    with pytest.raises(manifest.ManifestError):
        uncertainty.run(q, tmp_path / "unc", bootstraps=2, realizations=1)
