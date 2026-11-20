"""da-robustness / da-resolutions: shared feature set, stability arithmetic, concordance, leave-one-out, the
conclusion rule, the wording rule and the seal (R replaced by a fake), plus one run of the real R scripts."""
import os
import shutil
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import da, manifest, robustness as rb

SAMPLES = [f"A{i}" for i in range(1, 7)] + [f"B{i}" for i in range(1, 7)]
UP = "f05"                       # six-fold more abundant in group B
RARE = "f30"                     # in one sample only: below the shared prevalence filter


def synthetic(rng, n_feat=30, fold=8.0) -> pd.DataFrame:
    base = rng.uniform(20, 200, n_feat)
    rows = {}
    for i in range(n_feat):
        fid = f"f{i + 1:02d}"
        mu = np.array([base[i] * (fold if (fid == UP and s.startswith("B")) else 1.0) for s in SAMPLES])
        rows[fid] = rng.poisson(mu).astype(float)
    df = pd.DataFrame(rows, index=SAMPLES).T
    df.loc[RARE] = 0.0
    df.loc[RARE, "A1"] = 5.0
    df.index.name = "feature_id"
    return df


def write_counts(df: pd.DataFrame, path: Path, decimals=6) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.round(decimals).reset_index().to_csv(path, sep="\t", index=False)


def em_folder(tmp_path, n_boot=4, n_real=3, low_conf=(), seed=0) -> tuple[Path, Path, Path]:
    """A sealed em-uncertainty folder, a sealed features table and a metadata table."""
    rng = np.random.default_rng(seed)
    em = tmp_path / "emu"
    em.mkdir()
    primary = synthetic(rng) * rng.uniform(0.97, 1.03, (30, 12))     # fractional, like EM expected counts
    write_counts(primary, em / "em_expected_counts.tsv")
    write_counts(primary.round(), em / "unique_counts.tsv")
    write_counts(primary / primary.sum() * 1e6, em / "normalized_abundance.tsv")
    files = ["em_expected_counts.tsv", "unique_counts.tsv", "normalized_abundance.tsv"]
    for i in range(1, n_boot + 1):
        write_counts(primary * rng.uniform(0.9, 1.1, primary.shape), em / f"bootstrap/counts_{i:03d}.tsv")
        files.append(f"bootstrap/counts_{i:03d}.tsv")
    for i in range(1, n_real + 1):
        write_counts((primary * rng.uniform(0.95, 1.05, primary.shape)).round(), em / f"realizations/counts_{i:03d}.tsv",
                     0)
        files.append(f"realizations/counts_{i:03d}.tsv")
    sup = pd.DataFrame({"feature_id": primary.index, "unique_fraction": 0.9, "boot_cv_median": 0.05,
                        "em_driven": False,
                        "assignment_confidence": ["low" if f in low_conf else "high" for f in primary.index]})
    with open(em / "viral_unique_read_support.tsv", "w") as fo:
        fo.write("# assignment_confidence rule: ...\n")
        sup.to_csv(fo, sep="\t", index=False)
    ba = [{"feature_id": f, "label": f, "sample": s, "em": primary.loc[f, s], "boot_mean": primary.loc[f, s],
           "boot_sd": 1, "boot_p2_5": 0, "boot_p97_5": 1, "boot_cv": 0.05} for f in primary.index for s in SAMPLES]
    pd.DataFrame(ba).to_csv(em / "em_bootstrap_abundance.tsv", sep="\t", index=False)
    files += ["viral_unique_read_support.tsv", "em_bootstrap_abundance.tsv"]
    for name in ("aggregated_counts.tsv", "aggregated_unique_counts.tsv", "aggregated_groups.tsv",
                 "em_assignment_uncertainty.tsv"):                    # sealed by em-uncertainty, not read here
        (em / name).write_text("feature_id\n")
        files.append(name)
    manifest.seal(em, files, {"step": "em-uncertainty", "bootstraps": n_boot, "realizations": n_real, "seed": 7})
    quant = tmp_path / "quant"
    quant.mkdir()
    feats = pd.DataFrame({"feature_id": primary.index, "label": [f"virus {f}" for f in primary.index],
                          "parent_feature": [f"bin{int(f[1:]) % 10}" for f in primary.index]})
    feats.to_csv(quant / "features.tsv", sep="\t", index=False)
    manifest.seal(quant, ["features.tsv"], {"step": "panel-quant"})
    meta = tmp_path / "meta.tsv"
    pd.DataFrame({"sample": SAMPLES, "condition": ["ctrl"] * 6 + ["macho"] * 6,
                  "plant": [f"p{i // 2}" for i in range(12)],
                  "batch": (["b1"] * 3 + ["b2"] * 3) * 2}).to_csv(meta, sep="\t", index=False)
    return em, quant / "features.tsv", meta


class FakeR:
    """Stands in for the three R scripts: effect = log2 ratio of group means (centred over features); a feature
    with |effect| > 1.5 gets p = 1e-4. Records the features every call received; can fail chosen fits."""

    def __init__(self, fail=(), alter=None):
        self.calls, self.fail, self.alter = [], set(fail), alter

    def __call__(self, cmd, capture_output=True, text=True):
        if cmd[1] == "-e":                                         # package check
            return subprocess.CompletedProcess(cmd, 0, "", "")
        arg = lambda k: cmd[cmd.index(k) + 1] if k in cmd else None
        method = Path(cmd[1]).stem
        out = Path(arg("--outdir"))
        counts = pd.read_csv(arg("--counts"), sep="\t", dtype={"feature_id": str}).set_index("feature_id")
        meta = pd.read_csv(arg("--metadata"), sep="\t", dtype=str)
        self.calls.append({"method": method, "analysis": out.name.split("__", 1)[1], "features": list(counts.index),
                           "samples": list(meta["sample"]), "cmd": cmd})
        if (method, out.name.split("__", 1)[1]) in self.fail:
            return subprocess.CompletedProcess(cmd, 1, "", "Error: fake failure\n")
        g = meta.set_index("sample")["condition"]
        a = counts[[s for s in counts.columns if g.get(s) == "ctrl"]].mean(axis=1)
        b = counts[[s for s in counts.columns if g.get(s) == "macho"]].mean(axis=1)
        eff = np.log2((b + 0.5) / (a + 0.5))
        eff = eff - eff.median()
        p = np.where(eff.abs() > 1.5, 1e-4, 0.5)
        res = pd.DataFrame({"feature_id": counts.index, "term": "conditionmacho", "effect": eff.to_numpy(),
                            "se": 0.3, "ci_low": eff - 0.6, "ci_high": eff + 0.6, "stat": eff / 0.3, "p": p, "q": p,
                            "effect_definition": method})
        if self.alter:
            res = self.alter(method, out.name.split("__", 1)[1], res)
        out.mkdir(parents=True, exist_ok=True)
        if method == "ancombc2":
            assert arg("--prv-cut") == "0" and arg("--lib-cut") == "0"
            wide = pd.DataFrame({"feature_id": res["feature_id"], "label": "x", "lfc_(Intercept)": 0.0,
                                 "lfc_conditionmacho": res["effect"] * np.log(2), "se_conditionmacho": 0.3,
                                 "W_conditionmacho": res["stat"], "p_conditionmacho": res["p"],
                                 "q_conditionmacho": res["q"], "diff_conditionmacho": res["q"] < 0.05})
            with open(out / "ancombc2_results.tsv", "w") as fo:
                fo.write("# ANCOM-BC2 fake\n")
                wide.to_csv(fo, sep="\t", index=False)
        else:
            res.to_csv(out / "results.tsv", sep="\t", index=False)
        (out / "session_info.txt").write_text("R version 4.5.1\nother attached packages: ANCOMBC_2.14.0 "
                                              "ALDEx2_1.42.0 MicrobiomeStat_1.2\n")
        (out / "warnings.txt").write_text("none\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")


def parse(argv) -> object:
    p = ArgumentParser()
    rb.add_arguments_robustness(p)
    return p.parse_args(argv)


def run_fake(tmp_path, monkeypatch, fake=None, extra=(), low_conf=(), out="rob", em=None):
    if em is None:
        em = em_folder(tmp_path, low_conf=low_conf)
    emdir, feats, meta = em
    fake = fake or FakeR()
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(rb.subprocess, "run", fake)
    args = parse(["--em-uncertainty", str(emdir), "--features", str(feats), "--metadata", str(meta), "--group",
                  "condition", "--reference-level", "ctrl", "-o", str(tmp_path / out), "--threads", "3", *extra])
    assert rb.main_robustness(args) == 0
    return tmp_path / out, fake


def read(path) -> pd.DataFrame:
    return rb._read_tsv(path, ("feature_id", "term", "left_out", "drivers", "methods_significant", "significant_at",
                               "signs_97", "signs_99", "parent"))


# ------------------------------------------------------------------------------------------------ unit tests

def test_prevalence_filter_keeps_features_present_in_enough_samples():
    df = pd.DataFrame({"a": [1, 0, 0], "b": [1, 1, 0], "c": [1, 0, 0], "d": [1, 0, 0]}, index=["x", "y", "z"])
    assert rb.prevalence_filter(df, 0.25) == ["x", "y"]
    assert rb.prevalence_filter(df, 0.5) == ["x"]


def test_stability_counts_untested_features_as_neither_and_uses_successful_fits_only():
    prim = rb._index(pd.DataFrame({"feature_id": ["a", "b"], "term": "g", "effect": [1.0, -2.0], "q": [0.01, 0.01],
                                   "se": 0.2, "effect_definition": "x"}))
    fit = lambda eff, q, ids=("a", "b"): pd.DataFrame({"feature_id": list(ids), "term": "g", "effect": eff, "q": q})
    fits = [fit([0.8, -1.0], [0.01, 0.2]), fit([-0.1, -3.0], [0.5, 0.01]), fit([2.0], [0.001], ids=("a",))]
    out = rb.resample_stability(prim, fits, 0.05)
    assert out.loc[("a", "g"), "n"] == 3
    assert out.loc[("a", "g"), "direction_stability"] == pytest.approx(2 / 3)
    assert out.loc[("a", "g"), "significance_stability"] == pytest.approx(2 / 3)
    # b is not tested in the third fit: it counts as neither same-direction nor significant
    assert out.loc[("b", "g"), "direction_stability"] == pytest.approx(2 / 3)
    assert out.loc[("b", "g"), "significance_stability"] == pytest.approx(1 / 3)
    assert out.loc[("b", "g"), "effect_median"] == pytest.approx(-2.0)
    none = rb.resample_stability(prim, [], 0.05)
    assert (none["n"] == 0).all() and none["direction_stability"].isna().all()


@pytest.mark.parametrize("signs,sig,n,expected", [
    ({"a": 1, "b": 1, "c": 1}, {"a", "b", "c"}, 3, "all_significant_same_direction"),
    ({"a": 1, "b": 1, "c": -1}, {"a", "b"}, 3, "majority_significant_same_direction"),
    ({"a": 1, "b": 1}, {"a", "b"}, 3, "majority_significant_same_direction"),    # third method did not test it
    ({"a": 1, "b": -1, "c": 1}, {"a"}, 3, "single_method"),
    ({"a": 1, "b": -1, "c": 1}, {"a", "b"}, 3, "conflicting_direction"),
    ({"a": 1, "b": 1, "c": 1}, set(), 3, "none"),
])
def test_concordance_classes(signs, sig, n, expected):
    assert rb.classify_concordance(signs, sig, n) == expected


def test_concordance_table_reports_signs_and_rank_correlation_without_pooling():
    tab = lambda e, q: pd.DataFrame({"feature_id": ["a", "b", "c", "d"], "term": "g", "effect": e, "q": q})
    prim = {"ancombc2": tab([1.0, 0.5, -0.2, -1.0], [0.01, 0.5, 0.5, 0.01]),
            "aldex2": tab([2.0, 1.0, -0.4, -2.0], [0.01, 0.5, 0.5, 0.5]),
            "linda": tab([1.5, 0.2, -0.1, 1.0], [0.01, 0.5, 0.5, 0.01])}
    t, header = rb.concordance_table(prim, list(rb.METHODS), 0.05, {})
    t = t.set_index("feature_id")
    assert t.loc["a", "concordance"] == "all_significant_same_direction" and bool(t.loc["a", "direction_agreement"])
    assert t.loc["d", "concordance"] == "conflicting_direction" and t.loc["d", "methods_significant"] == "ancombc2,linda"
    assert t.loc["b", "concordance"] == "none"
    assert "ancombc2 vs aldex2: 1.000" in header[0] and not any("p" == c for c in t.columns)


def test_single_sample_driven_names_the_sample_that_carries_the_call():
    prim = pd.DataFrame({"feature_id": ["a", "b"], "term": "g", "effect": [2.0, 1.0], "se": [0.5, 0.4],
                         "q": [0.001, 0.001]})
    loo = []
    for s in ("S1", "S2", "S3"):
        q = [0.4 if s == "S2" else 0.001, 0.001]
        loo.append((s, pd.DataFrame({"feature_id": ["a", "b"], "term": "g", "effect": [1.5, 1.1], "q": q})))
    out = rb.replicate_rows("linda", prim, loo, 0.05, {"a": "A"})
    s = out[out["left_out"] == "SUMMARY"].set_index("feature_id")
    assert bool(s.loc["a", "single_sample_driven"]) and s.loc["a", "drivers"] == "S2"
    assert s.loc["a", "loo_significance_stability"] == pytest.approx(2 / 3)
    assert s.loc["a", "loo_direction_stability"] == 1.0
    assert not bool(s.loc["b", "single_sample_driven"]) and s.loc["b", "drivers"] == ""
    assert s.loc["a", "mde_80"] == pytest.approx((1.96 + 0.8416) * 0.5)
    assert len(out[out["left_out"] != "SUMMARY"]) == 6


def test_max_single_sample_share_is_the_largest_share_of_a_group():
    ab = pd.DataFrame({"s1": [9.0, 1], "s2": [1.0, 1], "s3": [2.0, 1], "s4": [2.0, 1]}, index=["a", "b"])
    g = pd.Series(["A", "A", "B", "B"], index=["s1", "s2", "s3", "s4"])
    out = rb.max_single_sample_share(ab, g)
    assert out["a"] == pytest.approx(0.9) and out["b"] == pytest.approx(0.5)


GOOD = {"boot_direction_stability": 1.0, "real_direction_stability": 0.95, "boot_significance_stability": 0.9,
        "real_significance_stability": 0.85}


def test_conclusion_rule():
    three = {"ancombc2": 1.0, "aldex2": 1.0, "linda": 1.0}
    stab = {m: GOOD for m in three}
    assert rb.conclude(three, stab, "high", "stable", True)[0] == "high_confidence"
    # significant in all methods, but the reads cannot be assigned reliably
    c, why = rb.conclude(three, stab, "low", "stable", True)
    assert c == "exploratory" and "assignment_confidence low" in why
    assert rb.conclude(three, stab, "high", "unstable", True)[0] == "exploratory"
    assert rb.conclude(three, stab, "high", "stable", False)[0] == "exploratory"
    assert rb.conclude({"linda": 1.0}, stab, "high", "stable", True)[0] == "exploratory"
    assert rb.conclude({"linda": 1.0, "aldex2": -1.0}, stab, "high", "stable", True)[0] == "exploratory"
    shaky = {**stab, "aldex2": {**GOOD, "boot_direction_stability": 0.7}}
    assert rb.conclude(three, shaky, "high", "stable", True)[0] == "exploratory"
    unassessed = {m: {**GOOD, "real_direction_stability": np.nan} for m in three}
    assert rb.conclude(three, unassessed, "high", "stable", True)[0] == "exploratory"
    assert rb.conclude({}, stab, "high", "stable", True)[0] == "not_significant"


def test_cv_classes():
    assert [rb.classify_cv(x) for x in (0.1, 0.2, 0.4, 0.9, np.nan)] == \
        ["stable", "stable", "moderate", "unstable", "not_assessed"]


def test_ancombc2_table_becomes_the_standard_table():
    wide = pd.DataFrame({"feature_id": ["01", "b"], "lfc_(Intercept)": [0, 0], "lfc_gB": [0.7, -0.2],
                         "se_gB": [0.1, 0.2], "W_gB": [7, -1], "p_gB": [0.001, 0.3], "q_gB": [0.002, 0.3]})
    std = rb.ancombc2_to_standard(wide)
    assert list(std.columns) == list(rb.STANDARD_COLUMNS) and list(std["term"]) == ["gB", "gB"]
    assert std.loc[0, "stat"] == 7 and std.loc[0, "ci_low"] == pytest.approx(0.7 - 1.959964 * 0.1)
    assert "natural log" in std.loc[0, "effect_definition"]


# ------------------------------------------------------------------------------------------------ pipeline (fake R)

def test_every_method_gets_the_same_prevalence_filtered_features(tmp_path, monkeypatch):
    out, fake = run_fake(tmp_path, monkeypatch, extra=["--bootstrap-fits", "2", "--realization-fits", "2"])
    by = {}
    for c in fake.calls:
        by.setdefault(c["analysis"], {})[c["method"]] = c["features"]
    assert set(by["primary"]) == set(rb.METHODS)
    for analysis, per in by.items():
        assert len({tuple(v) for v in per.values()}) == 1, analysis
    assert RARE not in by["primary"]["linda"] and UP in by["primary"]["linda"]
    used = pd.read_csv(out / "features.used.tsv", sep="\t", comment="#")
    assert list(used["feature_id"]) == by["primary"]["aldex2"]
    # 12 leave-one-sample-out analyses, 2 + 2 resampled, unique, 2 prevalence cuts
    assert sum(a.startswith("loo_") for a in by) == 12 and {"boot_001", "boot_002", "real_002", "unique",
                                                            "prv_0.25", "prv_0.5"} <= set(by)


def test_failed_fit_is_recorded_and_left_out_of_the_denominator(tmp_path, monkeypatch):
    fake = FakeR(fail={("linda", "boot_002")})
    out, _ = run_fake(tmp_path, monkeypatch, fake, extra=["--bootstrap-fits", "3", "--realization-fits", "2"])
    st = read(out / "fit_status.tsv")
    row = st[(st["method"] == "linda") & (st["analysis"] == "boot_002")].iloc[0]
    assert row["status"] == "failed" and "fake failure" in row["message"]
    stab = read(out / "differential_abundance_stability.tsv")
    assert set(stab.loc[stab["method"] == "linda", "boot_n"]) == {2}
    assert set(stab.loc[stab["method"] == "aldex2", "boot_n"]) == {3}
    up = stab[(stab["feature_id"] == UP)].set_index("method")
    assert (up["boot_direction_stability"] == 1.0).all() and (up["primary_significant"]).all()
    assert not (out / "fits/linda__boot_002").exists()
    report = (out / "robustness_report.md").read_text()
    assert "## Failed fits" in report and "boot_002" in report


def test_significant_everywhere_but_low_assignment_confidence_is_exploratory(tmp_path, monkeypatch):
    out, _ = run_fake(tmp_path, monkeypatch, low_conf=(UP,), extra=["--bootstrap-fits", "3"])
    summ = read(out / "robustness_summary.tsv").set_index("feature_id")
    assert summ.loc[UP, "cross_method_agreement"] == "all_significant_same_direction"
    assert summ.loc[UP, "assignment_confidence"] == "low"
    assert summ.loc[UP, "conclusion"] == "exploratory" and "assignment_confidence low" in summ.loc[UP, "unmet_conditions"]
    assert "method agreement" in summ.loc[UP, "robustness_basis"]
    assert "read-assignment stability" not in summ.loc[UP, "robustness_basis"]
    assert (summ.drop(UP)["conclusion"] == "not_significant").all()


def test_high_confidence_needs_every_dimension_and_report_wording(tmp_path, monkeypatch):
    out, _ = run_fake(tmp_path, monkeypatch, extra=["--bootstrap-fits", "3"])
    summ = read(out / "robustness_summary.tsv").set_index("feature_id")
    assert summ.loc[UP, "conclusion"] == "high_confidence"
    assert summ.loc[UP, "robustness_basis"].split(";") == [
        "method agreement", "effect-direction stability", "read-assignment stability",
        "biological replicate consistency"]
    summary_text = (out / "robustness_summary.tsv").read_text()
    assert rb.CONCLUSION_RULE in summary_text
    report = (out / "robustness_report.md").read_text()
    assert "Robust across three complementary differential abundance frameworks" in report
    for f in out.rglob("*"):
        if f.is_file() and f.suffix in (".md", ".tsv"):
            assert "validated by" not in f.read_text().lower(), f
    assert "## Low replication" not in report             # 6 per group
    conc = (out / "cross_method_concordance.tsv").read_text()
    assert "not independent validation" in conc and "Spearman" in conc


def test_single_sample_driven_call_is_not_high_confidence(tmp_path, monkeypatch):
    def alter(method, analysis, res):           # without sample B3 the call is lost in every method
        if analysis == "loo_B3":
            res.loc[res["feature_id"] == UP, ["p", "q"]] = 0.3
        return res
    out, _ = run_fake(tmp_path, monkeypatch, FakeR(alter=alter), extra=["--bootstrap-fits", "2"])
    rep = read(out / "replicate_sensitivity.tsv")
    s = rep[(rep["left_out"] == "SUMMARY") & (rep["feature_id"] == UP)]
    assert s["single_sample_driven"].astype(str).eq("True").all() and set(s["drivers"]) == {"B3"}
    assert set(s["n_per_group"]) == {"ctrl=6;macho=6"}
    summ = read(out / "robustness_summary.tsv").set_index("feature_id")
    assert summ.loc[UP, "conclusion"] == "exploratory"
    assert str(summ.loc[UP, "replicate_consistency"]) == "False"


def test_leave_one_unit_out_drops_every_sample_of_the_unit(tmp_path, monkeypatch):
    out, fake = run_fake(tmp_path, monkeypatch, extra=["--unit", "plant", "--bootstrap-fits", "1",
                                                       "--realization-fits", "1", "--prv-cuts", ""])
    loo = {c["analysis"] for c in fake.calls if c["analysis"].startswith("loo_")}
    assert loo == {f"loo_p{i}" for i in range(6)}
    kept = [c for c in fake.calls if c["analysis"] == "loo_p0"][0]["samples"]
    assert len(kept) == 10 and not {"A1", "A2"} & set(kept)


def test_manifest_seals_every_output_and_a_rerun_removes_stale_fits(tmp_path, monkeypatch):
    em = em_folder(tmp_path)
    out, _ = run_fake(tmp_path, monkeypatch, extra=["--bootstrap-fits", "3"], em=em)
    rec = manifest.verify(out)
    arts = rec["artifacts"]
    for f in ("differential_abundance_stability.tsv", "cross_method_concordance.tsv", "replicate_sensitivity.tsv",
              "robustness_summary.tsv", "robustness_report.md", "features.used.tsv", "fit_status.tsv",
              "figures/bootstrap_effects.png", "figures/direction_stability_heatmap.png", "figures/concordance.png",
              "figures/replicate_abundance.png", "fits/aldex2__boot_003/results.tsv",
              "fits/linda__primary/session_info.txt"):
        assert f in arts, f
    assert rec["step"] == "da-robustness" and rec["seed"] == 1 and rec["em_uncertainty"]["seed"] == 7
    assert "ALDEx2 1.42.0" in rec["versions"]["R (aldex2)"]
    assert not list(out.glob(".staging-*"))
    out, _ = run_fake(tmp_path, monkeypatch, extra=["--bootstrap-fits", "1"], em=em)
    assert not (out / "fits/aldex2__boot_003").exists() and (out / "fits/aldex2__boot_001/results.tsv").exists()
    manifest.verify(out)


def test_refuses_a_changed_em_folder(tmp_path, monkeypatch):
    em = em_folder(tmp_path)
    with open(em[0] / "bootstrap/counts_002.tsv", "a") as fo:
        fo.write("f99\t" + "\t".join(["1"] * 12) + "\n")
    with pytest.raises(Exception, match="changed after"):
        run_fake(tmp_path, monkeypatch, em=em)


def test_da_resolutions_compares_parents_across_runs(tmp_path, monkeypatch):
    em = em_folder(tmp_path)
    a, _ = run_fake(tmp_path, monkeypatch, extra=["--bootstrap-fits", "1"], em=em, out="r97")

    def flip(method, analysis, res):            # at the second resolution the change is not detected
        res.loc[res["feature_id"] == UP, ["p", "q"]] = 0.5
        return res
    b, _ = run_fake(tmp_path, monkeypatch, FakeR(alter=flip), extra=["--bootstrap-fits", "1"], em=em, out="r99")
    p = ArgumentParser()
    rb.add_arguments_resolutions(p)
    args = p.parse_args(["--run", f"97={a}", "--run", f"99={b}", "-o", str(tmp_path / "res")])
    assert rb.main_resolutions(args) == 0
    t = read(tmp_path / "res/multi_resolution_robustness.tsv")
    t = t.rename(columns=str)
    row = t[t["parent"] == "bin5"].iloc[0]                      # f05 and f15, f25 share bin5
    assert row["significant_at"] == "97" and row["resolution_stability"] == "resolution_dependent"
    assert row["n_features_97"] == 3 and row["signs_97"] == "+"
    assert set(t.loc[t["parent"] != "bin5", "resolution_stability"]) == {"not_significant"}
    manifest.verify(tmp_path / "res")


def test_resolution_table_stable_when_significant_everywhere_in_one_direction():
    conc = pd.DataFrame({"feature_id": ["x1"], "term": "g", "n_methods_significant": [2],
                         "methods_significant": ["linda,aldex2"], "sign_linda": [1.0], "sign_aldex2": [1.0],
                         "concordance": ["majority_significant_same_direction"]})
    used = pd.DataFrame({"feature_id": ["x1"], "parent_feature": [""]})
    t = rb.resolution_table({"97": (conc, used), "98": (conc, used)})
    assert t.loc[0, "parent"] == "x1" and t.loc[0, "resolution_stability"] == "stable"
    assert bool(t.loc[0, "consistent_direction"])
    neg = conc.assign(sign_linda=-1.0, sign_aldex2=-1.0)
    t = rb.resolution_table({"97": (conc, used), "98": (neg, used)})
    assert t.loc[0, "resolution_stability"] == "resolution_dependent" and not bool(t.loc[0, "consistent_direction"])


# ------------------------------------------------------------------------------------------------ real R

def _rscript():
    r = os.environ.get("UNCLASSIFIED_META_RSCRIPT") or shutil.which("Rscript")
    return r if r and Path(r).exists() else None


def _has_packages(r) -> bool:
    expr = 'q(status = if (all(vapply(c("ANCOMBC", "ALDEx2", "MicrobiomeStat"), requireNamespace, logical(1), ' \
           'quietly = TRUE))) 0 else 1)'
    return subprocess.run([r, "-e", expr], capture_output=True).returncode == 0


@pytest.mark.skipif(not _rscript() or not _has_packages(_rscript()),
                    reason="Rscript with ANCOMBC, ALDEx2 and MicrobiomeStat not available")
def test_three_r_methods_on_synthetic_data(tmp_path):
    # kept small for time: R start-up alone is ~13 s here, so 18 fits (3 methods x primary, unique, 1 bootstrap,
    # 1 realization, leave-one-batch-out x 2) on 12 threads
    em, feats, meta = em_folder(tmp_path, n_boot=1, n_real=1)
    args = parse(["--em-uncertainty", str(em), "--features", str(feats), "--metadata", str(meta), "--group",
                  "condition", "--reference-level", "ctrl", "-o", str(tmp_path / "rob"), "--threads", "12",
                  "--bootstrap-fits", "1", "--realization-fits", "1", "--mc-samples", "16", "--prv-cuts", "",
                  "--unit", "batch", "--rscript", _rscript(), "--p-adjust", "BH"])
    assert rb.main_robustness(args) == 0
    out = tmp_path / "rob"
    st = read(out / "fit_status.tsv")
    assert len(st) == 18 and (st["status"] == "ok").all(), st[st["status"] != "ok"].to_string()
    prim = {m: read(out / f"fits/{m}__primary/results.tsv") for m in rb.METHODS}
    tested = {m: set(r["feature_id"]) for m, r in prim.items()}
    assert tested["ancombc2"] == tested["aldex2"] == tested["linda"] and RARE not in tested["linda"]
    for m, r in prim.items():
        assert set(r["term"]) == {"conditionmacho"}, m              # model.matrix column name in every method
        assert list(r.columns) == list(rb.STANDARD_COLUMNS)
        up = r.set_index("feature_id").loc[UP]
        assert up["effect"] > 1 and up["q"] < 0.05, (m, up.to_dict())    # macho vs ctrl: up in every method
    conc = read(out / "cross_method_concordance.tsv").set_index("feature_id")
    assert conc.loc[UP, "concordance"] == "all_significant_same_direction"
    assert {conc.loc[UP, f"sign_{m}"] for m in rb.METHODS} == {1.0}
    assert (conc.drop(UP)["n_methods_significant"] == 0).all()
    rho = [ln for ln in (out / "cross_method_concordance.tsv").read_text().splitlines() if "Spearman" in ln]
    assert len(rho) == 3 and all(float(ln.split(": ")[1].split()[0]) > 0.5 for ln in rho), rho
    stab = read(out / "differential_abundance_stability.tsv")
    up = stab[stab["feature_id"] == UP].set_index("method")
    assert (up["boot_direction_stability"] == 1).all() and (up["real_direction_stability"] == 1).all()
    assert "rounded" in (out / "fits/aldex2__primary/warnings.txt").read_text()
    versions = manifest.verify(out)["versions"]
    assert "ALDEx2 1.42.0" in versions["R (aldex2)"] and "MicrobiomeStat 1.4" in versions["R (linda)"]
