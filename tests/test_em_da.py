"""da --imputations: per-imputation fits combined with Rubin's rules."""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import da, manifest, taxa

from test_em_merge import quant
from test_em_taxa import em_file
from test_ported_units import taxonomy_files


def test_p_adjust_matches_r():
    p = [0.01, 0.04, 0.03, 0.2, np.nan]
    assert np.allclose(da.p_adjust(p, "BH")[:4], [0.04, 0.0533333, 0.0533333, 0.2], atol=1e-6)
    assert np.allclose(da.p_adjust([0.01, 0.02, 0.03, 0.04, 0.05], "holm"), [0.05, 0.08, 0.09, 0.09, 0.09])
    assert np.isnan(da.p_adjust(p, "holm")[4])
    with pytest.raises(Exception, match="not available"):
        da.p_adjust(p, "hommel")


def result(lfc, se, ss=("TRUE", "TRUE"), robust=("TRUE", "FALSE"), ids=("f1", "f2")):
    return pd.DataFrame({"feature_id": list(ids), "label": ["A", "B"][:len(ids)], "lfc_(Intercept)": 0.0,
                         "se_(Intercept)": 1.0, "lfc_groupB": lfc, "se_groupB": se,
                         "diff_groupB": ["TRUE", "FALSE"][:len(ids)], "passed_ss_groupB": list(ss)[:len(ids)],
                         "diff_robust_groupB": list(robust)[:len(ids)]})


def test_robust_needs_robust_calls_in_the_fits_not_only_stable_non_significance():
    fits = [result([2.0, 0.1], [0.5, 0.5], robust=("FALSE", "FALSE")) for _ in range(4)]
    out, _ = da.combine_imputations(fits, 0.05, "holm")
    assert bool(out.loc[0, "diff_groupB"]) and not bool(out.loc[0, "diff_robust_groupB"])
    assert out.loc[0, "passed_ss_frac_groupB"] == 1.0 and out.loc[0, "robust_frac_groupB"] == 0.0


def test_a_feature_missing_from_one_fit_is_reported_with_p_1_and_counts_in_the_adjustment():
    fits = [result([2.0, 2.0], [0.5, 0.5]) for _ in range(3)] + [result([2.0], [0.5], ids=("f1",))]
    out, _ = da.combine_imputations(fits, 0.05, "bonferroni")
    row = out.set_index("feature_id")
    assert list(row.index) == ["f1", "f2"] and row.loc["f2", "p_groupB"] == 1.0
    assert row.loc["f2", "n_fits_groupB"] == 3 and not bool(row.loc["f2", "diff_groupB"])
    assert row.loc["f1", "q_groupB"] == pytest.approx(min(1, 2 * row.loc["f1", "p_groupB"]))


def test_a_fit_with_p_1_and_nonzero_w_is_a_failed_test_not_pooled_as_significant():
    # ANCOM-BC2 reports p = 1 when the feature's degrees of freedom are NA, even when W is large (ISSUE #0)
    fits = []
    for _ in range(4):
        r = result([5.0, 2.0], [0.2, 0.5])
        r["W_groupB"], r["p_groupB"] = [25.0, 4.0], [1.0, 0.0001]
        fits.append(r)
    out, _ = da.combine_imputations(fits, 0.05, "holm")
    row = out.set_index("feature_id")
    assert row.loc["f1", "p_groupB"] == 1.0 and row.loc["f1", "n_fits_groupB"] == 0
    assert not bool(row.loc["f1", "diff_groupB"])
    assert row.loc["f2", "p_groupB"] < 0.05 and row.loc["f2", "n_fits_groupB"] == 4


def test_combined_fit_widens_the_error_when_the_imputations_disagree():
    same = [result([2.0, 0.1], [0.5, 0.5]) for _ in range(5)]
    vary = [result([2.0 + d, 0.1], [0.5, 0.5]) for d in (-1.5, -0.7, 0, 0.7, 1.5)]
    a, terms = da.combine_imputations(same, 0.05, "holm")
    b, _ = da.combine_imputations(vary, 0.05, "holm")
    assert terms == ["groupB"]
    assert a.loc[0, "se_groupB"] == pytest.approx(0.5) and b.loc[0, "se_groupB"] > 0.5
    assert a.loc[0, "p_groupB"] < b.loc[0, "p_groupB"]
    assert a.loc[0, "fmi_groupB"] == 0 and b.loc[0, "fmi_groupB"] > 0.5
    assert bool(a.loc[0, "diff_robust_groupB"]) and a.loc[0, "diff_frac_groupB"] == 1.0


def test_da_imputations_fits_every_sealed_table_and_seals_the_combined_result(tmp_path, monkeypatch):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=3)
    taxa.merge_with_panel(out, quant(tmp_path, [1, 0, 0], 1), None, out / "merged")
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\ns\tA\n")
    seen = []

    def run(cmd, capture_output=True, text=True):
        o = Path(cmd[cmd.index("--outdir") + 1])
        counts = pd.read_csv(cmd[cmd.index("--counts") + 1], sep="\t")
        seen.append(counts["s"].sum())
        lfc = 1.0 + 0.1 * len(seen)
        pd.DataFrame({"feature_id": counts["feature_id"], "label": "x", "lfc_(Intercept)": 0.0, "se_(Intercept)": 1,
                      "lfc_groupB": lfc, "se_groupB": 0.2, "diff_groupB": "TRUE", "passed_ss_groupB": "TRUE",
                      "diff_robust_groupB": "TRUE"}) \
            .to_csv(o / "ancombc2_results.tsv", sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(da, "has_ancombc", lambda r: True)
    monkeypatch.setattr(da.subprocess, "run", run)
    from unclassified_meta.cli import build_parser
    args = build_parser().parse_args(["da", "--counts", str(out / "merged/counts.tsv"), "--features",
                                      str(out / "merged/features.tsv"), "--metadata", str(meta), "--group", "group",
                                      "--imputations", "-o", str(tmp_path / "da")])
    da.run_ancombc2(args)
    assert len(seen) == 4                                   # the expected table + 3 imputations
    rec = manifest.verify(tmp_path / "da")
    assert rec["imputations"] == 3 and "mi_results.tsv" in rec["artifacts"]
    mi = pd.read_csv(tmp_path / "da/mi_results.tsv", sep="\t", comment="#")
    assert mi["lfc_groupB"].iloc[0] == pytest.approx(np.mean([1.2, 1.3, 1.4]))
    assert (mi["se_groupB"] > 0.2).all() and (mi["fmi_groupB"] > 0).all()


def test_da_imputations_needs_imputation_tables(tmp_path, monkeypatch):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em")
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\ns\tA\n")

    def run(cmd, capture_output=True, text=True):
        (Path(cmd[cmd.index("--outdir") + 1]) / "ancombc2_results.tsv").write_text("feature_id\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(da, "has_ancombc", lambda r: True)
    monkeypatch.setattr(da.subprocess, "run", run)
    from unclassified_meta.cli import build_parser
    args = build_parser().parse_args(["da", "--counts", str(out / "counts.tsv"), "--features",
                                      str(out / "features.tsv"), "--metadata", str(meta), "--group", "group",
                                      "--imputations", "-o", str(tmp_path / "da")])
    with pytest.raises(Exception, match="holds 0 imputation"):
        da.run_ancombc2(args)


def test_a_fit_ancombc2_could_not_test_stays_untested_after_pooling():
    # ANCOM-BC2 reports p = 1 with a large W when a feature's degrees of freedom are NA (p_hat[is.na(p_hat)] = 1)
    fits = [result([4.3, 0.0], [0.2, 0.5]).assign(**{"W_groupB": [21.5, 0.0], "p_groupB": [1.0, 1.0]})
            for _ in range(4)]
    out, _ = da.combine_imputations(fits, 0.05, "holm")
    row = out.set_index("feature_id")
    assert row.loc["f1", "n_fits_groupB"] == 0 and row.loc["f1", "p_groupB"] == 1.0
    assert not bool(row.loc["f1", "diff_groupB"]) and not bool(row.loc["f1", "diff_robust_groupB"])
    assert row.loc["f2", "n_fits_groupB"] == 4               # W = 0, p = 1 is an ordinary test
