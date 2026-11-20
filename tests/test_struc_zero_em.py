"""Structural zeros with expected (fractional) counts, with multiple imputation, and the detection evidence."""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import da, detection, manifest, taxa
from unclassified_meta.cli import build_parser, main

from test_em_merge import quant
from test_em_taxa import em_file
from test_ported_units import CFR_HEADER, cfr, taxonomy_files
from test_review_fixes import counts_with_a_structural_zero, da_cli, rscript


def read_sz(path):
    return da.read_structural_zeros(path).set_index("feature_id")


def run_da(tmp_path):
    return da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv",
                  "--metadata", tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A",
                  "--struc-zero", "-o", tmp_path / "d", "--rscript", rscript())


@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_a_fraction_of_a_shared_read_is_presence_under_ancombc2s_own_rule_and_the_evidence_is_reported(tmp_path):
    counts_with_a_structural_zero(tmp_path)                       # t01 is absent from every B sample
    c = pd.read_csv(tmp_path / "counts.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id").astype(float)
    c.columns = c.columns.astype(str)
    c.loc["t01", "101"] = 0.33                                    # EM: a third of a shared read in one B sample
    c.to_csv(tmp_path / "counts.tsv", sep="\t")
    rc, log = run_da(tmp_path)
    assert rc == 0, log
    sz = read_sz(tmp_path / "d" / "structural_zeros.tsv")
    assert not sz.loc["t01", "flagged_by_ancombc2"] and not sz.loc["t01", "absent_in_a_group"]
    assert int(sz.loc["t01", "present_samples (group = B)"]) == 1 and int(sz.loc["t01", "samples (group = B)"]) == 6
    assert float(sz.loc["t01", "reads (group = B)"]) == pytest.approx(0.33)
    assert float(sz.loc["t01", "max_count (group = B)"]) == pytest.approx(0.33)
    assert int(sz.loc["t01", "present_samples (group = A)"]) == 6
    res = da.read_results(tmp_path / "d" / "ancombc2_results.tsv")
    assert "t01" in set(res["feature_id"].astype(str))            # tested: no value was changed or cut
    header = (tmp_path / "d" / "ancombc2_results.tsv").read_text()
    assert "count != 0" in header and "unclassified-meta detection" in header and "--neg-lb" in header


@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_a_group_without_any_read_in_a_whole_read_table_is_flagged_and_left_untested(tmp_path):
    counts_with_a_structural_zero(tmp_path)                       # integer counts, as from --assign lca
    rc, log = run_da(tmp_path)
    assert rc == 0, log
    sz = read_sz(tmp_path / "d" / "structural_zeros.tsv")
    assert sz.loc["t01", "flagged_by_ancombc2"] and sz.loc["t01", "absent_in_groups"] == "B"
    assert int(sz.loc["t01", "present_samples (group = B)"]) == 0 and float(sz.loc["t01", "reads (group = B)"]) == 0
    assert float(sz.loc["t01", "max_count (group = B)"]) == 0
    assert not sz.loc["t02", "flagged_by_ancombc2"] and int(sz.loc["t02", "present_samples (group = B)"]) == 6
    res = da.read_results(tmp_path / "d" / "ancombc2_results.tsv")
    assert "t01" not in set(res["feature_id"].astype(str))


def sz_frame(rows):
    return pd.DataFrame(rows, columns=["feature_id", "label", "absent_in_groups", "absent_in_a_group"])


def test_structural_zeros_of_the_imputations_are_summarised_per_feature():
    fits = [sz_frame([("f1", "A one", "B", True), ("f2", "two", "", False), ("f3", "three", "A", True)]),
            sz_frame([("f1", "A one", "B", True), ("f2", "two", "B", True), ("f3", "three", "B", True)]),
            sz_frame([("f1", "A one", "B", True), ("f2", "two", "", False), ("f3", "three", "A", True)])]
    out = da.combine_structural_zeros(fits).set_index("feature_id")
    assert out.loc["f1", "n_fits"] == 3 and out.loc["f1", "n_fits_absent"] == 3 and out.loc["f1", "frac_fits_absent"] == 1
    assert out.loc["f1", "absent_in_groups"] == "B" and out.loc["f1", "absent_in_every_fit"]
    assert not out.loc["f1", "inconsistent"]
    assert out.loc["f2", "n_fits_absent"] == 1 and not out.loc["f2", "absent_in_every_fit"]
    assert out.loc["f2", "inconsistent"] and out.loc["f2", "label"] == "two"
    assert out.loc["f3", "absent_in_groups"] == "A (2); B (1)" and not out.loc["f3", "inconsistent"]
    assert "n_fits_flagged_by_ancombc2" not in out.columns
    assert not any(c.startswith(("p_", "q_", "lfc_", "se_", "W_")) for c in out.columns)   # nothing is pooled


def em_folder(tmp_path, imputations=3):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=imputations)
    return out, nodes, names


def fake_r(calls, struc_zero_fit=1):
    """A stand-in for the R script: t|564|S is a structural zero (B) in the imputation fit `struc_zero_fit` only
    (call 0 is the fit of the expected table)."""
    def run(cmd, capture_output=True, text=True):
        o = Path(cmd[cmd.index("--outdir") + 1])
        f = list(pd.read_csv(cmd[cmd.index("--counts") + 1], sep="\t", dtype=str)["feature_id"])
        absent = len(calls) == struc_zero_fit
        calls.append(cmd)
        tested = [x for x in f if not (absent and x == "t|564|S")]
        pd.DataFrame({"feature_id": tested, "label": "x", "lfc_(Intercept)": 0.0, "se_(Intercept)": 1,
                      "lfc_groupB": 1.0, "se_groupB": 0.2, "diff_groupB": "TRUE", "passed_ss_groupB": "TRUE",
                      "diff_robust_groupB": "TRUE"}).to_csv(o / "ancombc2_results.tsv", sep="\t", index=False)
        if "--struc-zero" in cmd:
            sz_frame([(x, "x", "B" if absent and x == "t|564|S" else "", absent and x == "t|564|S") for x in f]) \
                .assign(flagged_by_ancombc2=lambda d: d["absent_in_a_group"]) \
                .to_csv(o / "structural_zeros.tsv", sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    return run


def da_args(out, meta, outdir, *extra):
    return build_parser().parse_args(["da", "--counts", str(out / "merged/counts.tsv"), "--features",
                                      str(out / "merged/features.tsv"), "--metadata", str(meta), "--group", "group",
                                      "--imputations", "-o", str(outdir), *extra])


def test_da_imputations_keeps_and_summarises_the_structural_zeros_of_every_fit(tmp_path, monkeypatch):
    out, _, _ = em_folder(tmp_path)
    taxa.merge_with_panel(out, quant(tmp_path, [1, 0, 0], 1), None, out / "merged")
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\ns\tA\n")
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(da, "has_ancombc", lambda r: True)
    calls = []
    monkeypatch.setattr(da.subprocess, "run", fake_r(calls))
    da.run_ancombc2(da_args(out, meta, tmp_path / "da", "--struc-zero"))
    rec = manifest.verify(tmp_path / "da")
    for name in ("structural_zeros.tsv", "mi_structural_zeros.tsv", "imputations/structural_zeros_001.tsv",
                 "imputations/structural_zeros_003.tsv"):
        assert name in rec["artifacts"]
    first = (tmp_path / "da/mi_structural_zeros.tsv").read_text().splitlines()[0]
    assert first.startswith("#") and "not pooled into any p-value" in first
    s = da.read_results(tmp_path / "da/mi_structural_zeros.tsv").set_index("feature_id")
    assert s.loc["t|564|S", "n_fits"] == 3 and s.loc["t|564|S", "n_fits_absent"] == 1
    assert bool(s.loc["t|564|S", "inconsistent"]) and not bool(s.loc["t|564|S", "absent_in_every_fit"])
    assert s.loc["t|564|S", "n_fits_flagged_by_ancombc2"] == 1 and s.loc["t|562|S", "n_fits_absent"] == 0
    # a feature left untested as a structural zero in one fit is not pooled: p = 1, no combined estimate
    mi = da.read_results(tmp_path / "da/mi_results.tsv").set_index("feature_id")
    assert mi.loc["t|564|S", "p_groupB"] == 1.0 and np.isnan(mi.loc["t|564|S", "lfc_groupB"])
    assert mi.loc["t|564|S", "n_fits_groupB"] == 2 and not bool(mi.loc["t|564|S", "diff_groupB"])
    assert mi.loc["t|562|S", "lfc_groupB"] == pytest.approx(1.0)

    calls.clear()                        # a later run without --struc-zero removes them through the manifest
    da.run_ancombc2(da_args(out, meta, tmp_path / "da"))
    rec = manifest.verify(tmp_path / "da")
    gone = ("structural_zeros.tsv", "mi_structural_zeros.tsv", "imputations/structural_zeros_001.tsv")
    assert not any(n in rec["artifacts"] or (tmp_path / "da" / n).exists() for n in gone)
    assert "mi_results.tsv" in rec["artifacts"]


def test_em_taxa_table_writes_and_seals_the_unique_read_counts(tmp_path):
    out, nodes, names = em_folder(tmp_path, imputations=0)
    rec = manifest.verify(out)
    assert "unique_counts.tsv" in rec["artifacts"]
    u = pd.read_csv(out / "unique_counts.tsv", sep="\t", dtype=str).set_index("feature_id")
    c = pd.read_csv(out / "counts.tsv", sep="\t", dtype=str).set_index("feature_id")
    assert list(u.index) == list(c.index) and list(u.columns) == list(c.columns)
    assert u["s"].str.fullmatch(r"[0-9]+").all()                  # whole reads
    # 9 reads only E. coli, 1 only E. fergusonii; the 10 shared reads and the phage-or-Salmonella reads are not unique
    assert u.loc["t|562|S", "s"] == "9" and u.loc["t|564|S", "s"] == "1"
    assert float(c.loc["t|562|S", "s"]) > 9 and float(c.loc["t|564|S", "s"]) > 1
    # an LCA table has no unique counts; the stale file of the earlier EM generation goes through the manifest
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0)
    assert "unique_counts.tsv" not in manifest.verify(out)["artifacts"] and not (out / "unique_counts.tsv").exists()


def two_group_files(tmp_path):
    """A: E. coli and E. fergusonii each have reads of their own. B: E. coli has its own reads; E. fergusonii only
    shares one read of s3 with Salmonella (no read of its own in B)."""
    d = tmp_path / "reads"
    d.mkdir()
    files = {}
    for s in ("s1", "s2", "s3", "s4"):
        body = CFR_HEADER + "".join(cfr(f"{s}c{i}", 562, 400, 60) for i in range(6))
        if s in ("s1", "s2"):
            body += "".join(cfr(f"{s}f{i}", 564, 400, 60) for i in range(3))
            body += "".join(cfr(f"{s}m{i}", 562, 400, 60) + cfr(f"{s}m{i}", 564, 400, 60) for i in range(4))
        if s == "s3":
            body += cfr("s3x", 564, 400, 60) + cfr("s3x", 28901, 400, 60)
        files[s] = d / f"{s}.cfr.tsv"
        files[s].write_text(body)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tcondition\ns1\tA\ns2\tA\ns3\tB\ns4\tB\n")
    return files, meta


def test_detection_shows_that_a_group_holds_only_shared_reads_and_seals_its_tables(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    files, meta = two_group_files(tmp_path)
    taxa.build_table_from_per_read(files, tmp_path / "em", nodes, names, "S", (), "t", 0, 0, assign="em",
                                   imputations=20, seed=4)
    taxa.build_table_from_per_read(files, tmp_path / "lca", nodes, names, "S", (), "t", 0, 0)
    rc = main(["detection", "--em", str(tmp_path / "em"), "--lca", str(tmp_path / "lca"), "--metadata", str(meta),
               "--group", "condition", "-o", str(tmp_path / "det")])
    assert rc == 0
    rec = manifest.verify(tmp_path / "det")
    assert set(detection.OUTPUTS) <= set(rec["artifacts"]) and rec["imputations"] == 20
    g = da.read_results(tmp_path / "det/detection_by_group.tsv").set_index(["feature_id", "level"])
    fer = g.loc[("t|564|S", "B")]
    assert bool(fer["struc_zero_unique"]) and not bool(fer["struc_zero_em_expected"]) and bool(fer["struc_zero_lca"])
    assert fer["em_expected_present_samples"] == 1 and fer["em_expected_reads"] == pytest.approx(0.5)
    assert 0 < fer["imputation_all_zero_frac"] < 1                # the shared read is drawn to it in some imputations
    assert g.loc[("t|564|S", "A"), "imputation_all_zero_frac"] == 0
    coli = g.loc[("t|562|S", "B")]
    assert not any(bool(coli[f"struc_zero_{k}"]) for k in ("lca", "unique", "em_expected"))
    summ = da.read_results(tmp_path / "det/detection_summary.tsv").set_index("feature_id")
    assert summ.loc["t|564|S", "interpretation"].startswith(
        "no unique-read support in B; EM expected count there comes only from shared reads")
    assert summ.loc["t|564|S", "struc_zero_levels_unique"] == "B" and bool(summ.loc["t|564|S", "inconsistent"])
    assert summ.loc["t|564|S", "imputation_all_zero_frac (condition = B)"] == fer["imputation_all_zero_frac"]
    assert summ.loc["t|562|S", "interpretation"] == "unique-read support in every level"
    assert pd.isna(summ.loc["t|562|S", "struc_zero_levels_unique"]) and not bool(summ.loc["t|562|S", "inconsistent"])
    s = da.read_results(tmp_path / "det/detection_by_sample.tsv").set_index(["feature_id", "sample"])
    assert s.loc[("t|564|S", "s3"), "unique_count"] == 0 and s.loc[("t|564|S", "s3"), "em_expected"] == 0.5
    assert 0 < s.loc[("t|564|S", "s3"), "imputation_nonzero_frac"] < 1 and s.loc[("t|564|S", "s3"), "lca_count"] == 0


def test_detection_reports_only_the_features_asked_for_by_id_or_label(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    files, meta = two_group_files(tmp_path)
    taxa.build_table_from_per_read(files, tmp_path / "em", nodes, names, "S", (), "t", 0, 0, assign="em")
    detection.run(tmp_path / "em", meta, "condition", tmp_path / "det",
                  features=["Escherichia fergusonii", "t|562|S"])
    summ = da.read_results(tmp_path / "det/detection_summary.tsv")
    assert list(summ["feature_id"]) == ["t|564|S", "t|562|S"]
    assert "inconsistent" not in summ.columns and "struc_zero_levels_lca" not in summ.columns   # no such source
    with pytest.raises(ValueError, match="no feature"):
        detection.run(tmp_path / "em", meta, "condition", tmp_path / "det2", features=["Homo sapiens"])


def test_detection_refuses_a_table_made_without_em(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    files, meta = two_group_files(tmp_path)
    taxa.build_table_from_per_read(files, tmp_path / "lca", nodes, names, "S", (), "t", 0, 0)
    with pytest.raises(ValueError, match="--assign em"):
        detection.run(tmp_path / "lca", meta, "condition", tmp_path / "det")
