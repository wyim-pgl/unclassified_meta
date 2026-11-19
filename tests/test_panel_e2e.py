"""End-to-end test of the panel workflow on synthetic two-group data: run -> panel-build -> panel-quant -> ANCOM-BC2.

The ANCOM-BC2 step needs the R environment (environment-r.yml); point UNCLASSIFIED_META_RSCRIPT at its Rscript or
put it on PATH. Without it that one test is skipped.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import da

HERE = Path(__file__).parent
TOOLS = ("spades.py", "minimap2", "blastn", "diamond", "bowtie2", "samtools", "seqtk", "seqkit")
pytestmark = [pytest.mark.e2e,
              pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS),
                                 reason="external tools not on PATH (activate the environment)")]
V4 = "refseq_viral_r999|NC_900004.1|like"


def cli(*args, ok=True):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)], capture_output=True, text=True)
    if ok:
        assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr


@pytest.fixture(scope="module")
def work(tmp_path_factory):
    d = tmp_path_factory.mktemp("da")
    subprocess.run([sys.executable, str(HERE / "make_synthetic_da.py"), str(d / "data")], check=True)
    cli("run", "--input-dir", d / "data" / "reads", "--refs", d / "data" / "refs",
        "--host", f"testhost={d / 'data' / 'host.fa'}", "--mode", "rna", "-t", 4, "--memory", 8,
        "--subsample", "all", "-o", d / "run")
    cli("panel-build", "--run", d / "run", "-o", d / "panel", "--categories", "virus,unclassified", "-t", 4)
    cli("panel-quant", "--panel", d / "panel", "--input-dir", d / "data" / "reads", "-o", d / "quant", "-t", 4)
    return d


def test_run_exports_raw_contig_counts(work):
    m = pd.read_csv(work / "run" / "tables" / "contig_counts.tsv", sep="\t").set_index("contig")
    t = pd.read_csv(work / "run" / "tables" / "read_category_by_sample.tsv", sep="\t").set_index("sample")
    assert list(m.columns) == [f"{g}{k}" for g in "AB" for k in range(1, 7)]
    assigned = t[[c for c in t.columns if c.startswith("n_") and c != "n_not_assembled"]].sum(axis=1)
    assert (m.sum(axis=0) == assigned[m.columns]).all()


def test_panel_features_keep_like_bins_apart_and_have_namespaced_ids(work):
    f = pd.read_csv(work / "panel" / "features.tsv", sep="\t")
    cls = f["identity_class"].value_counts().to_dict()
    assert cls == {"same": 3, "like": 3, "reference": 3, "protein_like": 1, "unmatched": 1}
    assert f["feature_id"].is_unique and f["feature_id"].str.count(r"\|").eq(2).all()
    assert set(f.loc[f.identity_class == "like", "accession"]) == {"NC_900004.1", "NC_900005.1", "NC_900006.1"}
    assert f.loc[f.identity_class == "protein_like", "label"].item() == "Test virus 7-like (protein only)"
    s2f = pd.read_csv(work / "panel" / "seq2feature.tsv", sep="\t")
    assert {"decoy_host:testhost", "decoy_rrna"} <= set(s2f["feature_or_decoy"])
    prov = (work / "panel" / "provenance.json").read_text()
    assert '"panel_id"' in prov and '"sha256"' in prov and "RefSeq release 999" in prov


def test_every_fragment_has_exactly_one_outcome_and_counts_match_the_truth(work):
    led = pd.read_csv(work / "quant" / "ledger.tsv", sep="\t").set_index("sample")
    outcome_cols = [c for c in led.columns if c not in ("paired", "input_fragments", "wall_s", "max_rss_gb")]
    assert (led[outcome_cols].sum(axis=1) == led["input_fragments"]).all()
    counts = pd.read_csv(work / "quant" / "counts.tsv", sep="\t").set_index("feature_id")
    assert (counts.dtypes == "int64").all()
    assert (counts.sum(axis=0) == led.loc[counts.columns, "assigned_to_features"]).all()
    truth = pd.read_csv(work / "data" / "truth_pairs.tsv", sep="\t").set_index("sample")
    assert ((counts.loc[V4] - truth.loc[counts.columns, "V4"]).abs() <= 0.05 * truth.loc[counts.columns, "V4"]).all()
    assert (led["decoy_host:testhost"] > 0).all() and (led["unmapped"] > 0).all()     # host and bacterial reads
    assert (counts[counts.index.str.endswith("|reference")].sum(axis=1) == 0).all()


def test_cached_counts_are_reused_only_for_the_same_panel(work):
    again = cli("panel-quant", "--panel", work / "panel", "--input-dir", work / "data" / "reads", "-o", work / "quant", "-t", 4)
    assert "RUN bowtie2" not in again
    panel2 = work / "panel2"
    cli("panel-build", "--run", work / "run", "-o", panel2, "--categories", "virus", "-t", 4)
    sheet = work / "sheet_of_the_first_count.tsv"          # inputs never lie inside the output folder
    sheet.write_text((work / "quant" / "samplesheet.used.tsv").read_text())
    out = cli("panel-quant", "--panel", panel2, "--samples", sheet,
              "-o", work / "quant", "-t", 4)
    assert "recounting" in out and "RUN bowtie2" in out
    # restore the counts of the first panel for the remaining tests
    cli("panel-quant", "--panel", work / "panel", "--input-dir", work / "data" / "reads", "-o", work / "quant", "-t", 4)


def rscript():
    try:
        r = da.find_rscript(None)
    except Exception:
        return None
    return r if da.has_ancombc(r) else None


def test_ancombc2_finds_the_one_differential_virus(work):
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package (create the environment-r.yml environment)")
    out = work / "da"
    log = cli("da", "--counts", work / "quant" / "counts.tsv", "--features", work / "quant" / "features.tsv",
              "--metadata", work / "data" / "metadata.tsv", "--group", "group", "--reference-level", "A",
              "-o", out, "--rscript", r)
    assert "fewer than 50" in log                        # small-panel warning is printed ...
    text = (out / "ancombc2_results.tsv").read_text()
    assert "# WARNING:" in text and "ANCOMBC " in text    # ... and written into the result header
    res = pd.read_csv(out / "ancombc2_results.tsv", sep="\t", comment="#").set_index("feature_id")
    assert bool(res.loc[V4, "diff_robust_groupB"]) and res.loc[V4, "lfc_groupB"] > 1.0
    others = res.drop(index=V4)
    assert not others["diff_robust_groupB"].astype(bool).any(), others["diff_robust_groupB"].to_dict()
    for f in ("volcano.png", "volcano.pdf", "feature_abundance.png", "feature_abundance.pdf", "session_info.txt"):
        assert (out / f).exists(), f


def test_da_refuses_a_group_with_one_sample(work, tmp_path):
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package")
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample\tgroup\nA1\tA\nB1\tB\nB2\tB\n")
    out = cli("da", "--counts", work / "quant" / "counts.tsv", "--features", work / "quant" / "features.tsv",
              "--metadata", meta, "--group", "group", "-o", tmp_path / "da", "--rscript", r, ok=False)
    assert "at least 2 samples" in out


# ------------------------------------------------------------ classified taxa as background features (taxa-table)

def run_da(work, counts_dir, out, r):
    log = cli("da", "--counts", counts_dir / "counts.tsv", "--features", counts_dir / "features.tsv",
              "--metadata", work / "data" / "metadata.tsv", "--group", "group", "--reference-level", "A",
              "-o", out, "--rscript", r)
    res = pd.read_csv(out / "ancombc2_results.tsv", sep="\t", comment="#").set_index("feature_id")
    return log, res


def test_merge_with_a_panel_counted_on_all_reads_requires_the_panel_lineage(work, tmp_path):
    out = cli("taxa-table", "--input-dir", work / "data" / "reports", "--suffix", ".kreport", "--report-format", "kraken2", "--rank", "S",
              "--merge-panel", work / "quant" / "counts.tsv", "-o", tmp_path / "taxa", ok=False)
    assert "--panel-lineage" in out and not (tmp_path / "taxa" / "merged").exists()


def test_background_taxa_lift_the_small_panel_limit_and_the_differential_virus_is_still_found(work):
    """Panel counted on ALL reads + the classifier's species (Viruses left out of the reports) in one table."""
    taxa_dir = work / "taxa_all"
    cli("taxa-table", "--input-dir", work / "data" / "reports", "--suffix", ".kreport", "--report-format", "kraken2", "--rank", "S",
        "--namespace", "k2_test", "--merge-panel", work / "quant" / "counts.tsv", "--panel-lineage", "Viruses",
        "-o", taxa_dir)
    led = pd.read_csv(taxa_dir / "ledger.tsv", sep="\t")
    assert (led[["unclassified", "dropped_lineages", "kept_at_rank", "unresolved_above_rank", "rounding_residual"]]
            .sum(axis=1) == led["input_fragments"]).all() and (led["dropped_lineages"] > 0).all()
    merged = pd.read_csv(taxa_dir / "merged" / "counts.tsv", sep="\t").set_index("feature_id")
    feats = pd.read_csv(taxa_dir / "merged" / "features.tsv", sep="\t")
    assert (feats["source"] == "background").sum() == 60 and not merged.index.str.contains(r"\|70001\|").any()
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package")
    log, res = run_da(work, taxa_dir / "merged", work / "da_all_bg", r)
    assert "fewer than 50" not in log and len(res) >= 60
    assert bool(res.loc[V4, "diff_robust_groupB"]) and res.loc[V4, "lfc_groupB"] > 1.0
    others = res.drop(index=V4)
    assert others["diff_robust_groupB"].astype(bool).sum() == 0, others.index[others["diff_robust_groupB"].astype(bool)]


def test_classified_taxa_and_unclassified_side_features_in_one_table(work):
    """Disjoint read sets: the classifier's taxa + a panel counted only on the reads the classifier left over."""
    quant_un = work / "quant_un"
    cli("panel-quant", "--panel", work / "panel", "--input-dir", work / "run" / "unclassified_reads",
        "-o", quant_un, "-t", 4)
    prov = __import__("json").loads((quant_un / "provenance.json").read_text())
    assert prov["reads_scope"] == "unclassified" and prov["input_summary"].endswith("input_summary.tsv")
    assert __import__("json").loads((work / "quant" / "provenance.json").read_text())["reads_scope"] == "all"
    taxa_dir = work / "taxa_un"
    cli("taxa-table", "--input-dir", work / "data" / "reports", "--suffix", ".kreport", "--report-format", "kraken2", "--rank", "S",
        "--merge-panel", quant_un / "counts.tsv", "-o", taxa_dir)                # no --panel-lineage needed
    rec = pd.read_csv(taxa_dir / "merged" / "reconciliation.tsv", sep="\t").set_index("sample")
    assert (rec["status"] == "ok").all() and (rec["difference"] == 0).all()
    truth = pd.read_csv(work / "data" / "reports" / "A1.kreport", sep="\t", header=None)
    assert rec.loc["A1", "total_input_reads"] == truth[1].iloc[:2].sum()            # unclassified + root rows
    parts = ["classified_kept_at_rank", "classified_dropped_lineages", "classified_unresolved_above_rank",
             "classified_rounding_residual", "panel_features", "panel_rest"]
    assert (rec[parts].sum(axis=1) == rec["total_input_reads"]).all()
    merged = pd.read_csv(taxa_dir / "merged" / "counts.tsv", sep="\t").set_index("feature_id")
    assert merged.index.str.contains(r"\|70001\|").any()       # the classifier's own virus species stays a feature
    all_reads = pd.read_csv(work / "quant" / "counts.tsv", sep="\t").set_index("feature_id")
    assert (merged.loc[V4] == all_reads.loc[V4]).all()          # this virus is entirely on the unclassified side
    r = rscript()
    if r is None:
        pytest.skip("no Rscript with the ANCOMBC package")
    log, res = run_da(work, taxa_dir / "merged", work / "da_un_bg", r)
    assert "fewer than 50" not in log
    assert bool(res.loc[V4, "diff_robust_groupB"]) and res.loc[V4, "lfc_groupB"] > 1.0
    others = res.drop(index=V4)
    assert others["diff_robust_groupB"].astype(bool).sum() == 0, others.index[others["diff_robust_groupB"].astype(bool)]
