"""End-to-end test of the robustness extension on the synthetic two-group data (docs/DESIGN_robustness.md):

run -> strain-evidence (with reads) -> panel-build --strain-clusters (98 and 99 %) -> panel-quant --ambiguous em
-> em-uncertainty -> da-robustness (ANCOM-BC2, ALDEx2, LinDA) -> da-resolutions.

Virus V4 is six times more abundant in group B, nothing else differs. Needs the external tools and the R
environment with ANCOMBC, ALDEx2 and MicrobiomeStat (UNCLASSIFIED_META_RSCRIPT); skipped otherwise. Slow (minutes).
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import manifest

HERE = Path(__file__).parent
TOOLS = ("spades.py", "minimap2", "blastn", "diamond", "bowtie2", "samtools", "seqtk", "seqkit")
V4 = "NC_900004.1"


def _rscript_with_all_methods():
    r = os.environ.get("UNCLASSIFIED_META_RSCRIPT") or shutil.which("Rscript")
    if not r:
        return None
    p = subprocess.run([r, "-e", "for (p in c('ANCOMBC','ALDEx2','MicrobiomeStat')) "
                                 "if (!requireNamespace(p, quietly=TRUE)) quit(status=1)"], capture_output=True)
    return r if p.returncode == 0 else None


RSCRIPT = _rscript_with_all_methods()
pytestmark = [pytest.mark.e2e,
              pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS),
                                 reason="external tools not on PATH (activate the environment)"),
              pytest.mark.skipif(RSCRIPT is None, reason="R environment with ANCOMBC, ALDEx2, MicrobiomeStat missing")]


def cli(*args):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr


def tsv(path):
    return pd.read_csv(path, sep="\t", comment=None, skiprows=sum(1 for _ in open(path) if _.startswith("#")),
                       keep_default_na=True)


@pytest.fixture(scope="module")
def work(tmp_path_factory):
    d = tmp_path_factory.mktemp("robust")
    data = d / "data"
    subprocess.run([sys.executable, str(HERE / "make_synthetic_da.py"), str(data)], check=True)
    cli("run", "--input-dir", data / "reads", "--refs", data / "refs", "--host", f"testhost={data / 'host.fa'}",
        "--mode", "rna", "-t", 4, "--memory", 8, "--subsample", "all", "-o", d / "run")
    cli("strain-evidence", "--run", d / "run", "--input-dir", data / "reads", "--max-fragments", 4000, "-t", 4,
        "-o", d / "strains")
    meta = d / "meta.tsv"
    rows = [f"{g}{k}\t{g}\t{'b1' if k <= 3 else 'b2'}" for g in "AB" for k in range(1, 7)]
    meta.write_text("sample\tgroup\tbatch\n" + "\n".join(rows) + "\n")
    for res in (98, 99):
        cli("panel-build", "--run", d / "run", "-o", d / f"panel{res}", "--categories", "virus,unclassified",
            "--strain-clusters", d / "strains" / "viral_strain_clusters.tsv", "--strain-resolution", res, "-t", 4)
        cli("panel-quant", "--panel", d / f"panel{res}", "--input-dir", data / "reads", "--ambiguous", "em",
            "-o", d / f"quant{res}", "-t", 4)
        cli("em-uncertainty", "--quant", d / f"quant{res}", "--bootstraps", 10, "--realizations", 5,
            "-o", d / f"unc{res}", "-t", 4)
        cli("da-robustness", "--em-uncertainty", d / f"unc{res}", "--features", d / f"quant{res}" / "features.tsv",
            "--metadata", meta, "--group", "group", "--unit", "batch", "--bootstrap-fits", 2,
            "--realization-fits", 2, "--prv-cuts", "0.5", "--mc-samples", 16, "--rscript", RSCRIPT,
            "-t", 8, "-o", d / f"da{res}")
    cli("da-resolutions", "--run", f"98={d / 'da98'}", "--run", f"99={d / 'da99'}", "-o", d / "resolutions")
    return d


def test_strain_clusters_are_labelled_separately_from_haplotypes(work):
    cl = tsv(work / "strains" / "viral_strain_clusters.tsv")
    allowed = {"operational_sequence_cluster", "candidate_strain_lineage", "phased_viral_haplotype",
               "reference_supported_viral_strain"}
    assert set(cl["class"]) <= allowed and {97, 98, 99} <= set(cl["resolution"].astype(int))
    phased = cl[cl["class"] == "phased_viral_haplotype"]
    assert phased["haplotype_blocks"].astype(str).str.len().gt(0).all()      # every haplotype cites its blocks
    for name in ("viral_cluster_overlap.tsv", "viral_haplotype_evidence.tsv", "phasing_limits.tsv"):
        assert (work / "strains" / name).exists()
    manifest.verify(work / "strains")


def test_observed_and_em_counts_are_kept_apart(work):
    unc = work / "unc98"
    manifest.verify(unc)
    uniq = tsv(unc / "unique_counts.tsv").set_index("feature_id")
    quant_uniq = tsv(work / "quant98" / "unique_counts.tsv").set_index("feature_id")
    em = tsv(unc / "em_expected_counts.tsv").set_index("feature_id")
    pd.testing.assert_frame_equal(uniq, quant_uniq)
    assert (em.loc[uniq.index] + 1e-6 >= uniq).all().all()                  # EM only adds shared fragments
    assert len(list((unc / "bootstrap").glob("counts_*.tsv"))) == 10
    sup = tsv(unc / "viral_unique_read_support.tsv")
    assert set(sup["assignment_confidence"]) <= {"high", "moderate", "low"}


def test_v4_is_found_by_all_three_frameworks_and_reported_without_validation_wording(work):
    out = work / "da98"
    manifest.verify(out)
    conc = tsv(out / "cross_method_concordance.tsv")
    v4 = conc[conc["feature_id"].str.contains(V4, regex=False)]
    assert not v4.empty and (v4["n_methods_significant"] == 3).any()
    assert (v4["concordance"] == "all_significant_same_direction").any()
    others = conc[~conc["feature_id"].str.contains(V4, regex=False)]
    assert (others["n_methods_significant"] < 3).all()                       # nothing else in all three
    stab = tsv(out / "differential_abundance_stability.tsv")
    for c in ("boot_direction_stability", "boot_significance_stability", "real_direction_stability"):
        assert c in stab.columns
    summ = tsv(out / "robustness_summary.tsv")
    for c in ("assignment_confidence", "abundance_stability", "da_support", "cross_method_agreement",
              "robustness_basis", "conclusion"):
        assert c in summ.columns
    rep = (out / "robustness_report.md").read_text()
    assert "complementary differential abundance frameworks" in rep and "validated by" not in rep.lower()
    assert (out / "replicate_sensitivity.tsv").exists()


def test_resolutions_are_compared_per_parent_feature(work):
    mr = tsv(work / "resolutions" / "multi_resolution_robustness.tsv")
    assert {"parent", "term", "resolution_stability"} <= set(mr.columns)
    assert mr["parent"].astype(str).str.contains(V4, regex=False).any()
