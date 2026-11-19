"""End-to-end runs on the synthetic data set (needs the micromamba environment; no downloads)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).parent
TOOLS = ("spades.py", "megahit", "minimap2", "blastn", "diamond", "bowtie2", "samtools", "seqtk")
pytestmark = [pytest.mark.e2e,
              pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS),
                                 reason="external tools not on PATH (activate the environment)")]


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("syn")
    subprocess.run([sys.executable, str(HERE / "make_synthetic.py"), str(d)], check=True)
    return d


def run_cli(*args):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr


@pytest.mark.parametrize("mode", ["rna", "dna"])
def test_pipeline_explains_unclassified_reads(data, tmp_path, mode):
    out = tmp_path / f"out_{mode}"
    log = run_cli("run", "--input-dir", data / "reads", "--refs", data / "refs", "--host", f"testhost={data / 'host.fa'}",
                  "--mode", mode, "-o", out, "-t", 4, "--memory", 8, "--subsample", "all")
    t = pd.read_csv(out / "tables" / "read_category_by_sample.tsv", sep="\t")
    frac = t[[c for c in t.columns if c.startswith("frac_")]]
    assert (frac.sum(axis=1) - 1).abs().max() < 1e-9
    assert set(t["sample"]) == {"S1", "S2", "S3"}
    # virus is the largest class in every sample, and both viruses (genome hit and protein-only hit) are found
    for _, r in frac.iterrows():
        assert r.idxmax() == "frac_virus", r.to_dict()
    assert (t["frac_virus"] > 0.5).all()
    assert (t["frac_host:testhost"] > 0.05).all()
    refs = pd.read_csv(out / "tables" / "reference_by_sample.tsv", sep="\t")
    titles = " ".join(refs.loc[refs.category == "virus", "best_hit_title"].astype(str))
    assert "Test virus A" in titles and "Test virus B" in titles
    inputs = pd.read_csv(out / "tables" / "input_summary.tsv", sep="\t")
    s1 = inputs.set_index("sample").loc["S1"]
    assert s1["classifier"] == "kraken2" and s1["classified_reads"] > 0
    assert (s1["total_reads"] - s1["classified_reads"]) == s1["unclassified_reads"]
    for f in ("unclassified_reads/S1_R1.fastq.gz", "unclassified_reads/S2_R2.fastq.gz", "unclassified_reads/S3.fastq.gz",
              "figures/read_categories.png", "figures/read_categories.pdf", "figures/sankey_all_samples.html",
              "figures/sankey/S1.html", "report.html", "samplesheet.used.tsv", "tables/contig_category.tsv",
              "tables/resources.tsv"):
        assert (out / f).exists(), f
    # a rerun resumes from the done markers: no external tool runs again, and no marker is rewritten
    markers = {m.name: m.stat().st_mtime_ns for m in (out / "work" / ".done").iterdir()}
    again = run_cli("run", "--input-dir", data / "reads", "--refs", data / "refs",
                    "--host", f"testhost={data / 'host.fa'}", "--mode", mode, "-o", out, "-t", 4,
                    "--memory", 8, "--subsample", "all")
    assert "done" in again
    assert "RUN " not in again, [l for l in again.splitlines() if "RUN " in l][:5]
    assert {m.name: m.stat().st_mtime_ns for m in (out / "work" / ".done").iterdir()} == markers


def test_extract_only_writes_unclassified_fastq(data, tmp_path):
    out = tmp_path / "ex"
    run_cli("extract", "--input-dir", data / "reads", "-o", out)
    inputs = pd.read_csv(out / "tables" / "input_summary.tsv", sep="\t").set_index("sample")
    assert (out / "unclassified_reads" / "S1_R1.fastq.gz").exists()
    assert inputs.loc["S2", "classifier"] == "centrifuge"
    assert inputs.loc["S2", "unclassified_reads"] < inputs.loc["S2", "total_reads"]


def test_extract_thresholds_move_low_score_reads_to_the_unclassified_side(data, tmp_path):
    """Centrifuge(r) assignments below --min-score / --min-hitlen are extracted and counted, never lost."""
    sheet = tmp_path / "s2.tsv"
    sheet.write_text("sample\tr1\tr2\tclassifier_output\tclassifier\n"
                     f"S2\t{data / 'reads' / 'S2_R1.fastq.gz'}\t{data / 'reads' / 'S2_R2.fastq.gz'}\t"
                     f"{data / 'reads' / 'S2.centrifuge'}\tcentrifuge\n")
    run_cli("extract", "--samples", sheet, "-o", tmp_path / "plain")
    run_cli("extract", "--samples", sheet, "-o", tmp_path / "thr", "--min-score", 100, "--min-hitlen", 30)
    a = pd.read_csv(tmp_path / "plain" / "tables" / "input_summary.tsv", sep="\t").iloc[0]
    b = pd.read_csv(tmp_path / "thr" / "tables" / "input_summary.tsv", sep="\t").iloc[0]
    assert a["low_score_reads"] == 0 and b["low_score_reads"] > 0
    assert b["classified_reads"] == a["classified_reads"] - b["low_score_reads"]
    for r in (a, b):
        assert r["classified_reads"] + r["low_score_reads"] + r["no_hit_reads"] == r["total_reads"]
        assert r["unclassified_reads"] == r["low_score_reads"] + r["no_hit_reads"]
    assert b["unclassified_reads"] == a["unclassified_reads"] + b["low_score_reads"]
    # the flags do not exist for Kraken2 output (sample S1 in the folder)
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", "extract", "--input-dir", str(data / "reads"),
                        "-o", str(tmp_path / "k"), "--min-score", "100"], capture_output=True, text=True)
    assert p.returncode == 1 and "Kraken2" in p.stderr


def test_run_without_references_fails_with_hint(data, tmp_path):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", "run", "--input-dir", str(data / "reads"),
                        "--host", str(data / "host.fa"), "-o", str(tmp_path / "x")], capture_output=True, text=True)
    assert p.returncode == 1 and "download-refs" in p.stderr


def test_changing_subsample_recounts_instead_of_reusing_old_counts(data, tmp_path):
    """Regression: read-back counts made with one --subsample were reused after switching to another."""
    out = tmp_path / "sub"
    base = ("run", "--input-dir", data / "reads", "--refs", data / "refs", "--host", f"testhost={data / 'host.fa'}",
            "--mode", "dna", "-o", out, "-t", 4, "--memory", 8)
    run_cli(*base, "--subsample", 800)
    first = pd.read_csv(out / "tables" / "contig_counts.tsv", sep="\t").set_index("contig").sum()
    assert (pd.read_csv(out / "tables" / "input_summary.tsv", sep="\t")["reads_used"] == 800).all()
    log = run_cli(*base, "--subsample", "all")
    inputs = pd.read_csv(out / "tables" / "input_summary.tsv", sep="\t").set_index("sample")
    assert (inputs["reads_used"] == inputs["unclassified_reads"]).all()
    assert "recomputing" in log
    second = pd.read_csv(out / "tables" / "contig_counts.tsv", sep="\t").set_index("contig").sum()
    t = pd.read_csv(out / "tables" / "read_category_by_sample.tsv", sep="\t").set_index("sample")
    assert (second > first[second.index]).all()
    assert (second + t.loc[second.index, "n_not_assembled"] == inputs.loc[second.index, "reads_used"]).all()
