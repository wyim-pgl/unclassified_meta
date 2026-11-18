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
    # a rerun resumes from the done markers
    assert "done" in run_cli("run", "--input-dir", data / "reads", "--refs", data / "refs",
                             "--host", f"testhost={data / 'host.fa'}", "--mode", mode, "-o", out, "-t", 4,
                             "--memory", 8, "--subsample", "all")


def test_extract_only_writes_unclassified_fastq(data, tmp_path):
    out = tmp_path / "ex"
    run_cli("extract", "--input-dir", data / "reads", "-o", out)
    inputs = pd.read_csv(out / "tables" / "input_summary.tsv", sep="\t").set_index("sample")
    assert (out / "unclassified_reads" / "S1_R1.fastq.gz").exists()
    assert inputs.loc["S2", "classifier"] == "centrifuge"
    assert inputs.loc["S2", "unclassified_reads"] < inputs.loc["S2", "total_reads"]


def test_run_without_references_fails_with_hint(data, tmp_path):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", "run", "--input-dir", str(data / "reads"),
                        "--host", str(data / "host.fa"), "-o", str(tmp_path / "x")], capture_output=True, text=True)
    assert p.returncode == 1 and "download-refs" in p.stderr
