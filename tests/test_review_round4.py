"""Regression tests for the fourth code review (2026-11-18): preflight of every write, concurrency, cached-step
verification, exact manifest coverage, producer markers, reference hashes, download caches, seal-refs, Sankey."""
import gzip
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import control, io_utils, panel, plots, refs, runner, taxa
from unclassified_meta.config import Sample

from test_ported_units import K2, fastq, panel_dir, write
from test_review_fixes import HAVE_PANEL_TOOLS, counts_with_a_structural_zero, da_cli, rscript, tiny_run
from test_review_round2 import fastq_pair


def cli(*args):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


# ------------------------------------------------------------------ 1 CRITICAL: preflight before the first write

def test_cli_written_files_never_touch_an_input(tmp_path):
    (tmp_path / "reads").mkdir()
    r1, r2 = fastq_pair(tmp_path / "reads", "s", "ab")
    out = tmp_path / "out"
    out.mkdir()
    os.symlink(r1, out / "samplesheet.used.tsv")                     # the CLI's first file links to an input
    before = r1.read_bytes()
    rc, log = cli("extract", "--input-dir", tmp_path / "reads", "-o", out)
    assert rc != 0 and r1.read_bytes() == before


def test_an_input_inside_an_owned_output_folder_is_refused(tmp_path):
    q = tmp_path / "quant"
    (q / "work").mkdir(parents=True)
    fastq(q / "work" / "s.fastq.gz", "ab")                           # input inside panel-quant's scratch folder
    with pytest.raises(ValueError, match="inside"):
        io_utils.preflight("panel-quant", [q / "work" / "s.fastq.gz"], owned_dirs=[q / "work"])


# ------------------------------------------------------------------ 2 HIGH: unique temporaries and folder locks

def test_staged_writers_never_share_or_delete_each_others_temporary_files(tmp_path):
    dst = tmp_path / "x.tsv"
    a, b = io_utils.staged(dst), io_utils.staged(dst)
    ta, tb = a.__enter__(), b.__enter__()
    assert ta != tb
    ta.write_text("A")
    tb.write_text("B")
    b.__exit__(None, None, None)
    a.__exit__(None, None, None)
    assert dst.read_text() == "A" and not ta.exists() and not tb.exists()


def test_a_second_command_on_the_same_folder_is_refused(tmp_path):
    with io_utils.folder_lock(tmp_path, "run"):
        with pytest.raises(RuntimeError, match="another run is using this folder"):
            with io_utils.folder_lock(tmp_path, "run"):
                pass


# ------------------------------------------------------------------ 3 HIGH: cached steps verify their outputs

def test_a_cached_step_whose_output_changed_is_recomputed(tmp_path):
    steps = runner.Steps(tmp_path)
    out = write(tmp_path / "contigs.fa", ">c1\nACGT\n")
    steps.mark("assemble", "k1", outputs=[out])
    assert steps.done("assemble", "k1")
    out.write_text(">c1\nTTTT\n")                                     # same identifier, new sequence
    assert not steps.done("assemble", "k1")


def test_classification_key_follows_the_contig_content(tmp_path):
    from unclassified_meta import pipeline
    fa = write(tmp_path / "contigs.fa", ">c1\nACGT\n")
    cfg = {"references": {}, "thresholds": {}, "priority": []}
    k1 = pipeline.classification_key(fa, cfg)
    fa.write_text(">c1\nTTTT\n")
    assert pipeline.classification_key(fa, cfg) != k1


# ------------------------------------------------------------------ 4 + 5 HIGH: exactly the sealed files are used

def test_merge_refuses_a_count_file_not_covered_by_the_manifest(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", in_frag=(20, 20), summary=str(summ))
    alt = pc.parent / "edited_counts.tsv"
    alt.write_text(pc.read_text())                                   # valid seal on counts.tsv, but another file used
    write(tmp_path / "s1.kreport", K2)
    write(tmp_path / "s2.kreport", K2)
    out = tmp_path / "taxa"
    taxa.build_table(taxa.find_reports(tmp_path, ".kreport"), out, "S", report_format="kraken2")
    with pytest.raises(ValueError, match="not covered|manifest"):
        taxa.merge_with_panel(out, alt, None, out / "merged")


def test_generated_folder_without_a_complete_seal_is_refused_by_da(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    (tmp_path / ".unclassified_meta_generation").write_text("taxa-table --merge-panel\n")   # producer started, no seal
    rc, log = da_cli("--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
                     tmp_path / "metadata.tsv", "--group", "group", "-o", tmp_path / "d",
                     "--rscript", rscript() or "Rscript")
    assert rc != 0 and "manifest" in log
    (tmp_path / ".unclassified_meta_generation").unlink()
    rc, log = da_cli("--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
                     tmp_path / "metadata.tsv", "--group", "group", "-o", tmp_path / "d2",
                     "--rscript", rscript() or "Rscript")
    assert rc != 0 and "--unsealed-input" in log                     # hand-made tables need the explicit flag


def test_a_manifest_without_artifacts_is_not_a_seal(tmp_path):
    from unclassified_meta import manifest
    with pytest.raises(ValueError):
        manifest.seal(tmp_path, [], {"step": "x"})                   # the producer cannot seal nothing
    (tmp_path / manifest.NAME).write_text(json.dumps({"step": "x", "artifacts": {}}))
    with pytest.raises(ValueError):
        manifest.verify(tmp_path)                                    # and a consumer does not accept it


# ------------------------------------------------------------------ 6 HIGH: references verified against the run

@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_panel_build_refuses_a_reference_changed_since_the_run(tmp_path):
    run, _ = tiny_run(tmp_path)
    vir = tmp_path / "refseq_viral.fa"
    vir.write_text(vir.read_text().replace("Virus one", "Virus one changed"))       # same accessions, new content
    with pytest.raises(ValueError, match="recorded"):
        panel.build_panel(run, tmp_path / "p", ["virus"], [], None, None, None, threads=1)


# ------------------------------------------------------------------ 7 HIGH: cached downloads belong to one release

def test_cached_download_from_another_release_is_fetched_again(tmp_path, monkeypatch):
    payload = gzip.compress(b">v\nACGT\n")
    calls = []

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=0):
        calls.append(req.full_url)
        return Resp(payload)
    monkeypatch.setattr(refs.urllib.request, "urlopen", fake_urlopen)
    dst = tmp_path / "src" / "viral.1.1.genomic.fna.gz"
    refs.download("https://example.org/viral.1.1.genomic.fna.gz", dst, release="237")
    refs.download("https://example.org/viral.1.1.genomic.fna.gz", dst, release="237")
    assert len(calls) == 1                                          # same release: cached part reused
    refs.download("https://example.org/viral.1.1.genomic.fna.gz", dst, release="238")
    assert len(calls) == 2                                          # other release: fetched again


# ------------------------------------------------------------------ 8 MEDIUM: failed sampling publishes nothing

def test_failed_control_sampling_leaves_the_previous_output(tmp_path, monkeypatch):
    fastq(tmp_path / "s.fastq.gz", "abc")
    dst = tmp_path / "ctrl.fastq.gz"
    dst.write_text("previous")

    class FakeProc:
        def __init__(self, *a, **k):
            self.stdout = io.StringIO("@a\nACGT\n+\nIIII\n")

        def wait(self):
            return 1
    monkeypatch.setattr(control.shutil, "which", lambda x: "/bin/seqtk")
    monkeypatch.setattr(control.subprocess, "Popen", FakeProc)
    with pytest.raises(Exception):
        control.reverse_fastq(tmp_path / "s.fastq.gz", dst, 1, 11)
    assert dst.read_text() == "previous"


# ------------------------------------------------------------------ 9 MEDIUM: adopting an existing reference folder

def test_seal_refs_adopts_a_manual_reference_folder(tmp_path):
    from unclassified_meta import manifest
    d = tmp_path / "refs"
    for sub, files in (("viral", {"refseq_viral.fna": ">v\nACGT\n", "refseq_viral_protein.faa": ">p\nMKV\n"}),
                       ("rrna", {"plant_rrna.fasta": ">r\nACGT\n"}), ("protein", {"sprot_viridiplantae.faa": ">p\nMKV\n"})):
        (d / sub).mkdir(parents=True)
        for n, t in files.items():
            write(d / sub / n, t)
    rc, log = cli("seal-refs", d)
    assert rc == 0, log
    for sub in ("viral", "rrna", "protein"):
        manifest.verify(d / sub)
    assert "unknown" in (d / "REFS_VERSION.tsv").read_text()


# ------------------------------------------------------------------ 10 LOW: Sankey flows are conserved

def test_reference_flows_sum_exactly_to_their_parent():
    assert sum(plots.allocate(5, [1, 1, 1])) == 5
    assert plots.allocate(10, [3, 1]) == [8, 2] or sum(plots.allocate(10, [3, 1])) == 10
