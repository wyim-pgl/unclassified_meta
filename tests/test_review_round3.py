"""Regression tests for the third code review (2026-11-18): output safety across samples and commands,
generation manifests, extraction validation, structural detection, decoys, provenance, DA labels, Sankey."""
import gzip
import json
import os
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import control, panel, plots, refs, summarize
from unclassified_meta.config import Sample, detect_samples
from unclassified_meta.pipeline import extract_all

from test_ported_units import CFR_HEADER, cfr, fastq, write
from test_review_fixes import HAVE_PANEL_TOOLS, counts_with_a_structural_zero, da_cli, rscript, tiny_run
from test_review_round2 import fastq_pair

ROOT = Path(__file__).parent.parent


# ------------------------------------------------------------------ 1 + 2 CRITICAL: no output may be any input

def test_one_samples_output_never_replaces_another_samples_input(tmp_path):
    fastq(tmp_path / "a.fastq.gz", "ab")
    out = tmp_path / "out"
    (out / "unclassified_reads").mkdir(parents=True)
    b_in = out / "unclassified_reads" / "A.fastq.gz"                  # B's input sits where A's output goes
    fastq(b_in, "xyz")
    before = b_in.read_bytes()
    with pytest.raises(ValueError, match="input"):
        extract_all([Sample("A", tmp_path / "a.fastq.gz"), Sample("B", b_in)], out)
    assert b_in.read_bytes() == before


def test_control_reads_never_write_over_their_input(tmp_path):
    out = tmp_path / "ctl"
    (out / "control_reads").mkdir(parents=True)
    src = out / "control_reads" / "s_ctrl.fastq.gz"                    # same path as the control output
    fastq(src, "abc")
    before = src.read_bytes()
    with pytest.raises(ValueError, match="input"):
        control.make_controls([Sample("s", src)], out, n_reads=None)
    assert src.read_bytes() == before
    other = tmp_path / "orig.fastq.gz"
    fastq(other, "abc")
    os.link(other, out / "control_reads" / "t_ctrl.fastq.gz")        # a hard link to the input
    keep = other.read_bytes()
    with pytest.raises(ValueError, match="input"):
        control.make_controls([Sample("t", other)], out, n_reads=None)
    assert other.read_bytes() == keep


# ------------------------------------------------------------------ 3: every FASTQ read must be a classifier read

def test_fastq_reads_must_be_exactly_the_classifier_reads(tmp_path):
    fastq(tmp_path / "s.fastq.gz", ["x", "x", "c"])
    k = write(tmp_path / "s.kraken", "C\ta\t562\t1\t0:1\nC\tb\t562\t1\t0:1\nU\tc\t0\t1\t0:1\n")
    with pytest.raises(ValueError, match="classifier"):
        extract_all([Sample("s", tmp_path / "s.fastq.gz", None, k, "kraken2")], tmp_path / "o")


def test_mates_are_checked_also_without_a_classifier(tmp_path):
    r1, r2 = fastq_pair(tmp_path, "s", ["a", "b"], ["b", "a"])
    with pytest.raises(ValueError, match="mate"):
        extract_all([Sample("s", r1, r2, None, "none")], tmp_path / "o")


# ------------------------------------------------------------------ 4: classifier files recognised by structure

def test_classifier_files_are_recognised_by_their_content(tmp_path):
    fastq(tmp_path / "reporter1.fastq.gz", "ab")
    write(tmp_path / "reporter1.cfr", CFR_HEADER + cfr("a", 562, 400, 60) + cfr("b", 0, 0, 0))
    fastq(tmp_path / "s2.fastq.gz", ["1", "2"])
    write(tmp_path / "s2.cfr", "1\tNC_1\t562\t400\t0\t60\t150\t1\n2\tunclassified\t0\t0\t0\t0\t150\t1\n")  # numeric IDs
    fastq(tmp_path / "s3.fastq.gz", "ab")
    with gzip.open(tmp_path / "s3.cfr.tsv.gz", "wt") as fh:
        fh.write(CFR_HEADER + cfr("a", 0, 0, 0) + cfr("b", 0, 0, 0))
    got = {s.name: s.classifier for s in detect_samples(tmp_path)}
    assert got == {"reporter1": "centrifuger", "s2": "centrifuger", "s3": "centrifuger"}


def test_plain_and_gzipped_copies_of_one_file_are_ambiguous(tmp_path):
    fastq(tmp_path / "s.fastq.gz", "ab")
    (tmp_path / "s.fastq").write_text("@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="s.fastq"):
        detect_samples(tmp_path)


# ------------------------------------------------------------------ 5: a run is used only as a sealed generation

@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_panel_build_refuses_a_run_whose_contigs_changed_after_it_was_sealed(tmp_path):
    from unclassified_meta import manifest
    run, _ = tiny_run(tmp_path)
    manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "run"})
    panel.build_panel(run, tmp_path / "p1", ["virus"], [], None, None, None, threads=1)
    fa = run / "contigs" / "contigs.fa"
    fa.write_text(fa.read_text().replace(">c1", ">c1_new"))           # an interrupted rerun replaced the contigs
    with pytest.raises(ValueError, match="changed|manifest"):
        panel.build_panel(run, tmp_path / "p2", ["virus"], [], None, None, None, threads=1)


# ------------------------------------------------------------------ 6: da refuses unpublished quant tables

def test_da_refuses_quant_tables_without_a_valid_seal(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    write(tmp_path / "ledger.tsv", "sample\tinput_fragments\n")      # a panel-quant folder, but never sealed
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
                     tmp_path / "metadata.tsv", "--group", "group", "-o", tmp_path / "d",
                     "--rscript", rscript() or "Rscript")
    assert rc != 0 and "manifest" in log and not (tmp_path / "d" / "ancombc2_results.tsv").exists()


# ------------------------------------------------------------------ 7: interrupted reference builds are redone

def test_interrupted_reference_artifacts_are_rebuilt(tmp_path, monkeypatch):
    src = write(tmp_path / "rrna_src.fa", ">r1\nACGT\n>r2\nGGCC\n")
    rd = tmp_path / "refs" / "rrna"
    rd.mkdir(parents=True)
    write(rd / "plant_rrna.fasta", ">r1\nACGT\n")                      # truncated copy from an interrupted run
    pd_ = tmp_path / "refs" / "protein"
    pd_.mkdir(parents=True)
    write(pd_ / "sprot_viridiplantae.dmnd", "7 bytes")                 # partial database at the final name
    calls = []

    def fake_sh(cmd, *a, **k):
        calls.append(cmd[:2])
        if cmd[0] == "makeblastdb":
            out = cmd[cmd.index("-out") + 1]
            for e in (".nhr", ".nin", ".nsq"):
                Path(out + e).write_text("x")
        if cmd[:2] == ["diamond", "makedb"]:
            Path(cmd[cmd.index("-d") + 1] + ".dmnd").write_text("db")
    monkeypatch.setattr(refs, "sh", fake_sh)
    monkeypatch.setattr("unclassified_meta.classify.sh", fake_sh)
    monkeypatch.setattr(refs, "require", lambda *a: None)
    monkeypatch.setattr(refs, "size_check", lambda *a, **k: None)
    monkeypatch.setattr(refs, "plan", lambda *a, **k: [])
    def fake_download(url, dst):
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        return write(Path(dst), ">p\nMKV\n")
    monkeypatch.setattr(refs, "download", fake_download)
    monkeypatch.setattr(refs, "concat_gz", lambda parts, out, transform=None: write(Path(out), ">p\nMKV\n") and 1)
    refs.download_refs(tmp_path / "refs", rrna_fasta=src, only=("rrna", "protein"))
    assert (rd / "plant_rrna.fasta").read_text() == src.read_text()
    assert ["diamond", "makedb"] in calls


# ------------------------------------------------------------------ 8: an empty decoy reference is an error

@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_empty_decoy_reference_is_refused(tmp_path):
    from unclassified_meta import manifest
    run, _ = tiny_run(tmp_path)
    manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "run"})
    empty = write(tmp_path / "empty_host.fa", "")
    with pytest.raises(ValueError, match="no sequence"):
        panel.build_panel(run, tmp_path / "p", ["virus"], [{"name": "h2", "fasta": str(empty)}], None, None, None,
                          threads=1)


# ------------------------------------------------------------------ 9: versions per reference, unknown stays unknown

def test_namespaces_use_the_version_of_each_reference():
    text = ("reference\tversion\tsource\tfile\tn_sequences\tbytes\tdate\n"
            "refseq_viral_genomic\tRefSeq release 100\tx\tf\t1\t1\t2025-01-01\n"
            "refseq_viral_protein\tRefSeq release 200\tx\tf\t1\t1\t2025-06-01\n")
    ns, _ = panel.refs_namespaces(None, None, text)
    assert ns["virus"] == "refseq_viral_r100" and ns["virus_protein"] == "refseq_viral_protein_r200"


def test_reused_database_without_a_record_is_recorded_as_unknown(tmp_path, monkeypatch):
    from unclassified_meta import manifest
    pd_ = tmp_path / "refs" / "protein"
    pd_.mkdir(parents=True)
    write(pd_ / "sprot_viridiplantae.dmnd", "db")
    manifest.seal(pd_, ["sprot_viridiplantae.dmnd"], {"step": "reference"})   # complete, but no version record
    monkeypatch.setattr(refs, "require", lambda *a: None)
    monkeypatch.setattr(refs, "size_check", lambda *a, **k: None)
    monkeypatch.setattr(refs, "plan", lambda *a, **k: [])
    vt = refs.download_refs(tmp_path / "refs", only=("protein",))
    row = [l for l in vt.read_text().splitlines() if l.startswith("sprot_viridiplantae")][0].split("\t")
    assert row[1] == "unknown" and row[-1] == "unknown"


# ------------------------------------------------------------------ 10: interaction terms are labelled as such

@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_interaction_coefficients_are_not_labelled_as_group_contrasts(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    meta = pd.read_csv(tmp_path / "metadata.tsv", sep="\t", dtype=str)
    meta["time"] = ["1", "2", "3"] * 4
    meta.to_csv(tmp_path / "metadata.tsv", sep="\t", index=False)
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
                     tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A", "--formula",
                     "group * time", "-o", tmp_path / "d", "--rscript", rscript())
    assert rc == 0, log
    header = [l for l in (tmp_path / "d" / "ancombc2_results.tsv").read_text().splitlines() if l.startswith("#")]
    labels = [l for l in header if l.startswith("# coefficients:")]
    assert labels and "groupB = B vs A" in labels[0] and "groupB:time = interaction" in labels[0]


# ------------------------------------------------------------------ 11: Sankey weights from unrounded shares

def test_sankey_weights_do_not_lose_small_but_heavily_scaled_shares():
    contigs = pd.DataFrame([{"contig": "c1", "category": "virus", "best_hit": "X", "best_hit_title": "X",
                             "method": "blastn", "identity": 90.0},
                            {"contig": "c2", "category": "virus", "best_hit": "Y", "best_hit_title": "Y",
                             "method": "blastn", "identity": 90.0}])
    t = summarize.reference_table(contigs, {"S1": Counter({"c1": 1}), "S2": Counter({"c2": 1})},
                                  {"S1": 1_000_000, "S2": 10})
    inputs = pd.DataFrame({"sample": ["S1", "S2"], "unclassified_reads": [4_000_000_000, 10], "reads_used": [1_000_000, 10]})
    w = dict(zip(t["best_hit"], plots.reference_weights(t, inputs)))
    assert w["X"] == pytest.approx(4000) and w["Y"] == pytest.approx(1)


# ------------------------------------------------------------------ 12: README commands match the CLI

def test_readme_report_commands_name_their_format():
    text = (ROOT / "README.md").read_text()
    cmds = re.findall(r"`(unclassified-meta taxa-table[^`]*)`", text)
    cmds += [l.strip() for l in text.splitlines() if l.strip().startswith("unclassified-meta taxa-table")]
    report_cmds = [c for c in cmds if "--per-read" not in c]
    assert report_cmds and all("--report-format" in c for c in report_cmds), \
        [c for c in report_cmds if "--report-format" not in c]
