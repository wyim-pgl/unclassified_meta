"""Regression tests for the second code review (2026-11-18): input safety, read accounting, caches, references,
ANCOM-BC2 design checks and report formats."""
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

from unclassified_meta import classify, config, extract, panel, plots, refs, summarize, taxa
from unclassified_meta.config import Sample, detect_samples, validate_samples
from unclassified_meta.pipeline import extract_all
from unclassified_meta.runner import StepError

from test_ported_units import CFR_HEADER, K2, cfr, fastq, panel_dir, write
from test_review_fixes import HAVE_PANEL_TOOLS, counts_with_a_structural_zero, da_cli, rscript, tiny_run


def fq_ids(path):
    with gzip.open(path, "rt") as fh:
        return [l[1:].split()[0] for i, l in enumerate(fh) if i % 4 == 0]


def fastq_pair(d: Path, name: str, ids1, ids2=None):
    fastq(d / f"{name}_R1.fastq.gz", ids1)
    fastq(d / f"{name}_R2.fastq.gz", ids2 if ids2 is not None else ids1)
    return d / f"{name}_R1.fastq.gz", d / f"{name}_R2.fastq.gz"


# ------------------------------------------------------------------ 1 CRITICAL: never write through to an input

def test_rerun_with_a_classifier_after_classifier_none_leaves_the_original_fastq_untouched(tmp_path):
    r1, r2 = fastq_pair(tmp_path, "s", "abcd")
    before = (r1.read_bytes(), r2.read_bytes())
    out = tmp_path / "out"
    extract_all([Sample("s", r1, r2, None, "none")], out)
    assert not (out / "unclassified_reads" / "s_R1.fastq.gz").is_symlink()
    c = write(tmp_path / "s.cfr", CFR_HEADER + cfr("a", 562, 400, 60) + cfr("b", 0, 0, 0) + cfr("c", 0, 0, 0)
              + cfr("d", 562, 400, 60))
    extract_all([Sample("s", r1, r2, c, "centrifuger")], out)
    assert (r1.read_bytes(), r2.read_bytes()) == before
    assert fq_ids(out / "unclassified_reads" / "s_R1.fastq.gz") == ["b", "c"]


def test_an_output_that_resolves_to_an_input_is_refused(tmp_path):
    r1, _ = fastq_pair(tmp_path, "s", "ab")
    out = tmp_path / "out"
    (out / "unclassified_reads").mkdir(parents=True)
    os.symlink(r1, out / "unclassified_reads" / "s.fastq.gz")      # a stale link pointing at the input
    before = r1.read_bytes()
    with pytest.raises(ValueError, match="input"):                   # an output that resolves to an input
        extract_all([Sample("s", r1, None, None, "none")], out)
    assert r1.read_bytes() == before


# ------------------------------------------------------------------ 2: exact ID coverage, no duplicates, mates in sync

def test_selected_reads_must_all_be_found_once(tmp_path):
    fastq(tmp_path / "s.fastq.gz", ["a", "a"])
    k = write(tmp_path / "s.kraken", "U\ta\t0\t1\t0:1\nU\tb\t0\t1\t0:1\n")
    with pytest.raises(ValueError, match="twice|duplicate|not found"):
        extract_all([Sample("s", tmp_path / "s.fastq.gz", None, k, "kraken2")], tmp_path / "o")


def test_mates_out_of_order_are_refused(tmp_path):
    r1, r2 = fastq_pair(tmp_path, "s", ["a", "b"], ["b", "a"])
    k = write(tmp_path / "s.kraken", "U\ta\t0\t1\t0:1\nU\tb\t0\t1\t0:1\n")
    with pytest.raises(ValueError, match="mate"):
        extract_all([Sample("s", r1, r2, k, "kraken2")], tmp_path / "o")


# ------------------------------------------------------------------ 3: generated names never collide

def test_generated_read_file_names_and_columns_must_be_unique(tmp_path):
    r1, r2 = fastq_pair(tmp_path, "A", "ab")
    fastq(tmp_path / "x.fastq.gz", "ab")
    with pytest.raises(ValueError, match="same output file"):
        validate_samples([Sample("A", r1, r2), Sample("A_R1", tmp_path / "x.fastq.gz")])
    with pytest.raises(ValueError, match="column"):
        validate_samples([Sample("of_reads_used_pooled", tmp_path / "x.fastq.gz")])


# ------------------------------------------------------------------ 4: detection of real naming schemes

def test_detection_of_star_unmapped_and_clean_read_names_with_exact_classifier_match(tmp_path):
    fastq(tmp_path / "19_RNA.Unmapped.out.mate1.gz", "ab")
    fastq(tmp_path / "19_RNA.Unmapped.out.mate2.gz", "ab")
    write(tmp_path / "19_RNA.cfr.tsv", CFR_HEADER + cfr("a", 0, 0, 0) + cfr("b", 0, 0, 0))
    fastq(tmp_path / "1_DNA_S21_clean_R1.fastq.gz", "ab")
    fastq(tmp_path / "1_DNA_S21_clean_R2.fastq.gz", "ab")
    write(tmp_path / "1_DNA_S21.kraken2", "U\ta\t0\t1\t0:1\nU\tb\t0\t1\t0:1\n")
    write(tmp_path / "1_DNA_S21.quant.kreport", "100.00\t2\t2\tR\t1\troot\n")       # a report is ignored
    got = {s.name: s for s in detect_samples(tmp_path)}
    assert set(got) == {"19_RNA", "1_DNA_S21"}
    assert got["19_RNA"].r2.name == "19_RNA.Unmapped.out.mate2.gz" and got["19_RNA"].classifier == "centrifuger"
    assert got["1_DNA_S21"].classifier == "kraken2"


def test_detection_refuses_orphan_and_ambiguous_classifier_outputs(tmp_path):
    fastq(tmp_path / "s1_R1.fastq.gz", "ab")
    fastq(tmp_path / "s1_R2.fastq.gz", "ab")
    write(tmp_path / "s9.cfr.tsv", CFR_HEADER + cfr("a", 0, 0, 0))
    with pytest.raises(ValueError, match="no reads"):
        detect_samples(tmp_path)
    (tmp_path / "s9.cfr.tsv").unlink()
    write(tmp_path / "s1.cfr", CFR_HEADER + cfr("a", 0, 0, 0))
    write(tmp_path / "s1.cfr.tsv", CFR_HEADER + cfr("a", 0, 0, 0))
    with pytest.raises(ValueError, match="more than one classifier output"):
        detect_samples(tmp_path)


# ------------------------------------------------------------------ 5: marker cleared before outputs change

def test_failed_rerun_with_other_parameters_cannot_leave_a_reusable_marker(tmp_path, monkeypatch):
    r1, r2 = fastq_pair(tmp_path, "s", "abcd")
    c = write(tmp_path / "s.cfr", CFR_HEADER + cfr("a", 562, 400, 60) + cfr("b", 562, 30, 60) + cfr("c", 0, 0, 0)
              + cfr("d", 562, 400, 60))
    out = tmp_path / "out"
    s = Sample("s", r1, r2, c, "centrifuger")
    extract_all([s], out, 0, 0)                                       # A: b (score 30) is classified
    a_r1 = fq_ids(out / "unclassified_reads" / "s_R1.fastq.gz")
    assert a_r1 == ["c"]

    class Boom:                                                       # B fails after its reads were written
        @staticmethod
        def dumps(*a, **k):
            raise RuntimeError("simulated failure after the reads were written")
        loads = staticmethod(json.loads)
    monkeypatch.setattr("unclassified_meta.pipeline.json", Boom)
    with pytest.raises(RuntimeError):
        extract_all([s], out, 100, 0)                                 # B: b falls below 100 -> b, c
    monkeypatch.undo()
    extract_all([s], out, 0, 0)                                       # A again must not reuse B's files
    assert fq_ids(out / "unclassified_reads" / "s_R1.fastq.gz") == a_r1


def test_merge_refuses_quant_tables_that_changed_or_were_not_published(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", in_frag=(20, 20), summary=str(summ))
    write(tmp_path / "s1.kreport", K2)
    write(tmp_path / "s2.kreport", K2)
    out = tmp_path / "taxa"
    taxa.build_table(taxa.find_reports(tmp_path, ".kreport"), out, "S", report_format="kraken2")
    taxa.merge_with_panel(out, pc, None, out / "m1")                  # sealed quant folder: accepted
    df = pd.read_csv(pc, sep="\t")
    df.loc[0, "s1"] += 1
    df.to_csv(pc, sep="\t", index=False)                              # changed after panel-quant wrote it
    with pytest.raises(ValueError, match="changed"):
        taxa.merge_with_panel(out, pc, None, out / "m2")
    (pc.parent / "provenance.json").unlink()                          # interrupted before publication
    with pytest.raises(ValueError, match="provenance"):
        taxa.merge_with_panel(out, pc, None, out / "m3")


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_interrupted_panel_quant_leaves_no_valid_provenance(tmp_path, monkeypatch):
    run, _ = tiny_run(tmp_path)
    p = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    fastq(tmp_path / "r.fastq.gz", "ab")
    samples = [Sample("s1", tmp_path / "r.fastq.gz")]
    q = tmp_path / "q"
    panel.quantify(p, samples, q, threads=1)
    assert (q / "provenance.json").exists()

    def boom(*a, **kw):
        raise StepError("simulated bowtie2 failure")
    monkeypatch.setattr(panel, "run_bowtie2", boom)
    with pytest.raises(StepError):
        panel.quantify(p, samples, q, threads=1, max_alignments=5)    # other options: cache miss, then failure
    assert not (q / "provenance.json").exists()


# ------------------------------------------------------------------ 6: databases bound to their source FASTA

def fake_makeblastdb(calls):
    def sh(cmd, *a, **kw):
        calls.append(cmd)
        out = cmd[cmd.index("-out") + 1]
        for ext in (".nhr", ".nin", ".nsq"):
            Path(out + ext).write_text("x")
    return sh


def test_blast_database_is_rebuilt_when_the_fasta_changes_in_the_middle_with_the_same_time(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(classify, "sh", fake_makeblastdb(calls))
    fa = tmp_path / "ref.fa"
    body = "A" * (3 << 20)
    fa.write_text(">r\n" + body + "\n")
    db1 = classify.blast_db(fa, tmp_path / "w" / "ref_db")
    t = fa.stat().st_mtime_ns
    fa.write_text(">r\n" + body[: 3 << 19] + "C" + body[(3 << 19) + 1:] + "\n")     # same size, middle changed
    os.utime(fa, ns=(t, t))
    classify.blast_db(fa, tmp_path / "w" / "ref_db")
    assert len(calls) == 2


def test_blast_database_with_a_mismatching_sidecar_or_a_partial_set_is_not_used(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(classify, "sh", fake_makeblastdb(calls))
    fa = write(tmp_path / "ref.fa", ">r\nACGT\n")
    db = classify.blast_db(fa, tmp_path / "w" / "ref_db")
    side = Path(str(db) + ".source.json")
    rec = json.loads(side.read_text())
    rec["sha256"] = "0" * 64
    side.write_text(json.dumps(rec))                                  # sidecar no longer matches the FASTA
    classify.blast_db(fa, tmp_path / "w" / "ref_db")
    assert len(calls) == 2
    Path(str(db) + ".nsq").unlink()                                   # incomplete database
    classify.blast_db(fa, tmp_path / "w" / "ref_db")
    assert len(calls) == 3
    for ext in (".nhr", ".nin", ".nsq"):                              # a database next to the FASTA without sidecar
        Path(str(fa) + ext).write_text("x")
    assert classify.blast_db(fa, tmp_path / "w" / "ref_db") != fa


def test_partial_diamond_database_is_not_reused(tmp_path, monkeypatch):
    built = []

    def fake_sh(cmd, *a, **kw):
        if cmd[:2] == ["diamond", "makedb"]:
            built.append(cmd)
            Path(cmd[cmd.index("-d") + 1] + ".dmnd").write_text("db")
        elif cmd[:2] == ["diamond", "blastx"]:
            Path(cmd[cmd.index("-o") + 1]).write_text("")
    monkeypatch.setattr(classify, "sh", fake_sh)
    monkeypatch.setattr(classify, "require", lambda *a: None)
    faa = write(tmp_path / "p.faa", ">p\nMKV\n")
    c = classify.Classifier({}, tmp_path / "w", 1)
    c._diamond(write(tmp_path / "q.fa", ">q\nACGT\n"), [str(faa), str(faa)], "prot", 1e-10)
    dmnd = next((tmp_path / "w").glob("prot*.dmnd"))
    dmnd.write_text("partial")                                       # interrupted write at the final name
    Path(str(dmnd) + ".source.json").unlink()
    c._diamond(write(tmp_path / "q.fa", ">q\nACGT\n"), [str(faa), str(faa)], "prot", 1e-10)
    assert len(built) == 2


# ------------------------------------------------------------------ 7: panel lineage also in --drop-lineage

def test_panel_lineage_listed_in_drop_lineage_is_still_checked_against_the_rank(tmp_path):
    (tmp_path / "reports").mkdir()
    write(tmp_path / "reports" / "s1.kreport", K2)
    write(tmp_path / "reports" / "s2.kreport", K2)
    pc = panel_dir(tmp_path, "all")
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", "taxa-table", "--input-dir", tmp_path / "reports",
                        "--suffix", ".kreport", "--report-format", "kraken2", "--rank", "G",
                        "--drop-lineage", "Escherichia coli", "--merge-panel", pc, "--panel-lineage",
                        "Escherichia coli", "-o", tmp_path / "t"], capture_output=True, text=True)
    assert p.returncode != 0 and "lies below rank G" in p.stdout + p.stderr


# ------------------------------------------------------------------ 8: requested reference genomes must exist

@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_missing_or_empty_reference_genome_is_an_error(tmp_path):
    from unclassified_meta import manifest, pipeline

    def run_with_virus(name, text):
        """The tiny run, classified with this virus file (recorded by the run, so panel-build accepts it)."""
        run, _ = tiny_run(tmp_path / name, extra_virus=False)
        vir = write(tmp_path / name / "virus.fa", text)
        cfg = yaml.safe_load((run / "config.used.yaml").read_text())
        cfg["references"]["virus"] = str(vir)
        (run / "config.used.yaml").write_text(yaml.safe_dump(cfg))
        (run / "references.used.json").write_text(json.dumps(
            {"refs_dir": None, "refs_version": "", "references": pipeline.reference_identity(cfg["references"])}))
        manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "test fixture"})
        return run
    run = run_with_virus("a", ">NC_9.1 other\nACGT\n>NC_1.1 Virus one\n\n")              # NC_1.1 has no sequence
    with pytest.raises(ValueError, match="NC_1.1"):
        panel.build_panel(run, tmp_path / "p", ["virus"], [], None, None, None, threads=1)
    run = run_with_virus("b", ">NC_9.1 other\nACGT\n")                                     # not in the file
    with pytest.raises(ValueError, match="NC_1.1"):
        panel.build_panel(run, tmp_path / "p2", ["virus"], [], None, None, None, threads=1)
    panel.build_panel(run, tmp_path / "p3", ["virus"], [], None, None, None, add_reference_genomes=False, threads=1)


# ------------------------------------------------------------------ 9 and 10: ANCOM-BC2 design checks (needs R)

@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_formula_variables_are_checked_and_samples_with_missing_values_are_reported(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    meta = pd.read_csv(tmp_path / "metadata.tsv", sep="\t", dtype=str)
    meta["batch"] = ["x", "y"] * 6
    meta.loc[2, "batch"] = ""                                         # sample '03' has no batch
    meta.to_csv(tmp_path / "metadata.tsv", sep="\t", index=False)
    base = ("--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
            tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A", "--rscript", rscript())
    rc, log = da_cli("--unsealed-input", *base, "--formula", "batch", "-o", tmp_path / "d1")
    assert rc != 0 and "group variable" in log
    rc, log = da_cli("--unsealed-input", *base, "--formula", "group + nosuchcolumn", "-o", tmp_path / "d2")
    assert rc != 0 and "nosuchcolumn" in log
    rc, log = da_cli("--unsealed-input", *base, "--formula", "group + batch", "-o", tmp_path / "d3")
    assert rc == 0, log
    header = [l for l in (tmp_path / "d3" / "ancombc2_results.tsv").read_text().splitlines() if l.startswith("#")]
    assert any("missing values" in l and "03" in l for l in header)
    assert any("samples per group: A=5, B=6" in l for l in header)


@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_struc_zero_does_not_call_a_feature_seen_in_one_sample_absent(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    c = pd.read_csv(tmp_path / "counts.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
    c.loc["t02", ["101", "102", "103", "104", "105"]] = 0             # t02 in 1 of 6 B samples
    c.loc["t02", "106"] = 3
    c.to_csv(tmp_path / "counts.tsv", sep="\t")
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv", "--metadata",
                     tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A", "--struc-zero",
                     "-o", tmp_path / "d", "--rscript", rscript())
    assert rc == 0, log
    sz = pd.read_csv(tmp_path / "d" / "structural_zeros.tsv", sep="\t")
    flagged = set(sz.loc[sz["absent_in_a_group"], "feature_id"])
    assert "t01" in flagged and "t02" not in flagged
    header = (tmp_path / "d" / "ancombc2_results.tsv").read_text()
    assert "neg_lb=FALSE" in header


# ------------------------------------------------------------------ 11: explicit report format, tree arithmetic

def test_report_format_is_explicit_and_quant_rounding_is_allowed_only_for_quant(tmp_path):
    rep = ("100.00\t10001\t1\tR\t1\troot\n"
           "50.00\t5001\t5001\tS\t11\t  Species one\n"
           "50.00\t5000\t5000\tS\t12\t  Species two\n")              # 1 + 5001 + 5000 = 10002: rounded rows
    f = write(tmp_path / "q.kreport", rep)
    with pytest.raises(TypeError):
        taxa.parse_report(f, "S")                                     # no guessing: the format is required
    r = taxa.parse_report(f, "S", report_format="centrifuger-quant")
    assert r.counts == {"11": 5001, "12": 5000} and r.ledger["rounding_residual"] == -1
    with pytest.raises(ValueError, match="clade"):
        taxa.parse_report(f, "S", report_format="kraken2")
    bad = K2.replace(" 60.00\t60\t6\tG\t561", " 60.00\t61\t6\tG\t561").replace(" 80.00\t80\t1\tR", " 80.00\t81\t1\tR") \
        .replace(" 70.00\t70\t4\tD", " 70.00\t71\t4\tD")              # genus clade 61 != own 6 + children 54
    with pytest.raises(ValueError, match="clade"):
        taxa.parse_report(write(tmp_path / "b.kreport", bad), "S", report_format="kraken2")


# ------------------------------------------------------------------ 12: reference versions and namespaces

def test_download_refs_keeps_the_version_of_a_reused_database(tmp_path, monkeypatch):
    vdir = tmp_path / "refs" / "viral"
    vdir.mkdir(parents=True)
    write(vdir / "refseq_viral.fna", ">v\nACGT\n")
    write(vdir / "refseq_viral_protein.dmnd", "db")
    write(tmp_path / "refs" / "REFS_VERSION.tsv",
          "reference\tversion\tsource\tfile\tn_sequences\tbytes\tdate\n"
          f"refseq_viral_genomic\tRefSeq release 200\tx\t{vdir / 'refseq_viral.fna'}\t1\t7\t2025-01-01\n"
          f"refseq_viral_protein\tRefSeq release 200\tx\t{vdir / 'refseq_viral_protein.dmnd'}\t\t2\t2025-01-01\n")
    from unclassified_meta import manifest
    manifest.seal(vdir, ["refseq_viral.fna", "refseq_viral_protein.dmnd"], {"step": "earlier download"})
    monkeypatch.setattr(refs, "fetch_text", lambda url: "238")
    monkeypatch.setattr(refs, "size_check", lambda *a, **k: None)
    monkeypatch.setattr(refs, "plan", lambda *a, **k: [])
    monkeypatch.setattr(refs, "require", lambda *a: None)
    monkeypatch.setattr(refs, "sh", lambda cmd, *a, **k: [Path(cmd[cmd.index("-out") + 1] + e).write_text("x")
                                                          for e in (".nhr", ".nin", ".nsq")])
    vt = refs.download_refs(tmp_path / "refs", only=("viral",))
    text = vt.read_text()
    assert "RefSeq release 200" in text and "RefSeq release 238" not in text and "2025-01-01" in text


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_panel_namespaces_come_from_the_references_the_run_used(tmp_path):
    run, _ = tiny_run(tmp_path, extra_virus=False)
    ra, rb = tmp_path / "refsA", tmp_path / "refsB"
    for d, rel in ((ra, 100), (rb, 200)):
        d.mkdir()
        write(d / "REFS_VERSION.tsv", f"reference\tversion\nrefseq_viral_genomic\tRefSeq release {rel}\n")
    cfg = yaml.safe_load((run / "config.used.yaml").read_text())
    cfg["references"]["refs_dir"] = str(ra)
    (run / "config.used.yaml").write_text(yaml.safe_dump(cfg))
    write(run / "references.used.json", json.dumps({"refs_dir": str(ra), "refs_version": (ra / "REFS_VERSION.tsv").read_text()}))
    from unclassified_meta import manifest
    manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "test fixture"})
    write(ra / "REFS_VERSION.tsv", "reference\tversion\nrefseq_viral_genomic\tRefSeq release 300\n")   # updated later
    p = panel.build_panel(run, tmp_path / "p", ["virus"], [], None, None, None, threads=1)
    assert pd.read_csv(p / "features.tsv", sep="\t")["namespace"].str.contains("refseq_viral_r100").any()
    with pytest.raises(ValueError, match="--refs"):
        panel.build_panel(run, tmp_path / "p2", ["virus"], [], None, None, str(rb), threads=1)


# ------------------------------------------------------------------ 13: identity per method

def test_reference_table_keeps_identity_per_method():
    contigs = pd.DataFrame([
        {"contig": "c1", "category": "virus", "best_hit": "X", "best_hit_title": "X", "method": "blastn", "identity": 80.0},
        {"contig": "c2", "category": "virus", "best_hit": "X", "best_hit_title": "X", "method": "tblastx_extra", "identity": 99.0},
        {"contig": "c3", "category": "virus", "best_hit": "X", "best_hit_title": "X", "method": "diamond_blastx", "identity": 60.0}])
    from collections import Counter
    t = summarize.reference_table(contigs, {"s": Counter({"c1": 10, "c2": 10, "c3": 10})}, {"s": 30})
    r = t.iloc[0]
    assert (r.blastn_identity, r.tblastx_identity, r.diamond_identity) == (80.0, 99.0, 60.0)
    assert "read_weighted_identity" not in t.columns


# ------------------------------------------------------------------ 14: pooled Sankey weights are scaled per sample

def test_pooled_reference_weights_are_scaled_to_each_samples_unclassified_reads():
    refs_t = pd.DataFrame([{"category": "virus", "best_hit": "X", "reads_pooled": 100, "pct_S1": 10.0, "pct_S2": 0.0},
                           {"category": "virus", "best_hit": "Y", "reads_pooled": 10, "pct_S1": 0.0, "pct_S2": 100.0}])
    inputs = pd.DataFrame({"sample": ["S1", "S2"], "unclassified_reads": [100000, 10], "reads_used": [1000, 10]})
    w = plots.reference_weights(refs_t, inputs)
    assert w.tolist() == [10000.0, 10.0]                             # not the 100:10 of the subsampled reads


# ------------------------------------------------------------------ extras: NA-like sample names stay text

def test_sample_names_like_NA_and_NULL_stay_text(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["NA", "NULL"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", samples=("NA", "NULL"), in_frag=(20, 20), summary=str(summ))
    out = tmp_path / "taxa"
    taxa.build_table({"NA": write(tmp_path / "a.kreport", K2), "NULL": write(tmp_path / "b.kreport", K2)}, out, "S",
                     report_format="kraken2")
    taxa.merge_with_panel(out, pc, None, out / "merged")
    rec = pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t", keep_default_na=False)
    assert rec["sample"].tolist() == ["NA", "NULL"] and (rec["status"] == "ok").all()
