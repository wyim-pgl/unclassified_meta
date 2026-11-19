"""Regression tests for the code-review findings of 2026-11-18 (counting, caching, panel identity, DA output)."""
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

from unclassified_meta import classify, da, extract, panel, pipeline, runner, taxa
from unclassified_meta.config import Sample, check_sample_names, validate_samples
from unclassified_meta.runner import StepError, file_id

from test_ported_units import CFR_HEADER, K2, NAMES, NODES, cfr, fastq, panel_dir, write

HAVE_PANEL_TOOLS = all(shutil.which(t) for t in ("bowtie2-build", "seqkit", "bowtie2"))


# ------------------------------------------------------------------ 1: a read ID must not reopen

def test_reopened_read_id_is_refused_by_extract_and_taxa_table(tmp_path):
    f = write(tmp_path / "s.cfr.tsv", CFR_HEADER + cfr("a", 562, 10, 10) + cfr("b", 562, 400, 60)
              + cfr("a", 562, 400, 60))                      # 'a' low score first, accepted later
    with pytest.raises(ValueError, match="separate blocks"):
        list(extract.centrifuge_reads(f, 49, 22))
    with pytest.raises(ValueError, match="separate blocks"):
        extract.parse_centrifuge(f, 49, 22)
    nodes, names = write(tmp_path / "nodes.dmp", NODES), write(tmp_path / "names.dmp", NAMES)
    with pytest.raises(ValueError, match="separate blocks"):
        taxa.build_table_from_per_read({"s": f}, tmp_path / "t", nodes, names, "S", (), "ns", 49, 22)
    k = write(tmp_path / "s.kraken", "U\tr1\t0\t100\t0:66\nC\tr2\t562\t100\t562:66\nC\tr1\t562\t100\t562:66\n")
    with pytest.raises(ValueError, match="occurs twice"):
        extract.parse_kraken2(k)


# ------------------------------------------------------------------ 2: reconciliation and report trees are exact

def test_reconciliation_does_not_accept_a_one_read_difference(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [21, 21]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", in_frag=(21, 21), summary=str(summ))    # 80 + 21 = 101 != 100
    write(tmp_path / "s1.kreport", K2)
    write(tmp_path / "s2.kreport", K2)
    out = tmp_path / "taxa"
    taxa.build_table(taxa.find_reports(tmp_path, ".kreport"), out, "S", report_format="kraken2")
    with pytest.raises(ValueError, match="do not add up"):
        taxa.merge_with_panel(out, pc, None, out / "merged")
    rec = pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t")
    assert (rec["difference"] == 1).all() and (rec["status"] == "MISMATCH").all()


def test_read_share_reports_must_add_up_exactly(tmp_path):
    broken = K2.replace(" 40.00\t40\t30\tS\t562", " 40.00\t40\t31\tS\t562")      # one read too many in a leaf
    broken = broken.replace(" 60.00\t60\t6\tG\t561", " 60.00\t60\t5\tG\t561")      # tree no longer sums
    broken = broken.replace("  6.00\t6\t6\tG\t590", "  6.00\t6\t5\tG\t590")
    r_ok = taxa.parse_report(write(tmp_path / "ok.kreport", K2), "S", report_format="kraken2")
    assert r_ok.ledger["report_kind"] == "kraken2" and r_ok.ledger["rounding_residual"] == 0
    with pytest.raises(ValueError, match="clade reads|does not add up"):
        taxa.parse_report(write(tmp_path / "b.kreport", broken), "S", report_format="kraken2")


# ------------------------------------------------------------------ 4 and 6: step markers and input identity

def test_file_id_tells_apart_same_name_size_and_second(tmp_path):
    a, b = tmp_path / "x" / "r.fq", tmp_path / "y" / "r.fq"
    a.parent.mkdir(), b.parent.mkdir()
    a.write_text("@1\nACGT\n+\nIIII\n")
    b.write_text("@1\nTTTT\n+\nIIII\n")
    os.utime(b, ns=(a.stat().st_atime_ns, a.stat().st_mtime_ns))
    assert file_id(a) != file_id(b)                                   # different path
    before = file_id(a)
    t = a.stat().st_mtime_ns
    a.write_text("@1\nGGGG\n+\nIIII\n")                                 # same size, replaced in place
    os.utime(a, ns=(t, t))
    assert file_id(a) != before                                       # the content digest notices it


# ------------------------------------------------------------------ 5: classification keyed on reference content

def test_replacing_a_reference_in_place_changes_the_classification_key_and_its_database(tmp_path, monkeypatch):
    ref = tmp_path / "virus.fa"
    ref.write_text(">v1\nACGTACGTAC\n")
    refs = {"refs_dir": None, "virus": str(ref), "host": [{"name": "h", "fasta": str(ref)}]}
    k1 = pipeline.reference_identity(refs)
    built = []

    def fake_sh(cmd, **kw):
        out = cmd[cmd.index("-out") + 1]
        built.append(out)
        for ext in (".nhr", ".nin", ".nsq"):
            Path(out + ext).write_text("x")
    monkeypatch.setattr(classify, "sh", fake_sh)
    classify.blast_db(ref, tmp_path / "virus_db")
    t = ref.stat().st_mtime_ns
    ref.write_text(">v1\nTTTTTTTTTT\n")
    os.utime(ref, ns=(t, t))
    assert pipeline.reference_identity(refs) != k1
    classify.blast_db(ref, tmp_path / "virus_db")
    assert len(built) == 2                                            # rebuilt for the new content


# ------------------------------------------------------------------ 7: translated identity never makes a bin 'same'

def test_tblastx_identity_is_not_mixed_into_nucleotide_identity():
    contigs = pd.DataFrame([
        {"contig": "c1", "length": 1000, "category": "virus", "layer": "virus", "method": "blastn_extra",
         "best_hit": "X1", "best_hit_title": "X1 Test virus", "identity": 90.0},
        {"contig": "c2", "length": 1000, "category": "virus", "layer": "virus", "method": "tblastx_extra",
         "best_hit": "X1", "best_hit_title": "X1 Test virus", "identity": 100.0}])
    feats, _ = panel.define_features(contigs, {"c1": 1, "c2": 99}, {"virus"}, {"extra_virus": "extra_virus:x"})
    row = feats.iloc[0]
    assert row.identity_class == "like" and row.weighted_identity == 90.0 and row.identity_max == 90.0


# ------------------------------------------------------------------ 10: a panel lineage below the rank is refused

def test_panel_lineage_below_the_rank_is_refused(tmp_path):
    nodes, names = write(tmp_path / "nodes.dmp", NODES), write(tmp_path / "names.dmp", NAMES)
    f = write(tmp_path / "s.cfr.tsv", CFR_HEADER + cfr("r1", 562, 400, 60) + cfr("r1", 564, 400, 60) + cfr("r2", 564, 400, 60))
    with pytest.raises(ValueError, match="lies below rank G"):
        taxa.build_table_from_per_read({"s": f}, tmp_path / "t", nodes, names, "G", (), "ns", panel_lineage=("562",))
    taxa.build_table_from_per_read({"s": f}, tmp_path / "t2", nodes, names, "S", (), "ns", panel_lineage=("562",))
    rep = write(tmp_path / "a.kreport", K2)
    with pytest.raises(ValueError, match="lies below rank G"):
        taxa.parse_report(rep, "G", panel_lineage=("Escherichia coli",), report_format="kraken2")
    assert taxa.parse_report(rep, "S", panel_lineage=("Escherichia coli",), report_format="kraken2").ledger["dropped_lineages"] == 40
    assert taxa.parse_report(rep, "G", drop=("Escherichia coli",), report_format="kraken2").ledger["dropped_lineages"] == 40   # plain drop


# ------------------------------------------------------------------ 14 and 15: sample names

def test_numeric_looking_sample_names_stay_text_in_the_merge(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["01", "101"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", samples=("01", "101"), in_frag=(20, 20), summary=str(summ))
    out = tmp_path / "taxa"
    taxa.build_table({"01": write(tmp_path / "a.kreport", K2), "101": write(tmp_path / "b.kreport", K2)}, out, "S", report_format="kraken2")
    taxa.merge_with_panel(out, pc, None, out / "merged")
    rec = pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t", dtype={"sample": str})
    assert rec["sample"].tolist() == ["01", "101"] and (rec["status"] == "ok").all()


def test_reserved_sample_names_are_refused(tmp_path):
    with pytest.raises(ValueError, match="reserved"):
        check_sample_names(["s1", "feature_id"])
    with pytest.raises(ValueError, match="duplicate"):
        check_sample_names(["s1", "s1"])
    fastq(tmp_path / "f.fastq.gz", "ab")
    with pytest.raises(ValueError, match="reserved"):
        validate_samples([Sample("feature_id", tmp_path / "f.fastq.gz")])
    with pytest.raises(ValueError, match="reserved"):
        taxa.build_table({"label": write(tmp_path / "a.kreport", K2)}, tmp_path / "t", "S", report_format="kraken2")


# ------------------------------------------------------------------ 3, 11, 12: panel-build (needs bowtie2, seqkit)

def tiny_run(tmp_path, extra_virus=True):
    """A finished-run folder: two virus contigs (one BLASTN hit on a RefSeq genome, one on an extra virus FASTA)."""
    import random
    rng = random.Random(1)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    run = tmp_path / "run"
    (run / "tables").mkdir(parents=True)
    (run / "contigs").mkdir()
    g1, g2, host = seq(3000), seq(3000), seq(5000)
    (run / "contigs" / "contigs.fa").write_text(f">c1\n{g1[:2000]}\n>c2\n{g2[:2000]}\n")
    pd.DataFrame([
        {"contig": "c1", "length": 2000, "category": "virus", "layer": "virus", "method": "blastn",
         "best_hit": "NC_1.1", "best_hit_title": "NC_1.1 Virus one, complete genome", "identity": 80.0, "cover": 1.0},
        {"contig": "c2", "length": 2000, "category": "virus", "layer": "virus", "method": "blastn_extra",
         "best_hit": "X2", "best_hit_title": "X2 Extra virus two", "identity": 80.0, "cover": 1.0},
    ]).to_csv(run / "tables" / "contig_category.tsv", sep="\t", index=False)
    vir = write(tmp_path / "refseq_viral.fa", f">NC_1.1 Virus one\n{g1}\n>NC_9.1 other\n{seq(500)}\n")
    xv = write(tmp_path / "extra_virus.fa", f">X2 Extra virus two\n{g2}\n")
    hf = write(tmp_path / "host.fa", f">chr1\n{host}\n")
    refs = {"virus": str(vir), "host": [{"name": "h", "fasta": str(hf)}]}
    if extra_virus:
        refs["extra_virus_fasta"] = str(xv)
    (run / "config.used.yaml").write_text(yaml.safe_dump({"references": refs}))
    from unclassified_meta import manifest, pipeline
    (run / "references.used.json").write_text(json.dumps(
        {"refs_dir": None, "refs_version": "", "references": pipeline.reference_identity(refs)}))
    manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "test fixture"})     # a finished run
    return run, hf


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_extra_virus_genomes_are_competitors_and_missing_references_are_errors(tmp_path):
    run, _ = tiny_run(tmp_path)
    out = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    s2f = pd.read_csv(out / "seq2feature.tsv", sep="\t").set_index("sequence")["feature_or_decoy"]
    assert s2f["ref|NC_1.1"].endswith("|NC_1.1|reference")
    assert s2f["ref|extra_virus|X2"] == "extra_virus:extra_virus.fa|X2|reference"
    prov = json.loads((out / "provenance.json").read_text())
    assert {"virus_reference", "extra_virus_reference", "host:h"} <= set(prov["inputs"])
    with pytest.raises(FileNotFoundError, match="rRNA reference"):
        panel.build_panel(run, tmp_path / "p2", ["virus"], [], str(tmp_path / "missing_rrna.fa"), None, None, threads=1)
    with pytest.raises(FileNotFoundError, match="host genome"):
        panel.build_panel(run, tmp_path / "p3", ["virus"], [{"name": "x", "fasta": str(tmp_path / "nope.fa")}],
                          None, None, None, threads=1)
    with pytest.raises(FileNotFoundError, match="virus reference"):
        panel.build_panel(run, tmp_path / "p4", ["virus"], [], None, str(tmp_path / "nope_virus.fa"), None, threads=1)


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_interrupted_or_altered_panel_is_never_used(tmp_path, monkeypatch):
    run, _ = tiny_run(tmp_path)
    out = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    fastq(tmp_path / "r.fastq.gz", "ab")
    samples = [Sample("s1", tmp_path / "r.fastq.gz")]
    assert panel.verify_panel(out)["panel_id"]
    (out / "panel.fa").write_text((out / "panel.fa").read_text() + ">extra\nACGT\n")      # altered after build
    with pytest.raises(StepError, match="changed since panel-build"):
        panel.quantify(out, samples, tmp_path / "q", threads=1)
    real_sh = panel.sh

    def failing_sh(cmd, *a, **kw):
        if "bowtie2-build" in str(cmd):
            raise StepError("simulated failure")
        return real_sh(cmd, *a, **kw)
    monkeypatch.setattr(panel, "sh", failing_sh)
    with pytest.raises(StepError, match="simulated failure"):
        panel.build_panel(run, out, ["virus"], [], None, None, None, threads=1)
    assert not (out / "provenance.json").exists()
    with pytest.raises(StepError, match="incomplete"):
        panel.quantify(out, samples, tmp_path / "q2", threads=1)


# ------------------------------------------------------------------ 8 and 13: ANCOM-BC2 output (needs R)

def rscript():
    try:
        r = da.find_rscript(None)
    except Exception:
        return None
    return r if da.has_ancombc(r) else None


def da_cli(*args):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", "da", *map(str, args)], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def counts_with_a_structural_zero(d: Path, low_sample_depth: int | None = None):
    import random
    rng = random.Random(3)
    samples = [f"{i:02d}" for i in range(1, 7)] + [f"1{i:02d}" for i in range(1, 7)]   # '01'..'06', '101'..'106'
    group = ["A"] * 6 + ["B"] * 6
    feats = [f"t{i:02d}" for i in range(1, 61)]
    rows = []
    for f in feats:
        row = []
        for s, g in zip(samples, group):
            lam = 200 * (0.5 + rng.random())
            if f == "t01" and g == "B":
                lam = 0                                          # absent from every B sample
            if low_sample_depth is not None and s == "01":
                lam = low_sample_depth / len(feats)
            row.append(int(lam))
        rows.append(row)
    pd.DataFrame(rows, columns=samples, index=pd.Index(feats, name="feature_id")).to_csv(d / "counts.tsv", sep="\t")
    pd.DataFrame({"feature_id": feats, "label": feats}).to_csv(d / "features.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": samples, "group": group}).to_csv(d / "metadata.tsv", sep="\t", index=False)


@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_structural_zeros_are_exported_and_numeric_sample_names_are_matched_by_name(tmp_path):
    counts_with_a_structural_zero(tmp_path)
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv",
                     "--metadata", tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A",
                     "--struc-zero", "-o", tmp_path / "da", "--rscript", rscript())
    assert rc == 0, log
    sz = pd.read_csv(tmp_path / "da" / "structural_zeros.tsv", sep="\t")
    assert "t01" in set(sz.loc[sz["absent_in_a_group"], "feature_id"])
    header = [l for l in (tmp_path / "da" / "ancombc2_results.tsv").read_text().splitlines() if l.startswith("#")]
    assert any("structural zeros: 1 feature(s)" in l for l in header)
    assert any("retained after ANCOM-BC2's filters: A=6, B=6" in l for l in header)


@pytest.mark.skipif(rscript() is None, reason="R environment with ANCOMBC not found")
def test_design_is_reported_and_checked_after_the_library_size_filter(tmp_path):
    counts_with_a_structural_zero(tmp_path, low_sample_depth=60)       # sample '01' has ~60 reads
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv",
                     "--metadata", tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A",
                     "--lib-cut", 1000, "-o", tmp_path / "da", "--rscript", rscript())
    assert rc == 0, log
    header = [l for l in (tmp_path / "da" / "ancombc2_results.tsv").read_text().splitlines() if l.startswith("#")]
    assert any("samples per group: A=6, B=6; retained after ANCOM-BC2's filters: A=5, B=6" in l for l in header)
    assert any("removed by ANCOM-BC2's filters" in l and "01" in l for l in header)
    # a library-size cut that leaves fewer than 2 samples in a group is refused
    c = pd.read_csv(tmp_path / "counts.tsv", sep="\t", dtype={"feature_id": str}).set_index("feature_id")
    c[["02", "03", "04", "05"]] = 1
    c.to_csv(tmp_path / "counts.tsv", sep="\t")
    rc, log = da_cli("--unsealed-input", "--counts", tmp_path / "counts.tsv", "--features", tmp_path / "features.tsv",
                     "--metadata", tmp_path / "metadata.tsv", "--group", "group", "--reference-level", "A",
                     "--lib-cut", 1000, "-o", tmp_path / "da2", "--rscript", rscript())
    assert rc != 0 and "at least 2 retained samples" in log
