"""Regression tests for the fifth code review (2026-11-18): indirect inputs, links in output trees, private
copies for DA, sealed extraction, step markers, README recipes, SILVA checksums, stale DA outputs."""
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

from unclassified_meta import control, da, io_utils, manifest, pipeline, refs, runner
from unclassified_meta.config import Sample

from test_ported_units import CFR_HEADER, K2, cfr, fastq, panel_dir, write
from test_review_fixes import tiny_run

HERE = Path(__file__).parent
README = HERE.parent / "README.md"


def cli(*args, cwd=None):
    p = subprocess.run([sys.executable, "-m", "unclassified_meta.cli", *map(str, args)], capture_output=True,
                       text=True, cwd=cwd)
    return p.returncode, p.stdout + p.stderr


# ------------------------------------------------------------------ 1: indirect inputs are checked before any write

def test_reports_listed_in_a_sheet_are_checked_before_the_log_is_opened(tmp_path):
    out = tmp_path / "taxa"
    out.mkdir()
    write(out / "s1.kreport", K2)                                     # a listed report inside the output folder
    write(tmp_path / "s2.kreport", K2)
    sheet = write(tmp_path / "reports.tsv", f"sample\treport\ns1\t{out / 's1.kreport'}\ns2\t{tmp_path / 's2.kreport'}\n")
    rc, log = cli("taxa-table", "--reports", sheet, "--report-format", "kraken2", "--rank", "S", "-o", out)
    assert rc != 0 and "inside" in log
    assert not (out / "taxa-table.log").exists() and not (out / ".unclassified_meta.lock").exists()


def test_references_inherited_from_the_run_are_checked_before_the_log_is_opened(tmp_path):
    run, hf = tiny_run(tmp_path)
    p = tmp_path / "p"
    p.mkdir()
    inside = p / "host.fa"
    inside.write_text(hf.read_text())                                # the run's host genome lies in the panel folder
    cfg = yaml.safe_load((run / "config.used.yaml").read_text())
    cfg["references"]["host"] = [{"name": "h", "fasta": str(inside)}]
    (run / "config.used.yaml").write_text(yaml.safe_dump(cfg))
    (run / "references.used.json").write_text(json.dumps(
        {"refs_dir": None, "refs_version": "", "references": pipeline.reference_identity(cfg["references"])}))
    manifest.seal(run, manifest.RUN_ARTIFACTS_PRESENT(run), {"step": "test fixture"})
    rc, log = cli("panel-build", "--run", run, "-o", p, "--categories", "virus")
    assert rc != 0 and "inside" in log
    assert not (p / "panel-build.log").exists()


# ------------------------------------------------------------------ 2: no links inside output trees

def test_a_symlinked_folder_inside_an_output_folder_is_refused(tmp_path):
    (tmp_path / "elsewhere").mkdir()
    out = tmp_path / "out"
    out.mkdir()
    os.symlink(tmp_path / "elsewhere", out / "logs")
    with pytest.raises(ValueError, match="link"):
        io_utils.preflight("extract", [], owned_dirs=[out])


def test_a_hard_linked_file_inside_an_output_folder_is_refused(tmp_path):
    other = write(tmp_path / "other.tsv", "x\n")
    out = tmp_path / "out"
    (out / "tables").mkdir(parents=True)
    os.link(other, out / "tables" / "input_summary.tsv")
    with pytest.raises(ValueError, match="link"):
        io_utils.preflight("extract", [], owned_dirs=[out])


def test_the_log_is_never_opened_through_a_link(tmp_path):
    target = write(tmp_path / "target.txt", "keep\n")
    (tmp_path / "logs").mkdir()
    os.symlink(target, tmp_path / "logs" / "pipeline.log")
    with pytest.raises((OSError, ValueError)):
        runner.setup_logging(tmp_path / "logs" / "pipeline.log")
    assert target.read_text() == "keep\n"


# ------------------------------------------------------------------ 3: DA reads private, verified copies

def fake_r(monkeypatch, write_zeros=True):
    seen = {}

    def run(cmd, capture_output=True, text=True):
        seen["cmd"] = list(map(str, cmd))
        o = Path(cmd[cmd.index("--outdir") + 1])
        (o / "ancombc2_results.tsv").write_text("# fake\nfeature_id\n")
        if write_zeros:
            (o / "structural_zeros.tsv").write_text("feature_id\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(da, "has_ancombc", lambda r: True)
    monkeypatch.setattr(da.subprocess, "run", run)
    return seen


def da_args(pc, meta, out):
    from unclassified_meta.cli import build_parser
    return build_parser().parse_args(["da", "--counts", str(pc), "--features", str(pc.parent / "features.tsv"),
                                      "--metadata", str(meta), "--group", "group", "-o", str(out)])


def test_r_reads_private_copies_of_the_verified_tables(tmp_path, monkeypatch):
    pc = panel_dir(tmp_path, "all")
    meta = write(tmp_path / "meta.tsv", "sample\tgroup\ns1\tA\ns2\tB\n")
    seen = fake_r(monkeypatch)
    out = tmp_path / "da"
    da.run_ancombc2(da_args(pc, meta, out))
    cmd = seen["cmd"]
    for flag, orig in (("--counts", pc), ("--features", pc.parent / "features.tsv"), ("--metadata", meta)):
        used = Path(cmd[cmd.index(flag) + 1])
        assert used != orig and out.resolve() in used.resolve().parents, (flag, used)


def test_a_rerun_of_da_removes_outputs_it_did_not_make_again(tmp_path, monkeypatch):
    pc = panel_dir(tmp_path, "all")
    meta = write(tmp_path / "meta.tsv", "sample\tgroup\ns1\tA\ns2\tB\n")
    out = tmp_path / "da"
    fake_r(monkeypatch, write_zeros=True)
    da.run_ancombc2(da_args(pc, meta, out))
    assert (out / "structural_zeros.tsv").exists()
    fake_r(monkeypatch, write_zeros=False)
    da.run_ancombc2(da_args(pc, meta, out))
    assert not (out / "structural_zeros.tsv").exists()
    manifest.verify(out)


# ------------------------------------------------------------------ 4: extracted reads are one sealed generation

def test_extracted_reads_are_one_sealed_generation(tmp_path, monkeypatch):
    for n in ("a", "b"):
        fastq(tmp_path / f"{n}.fastq.gz", "xyz")
    samples = [Sample(n, tmp_path / f"{n}.fastq.gz") for n in ("a", "b")]
    out = tmp_path / "res"
    pipeline.extract_all(samples, out)
    ur = out / "unclassified_reads"
    rec = manifest.verify(ur)
    assert {"a.fastq.gz", "b.fastq.gz"} <= set(rec["artifacts"])
    # a consumer refuses reads changed after the seal
    fastq(ur / "a.fastq.gz", "xy")
    with pytest.raises(ValueError, match="manifest|sealed"):
        control.make_controls([Sample("a", ur / "a.fastq.gz")], tmp_path / "ctrl", None, 11)
    # a rerun that fails midway leaves no consumable generation
    fastq(tmp_path / "b.fastq.gz", "xyzw")                           # b changes, so it is extracted again
    monkeypatch.setattr(pipeline, "extract_sample", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    with pytest.raises(RuntimeError):
        pipeline.extract_all(samples, out)
    with pytest.raises(ValueError, match="manifest"):
        control.make_controls([Sample("b", ur / "b.fastq.gz")], tmp_path / "ctrl2", None, 11)


# ------------------------------------------------------------------ 5: a marker must record the step's outputs

def test_a_step_marker_without_outputs_is_not_trusted(tmp_path):
    steps = runner.Steps(tmp_path)
    (steps.dir / "pool").write_text(json.dumps({"key": "k", "outputs": {}}))
    assert not steps.done("pool", "k")


# ------------------------------------------------------------------ 6: README recipes run as written

def readme_blocks(tag):
    text = README.read_text()
    blocks = re.findall(rf"<!-- readme-test: {tag} -->\s*```bash\n(.*?)```", text, flags=re.S)
    assert blocks, f"no README block tagged readme-test: {tag}"
    cmds = []
    for b in blocks:
        for line in re.sub(r"\\\n\s*", " ", b).splitlines():
            line = line.split(" #")[0].strip()
            if line and not line.startswith("#") and not line.startswith("micromamba"):
                cmds.append(line)
    return text, cmds


def run_readme(cmds, cwd, fake=None):
    for c in cmds:
        argv = shlex.split(c)
        if argv[0] == "unclassified-meta":
            rc, log = cli(*argv[1:], cwd=cwd)
            assert rc == 0, f"{c}\n{log}"
        else:
            fake(argv)


def synthetic(tmp_path):
    subprocess.run([sys.executable, str(HERE / "make_synthetic.py"), str(tmp_path)], check=True)
    return tmp_path


@pytest.mark.skipif(shutil.which("seqtk") is None, reason="seqtk not on PATH (activate the environment)")
def test_readme_control_recipe_runs_as_written(tmp_path):
    synthetic(tmp_path)
    text, cmds = readme_blocks("control")

    def centrifuger(argv):                     # the classifier itself is not run: every control read is unclassified
        r1 = Path(tmp_path / argv[argv.index("-1") + 1])
        with gzip.open(r1, "rt") as fh:
            ids = [l[1:].split()[0] for i, l in enumerate(fh) if i % 4 == 0]
        out = tmp_path / c_out
        out.write_text(CFR_HEADER + "".join(cfr(i, 0, 0, 0) for i in ids))
    c_out = re.search(r"> (\S+\.cfr\.tsv)", text).group(1)
    run_readme(cmds, tmp_path, centrifuger)
    for documented in re.findall(r"`(control_summary/[\w.]+\.(?:tsv|png))`", text):
        assert (tmp_path / documented).exists(), documented


E2E_TOOLS = ("spades.py", "megahit", "minimap2", "blastn", "diamond", "bowtie2", "samtools", "seqtk")


@pytest.mark.skipif(any(shutil.which(t) is None for t in E2E_TOOLS), reason="external tools not on PATH")
def test_readme_quick_start_runs_as_written(tmp_path):
    synthetic(tmp_path)
    shutil.copy(tmp_path / "host.fa", tmp_path / "genome.fa")
    _, cmds = readme_blocks("quickstart")
    run_readme(cmds, tmp_path)
    assert (tmp_path / "results" / "report.html").exists()


# ------------------------------------------------------------------ 7: SILVA must be verified

def test_silva_without_a_checksum_is_refused_unless_allowed(tmp_path, monkeypatch):
    def no_md5(url):
        raise OSError("offline")
    monkeypatch.setattr(refs, "fetch_text", no_md5)
    for name in ("require", "size_check", "bind_blast_db"):
        monkeypatch.setattr(refs, name, lambda *a, **k: None)
    monkeypatch.setattr(refs, "plan", lambda *a, **k: [])
    monkeypatch.setattr(refs, "db_files", lambda f: [])
    monkeypatch.setattr(refs, "download", lambda url, dest, **k: (dest.parent.mkdir(parents=True, exist_ok=True),
                                                                  dest.write_bytes(gzip.compress(b">r\nACGU\n")), dest)[2])
    with pytest.raises(Exception, match="--allow-unverified-silva"):
        refs.download_refs(tmp_path / "refs", only=("rrna",))
    vt = refs.download_refs(tmp_path / "refs2", only=("rrna",), allow_unverified_silva=True)
    assert "unverified" in vt.read_text()


def test_a_cached_download_is_checked_against_the_md5_before_reuse(tmp_path, monkeypatch):
    payload = gzip.compress(b">v\nACGT\n")
    calls = []

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    monkeypatch.setattr(refs.urllib.request, "urlopen", lambda req, timeout=0: (calls.append(1), Resp(payload))[1])
    dst = tmp_path / "src" / "SILVA.fasta.gz"
    good = hashlib.md5(payload).hexdigest()
    refs.download("https://example.org/SILVA.fasta.gz", dst, release="138.2", md5=good)
    refs.download("https://example.org/SILVA.fasta.gz", dst, release="138.2", md5=good)
    assert len(calls) == 1
    with pytest.raises(RuntimeError, match="checksum"):                # the provider's checksum changed
        refs.download("https://example.org/SILVA.fasta.gz", dst, release="138.2", md5="0" * 32)
    assert len(calls) == 2


# ------------------------------------------------------------------ final round: provenance of SILVA and summaries

def silva_env(monkeypatch, md5_reachable):
    def fetch(url):
        if not md5_reachable:
            raise OSError("offline")
        return "d41d8cd98f00b204e9800998ecf8427e  file\n"
    monkeypatch.setattr(refs, "fetch_text", fetch)
    for name in ("require", "size_check", "bind_blast_db"):
        monkeypatch.setattr(refs, name, lambda *a, **k: None)
    monkeypatch.setattr(refs, "plan", lambda *a, **k: [])
    monkeypatch.setattr(refs, "db_files", lambda f: [])
    monkeypatch.setattr(refs, "download", lambda url, dest, **k: (dest.parent.mkdir(parents=True, exist_ok=True),
                                                                  dest.write_bytes(gzip.compress(b">r\nACGU\n")), dest)[2])


def silva_status(d):
    return json.loads((d / "rrna" / manifest.NAME).read_text())["silva_md5"]


def rrna_version(vt):
    """The version column of the plant_rrna row of REFS_VERSION.tsv."""
    return [l.split("\t") for l in vt.read_text().splitlines() if l.startswith("plant_rrna\t")][0][1]


def test_an_unverified_silva_is_never_upgraded_to_verified_by_a_rerun(tmp_path, monkeypatch):
    d = tmp_path / "refs"
    silva_env(monkeypatch, md5_reachable=False)
    refs.download_refs(d, only=("rrna",), allow_unverified_silva=True)
    assert silva_status(d) == "unverified"
    # an ordinary rerun cannot verify the cached files: refused, and the record stays unverified
    with pytest.raises(RuntimeError, match="unverified"):
        refs.download_refs(d, only=("rrna",))
    assert silva_status(d) == "unverified"
    # with the opt-out the cached files are kept, still unverified (also in REFS_VERSION.tsv)
    vt = refs.download_refs(d, only=("rrna",), allow_unverified_silva=True)
    assert silva_status(d) == "unverified" and "unverified" in rrna_version(vt)
    # once the published checksum can be fetched, the files are fetched and checked again: verified
    silva_env(monkeypatch, md5_reachable=True)
    vt = refs.download_refs(d, only=("rrna",))
    assert silva_status(d) == "verified" and "unverified" not in rrna_version(vt)


def test_merge_refuses_an_input_summary_changed_after_panel_quant(tmp_path):
    from unclassified_meta import taxa
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", in_frag=(20, 20), summary=str(summ))
    with open(summ, "a") as fh:
        fh.write("s3\t1\t1\t0\n")                                    # changed after the panel was counted
    write(tmp_path / "s1.kreport", K2)
    write(tmp_path / "s2.kreport", K2)
    out = tmp_path / "taxa"
    taxa.build_table(taxa.find_reports(tmp_path, ".kreport"), out, "S", report_format="kraken2")
    with pytest.raises(ValueError, match="changed"):
        taxa.merge_with_panel(out, pc, None, out / "merged")
