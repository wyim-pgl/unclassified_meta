"""Defects found in the P257 real-data test (2026-11-18)."""
import os
import stat
import subprocess
from pathlib import Path

import pytest

from unclassified_meta import da, manifest
from unclassified_meta.io_utils import staged

from test_review_round5 import da_args, panel_dir, write


@pytest.fixture
def umask():
    def set_(mask):
        old = os.umask(mask)
        restore.append(old)
    restore = []
    yield set_
    for old in restore[:1]:
        os.umask(old)


def mode(p):
    return stat.S_IMODE(Path(p).stat().st_mode)


@pytest.mark.parametrize("mask, expected", [(0o022, 0o644), (0o002, 0o664), (0o077, 0o600)])
def test_staged_outputs_follow_the_umask_like_any_other_file(tmp_path, umask, mask, expected):
    umask(mask)
    with staged(tmp_path / "table.tsv") as tmp:
        tmp.write_text("a\n")
    assert mode(tmp_path / "table.tsv") == expected


def test_a_sealed_manifest_is_as_readable_as_the_artifacts_it_lists(tmp_path, umask):
    umask(0o022)
    (tmp_path / "counts.tsv").write_text("x\n")
    manifest.seal(tmp_path, ["counts.tsv"], {"step": "test"})
    assert mode(tmp_path / manifest.NAME) == mode(tmp_path / "counts.tsv") == 0o644


def test_da_reports_the_published_results_path_not_the_removed_staging_folder(tmp_path, monkeypatch, capsys):
    pc = panel_dir(tmp_path, "all")
    meta = write(tmp_path / "meta.tsv", "sample\tgroup\ns1\tA\ns2\tB\n")
    out = tmp_path / "da"

    def run(cmd, capture_output=True, text=True):
        o = Path(cmd[cmd.index("--outdir") + 1])
        (o / "ancombc2_results.tsv").write_text("# fake\nfeature_id\n")
        return subprocess.CompletedProcess(cmd, 0, "", f"fit done\nresults: {o / 'ancombc2_results.tsv'}\n")
    monkeypatch.setattr(da, "find_rscript", lambda x: "Rscript")
    monkeypatch.setattr(da, "has_ancombc", lambda r: True)
    monkeypatch.setattr(da.subprocess, "run", run)
    da.run_ancombc2(da_args(pc, meta, out))
    printed = capsys.readouterr().out
    assert ".staging-" not in printed
    assert f"results: {out / 'ancombc2_results.tsv'}" in printed
