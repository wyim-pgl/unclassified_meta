"""Generation manifests: a producing step writes its files first and a manifest with the sha256 of every file
last; a consuming step verifies the manifest before it reads anything. An interrupted or partly replaced
generation therefore has no valid manifest and is never used.

Producers: `run` (MANIFEST.json in the output folder), `panel-build` and `panel-quant` (provenance.json),
`taxa-table` (MANIFEST.json), `download-refs` (MANIFEST.json in each reference folder).
Consumers: panel-build (run), panel-quant (panel), taxa-table --merge-panel and da (quant / taxa tables),
run (references), download-refs (reuse of an earlier download).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import time
from pathlib import Path

NAME = "MANIFEST.json"
MARKER = ".unclassified_meta_generation"        # written FIRST by every producer; a folder with it needs a seal

# files of a finished `run` that later steps read
RUN_ARTIFACTS = ("contigs/contigs.fa", "tables/contig_category.tsv", "tables/contig_counts.tsv",
                 "tables/input_summary.tsv", "tables/read_category_by_sample.tsv", "tables/reference_by_sample.tsv",
                 "config.used.yaml", "references.used.json")


class ManifestError(ValueError):
    pass


def RUN_ARTIFACTS_PRESENT(run_dir: str | Path) -> list[str]:
    """The run artifacts that exist (a run without contigs has no contig tables)."""
    return [f for f in RUN_ARTIFACTS if (Path(run_dir) / f).exists()]


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def invalidate(folder: str | Path, name: str = NAME) -> None:
    """Remove the manifest BEFORE the files it covers are changed."""
    (Path(folder) / name).unlink(missing_ok=True)


def begin(folder: str | Path, producer: str, name: str = NAME) -> None:
    """Start a generation: write the producer marker (it stays) and remove the old manifest, before any other
    file of the folder changes. From then on the folder is only usable with a complete seal."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / MARKER).write_text(f"{producer}\n")
    invalidate(folder, name)


def is_generated(folder: str | Path) -> bool:
    """A folder written by this program (it has a producer marker or a manifest)."""
    folder = Path(folder)
    return (folder / MARKER).exists() or (folder / NAME).exists() or (folder / "provenance.json").exists()


def previous_artifacts(folder: str | Path, name: str = NAME) -> dict:
    """{file: sha256} of the folder's current manifest (read before a producer starts a new generation), so that
    stale outputs of the previous generation can be removed without touching files this program never wrote."""
    try:
        arts = json.loads((Path(folder) / name).read_text()).get("artifacts") or {}
        return {str(k): str(v) for k, v in arts.items()} if isinstance(arts, dict) else {}
    except (OSError, ValueError):
        return {}


STALE_NAME = re.compile(r"(imputations/(counts|ledger|ancombc2_results|structural_zeros)_\d{3,}\.tsv|collapsed\.tsv"
                        r"|mi_results\.tsv|mi_structural_zeros\.tsv|unique_counts\.tsv)")


def remove_stale(folder: str | Path, previous: dict, keep) -> None:
    """Delete outputs the previous generation sealed (imputation tables, collapsed.tsv, unique_counts.tsv,
    mi_results.tsv, mi_structural_zeros.tsv) that this one does not write
    again: only names this program writes, only regular files inside the folder (never through a link), and only
    when the file is still exactly the sealed one (a file replaced since is somebody else's and stays)."""
    folder = Path(folder)
    base = folder.resolve()
    for rel in sorted(set(previous) - set(keep)):
        if not STALE_NAME.fullmatch(rel):
            continue
        parts = Path(rel).parts
        f = folder / rel
        if f.is_symlink() or any((folder / Path(*parts[:i])).is_symlink() for i in range(1, len(parts))):
            continue
        if f.is_file() and base in f.resolve().parents and sha256(f) == previous[rel]:
            f.unlink()


def seal(folder: str | Path, files, record: dict | None = None, name: str = NAME) -> dict:
    """Write the manifest last: the record plus the sha256 of every listed file (paths relative to `folder`)."""
    from .io_utils import staged
    folder = Path(folder)
    files = list(files)
    if not files:
        raise ManifestError(f"nothing to seal in {folder}: a manifest must list its artifacts")
    rec = {**(record or {}), "sealed": time.strftime("%Y-%m-%d %H:%M:%S"),
           "artifacts": {str(f): sha256(folder / f) for f in files}}
    with staged(folder / name) as tmp:
        tmp.write_text(json.dumps(rec, indent=2) + "\n")
    return rec


def verify(folder: str | Path, name: str = NAME, required=(), what: str = "files") -> dict:
    """The manifest record when every listed file is present with the recorded checksum and every `required`
    file is listed; ManifestError otherwise."""
    folder = Path(folder)
    mf = folder / name
    if not mf.exists():
        raise ManifestError(f"{mf} not found (no manifest): the {what} in {folder} are not a complete, published generation "
                            "(an interrupted step, or made by an older version); run the step that makes them again")
    try:
        rec = json.loads(mf.read_text())
    except ValueError as e:
        raise ManifestError(f"{mf} is not a readable manifest ({e}); run the step again") from None
    arts = rec.get("artifacts") or {}
    if not arts:
        raise ManifestError(f"{mf} lists no artifacts: not a seal; run the step again")
    missing = [f for f in required if f not in arts]
    if missing:
        raise ManifestError(f"{mf} does not cover {', '.join(missing)}; run the step again")
    changed = sorted(f for f, h in arts.items() if not (folder / f).exists() or sha256(folder / f) != h)
    if changed:
        raise ManifestError(f"the {what} in {folder} changed after they were sealed by the manifest "
                            f"({', '.join(changed[:5])}); run the step that makes them again")
    return rec


def verify_file(path: str | Path, what: str = "files") -> dict:
    """The manifest record of the folder holding `path`, which must list exactly this file (by name) with an
    unchanged checksum: a consumer uses only files covered by a verified seal, never another file next to them."""
    path = Path(path)
    folder = path.parent
    for name in ("provenance.json", NAME):
        if (folder / name).exists():
            rec = verify(folder, name, what=what)
            if path.name not in (rec.get("artifacts") or {}):
                raise ManifestError(f"{path} is not covered by the manifest {folder / name}; use the sealed files "
                                    f"({', '.join(sorted(rec['artifacts']))})")
            return rec
    raise ManifestError(f"{folder} has no manifest (provenance.json / {NAME}): {path} is not part of a published "
                        "generation; use the tables written by panel-quant / taxa-table")


def private_copies(files: dict, dest: str | Path) -> dict:
    """Copy each input into `dest` (the consumer's own staging folder) and check every copy against the sha256
    its manifest recorded (None: an input without a manifest), so the consumer reads files that cannot change
    underneath it. Returns {original path (str): copy}."""
    import shutil
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    out = {}
    for i, (f, want) in enumerate(files.items()):
        copy = dest / f"{i}_{Path(f).name}"
        shutil.copyfile(f, copy)
        if want is not None and sha256(copy) != want:
            raise ManifestError(f"{f} changed after its manifest was verified (the copy's checksum differs); run the "
                                "step that makes it again")
        out[str(f)] = copy
    return out


def verify_consumed(paths, what: str = "files") -> None:
    """Inputs that this program wrote are used only as part of their sealed generation: for every input whose
    folder carries a producer marker or a manifest, that folder's manifest must verify and list the file. A
    folder with a marker but no manifest (a producer that started and did not finish) is refused. Each folder is
    verified once; files in folders this program did not write are not affected."""
    by_folder: dict[Path, list[Path]] = {}
    for p in paths:
        if p is not None:
            by_folder.setdefault(Path(p).absolute().parent, []).append(Path(p))
    for folder, files in by_folder.items():
        if not is_generated(folder):
            continue
        name = "provenance.json" if (folder / "provenance.json").exists() else NAME
        rec = verify(folder, name, what=what)
        missing = sorted({f.name for f in files} - set(rec.get("artifacts") or {}))
        if missing:
            raise ManifestError(f"{', '.join(missing)} in {folder} not covered by the manifest {folder / name}: not "
                                "part of the published generation; run the step that makes them again")
