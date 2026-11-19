"""Thin wrapper around the ANCOM-BC2 R script (r/ancombc2.R)."""
from __future__ import annotations

import os
import shutil
import json
import logging
import subprocess
from pathlib import Path

from .runner import StepError

log = logging.getLogger("unclassified_meta")
SCRIPT = Path(__file__).parent / "r" / "ancombc2.R"


def find_rscript(explicit: str | None = None) -> str:
    cand = explicit or os.environ.get("UNCLASSIFIED_META_RSCRIPT") or shutil.which("Rscript")
    if not cand or not (Path(cand).exists() or shutil.which(cand)):
        raise StepError("Rscript not found. Create the R environment (micromamba create -f environment-r.yml) and "
                        "pass --rscript <env>/bin/Rscript or set UNCLASSIFIED_META_RSCRIPT.")
    return cand


def has_ancombc(rscript: str) -> bool:
    p = subprocess.run([rscript, "-e", "suppressPackageStartupMessages(library(ANCOMBC))"],
                       capture_output=True, text=True)
    return p.returncode == 0


def verify_inputs(counts: str | Path, features: str | Path, unsealed_ok: bool = False) -> dict:
    """Tables made by panel-quant or taxa-table are used only as their sealed generation: each consumed file must
    be listed, with an unchanged checksum, in its folder's manifest (provenance.json of panel-quant, MANIFEST.json
    of taxa-table). A folder this program wrote (producer marker, ledger.tsv) without a complete seal is refused.
    Tables made by hand carry no manifest and are used only with --unsealed-input.
    Returns {file: sha256 recorded by its manifest, or None for a table used unsealed}."""
    from . import manifest
    expected = {}
    for f in (counts, features):
        folder = Path(f).parent
        expected[str(f)] = None
        if any((folder / n).exists() for n in ("provenance.json", manifest.NAME)):
            rec = manifest.verify_file(f, what="count tables")
            expected[str(f)] = rec["artifacts"][Path(f).name]
        elif manifest.is_generated(folder) or (folder / "ledger.tsv").exists():
            raise StepError(f"{folder} was written by panel-quant or taxa-table but has no complete manifest "
                            "(provenance.json / MANIFEST.json): the tables were not published completely; run that "
                            "step again")
        elif not unsealed_ok:
            raise StepError(f"{f} has no manifest (not made by panel-quant or taxa-table); pass --unsealed-input to "
                            "use a table made by hand as it is")
        else:
            log.warning("using %s without a manifest (--unsealed-input): no completeness or checksum check", f)
    return expected


# every file name the R script can write; a rerun removes the ones it did not write again
DA_OUTPUT_PATTERNS = ("ancombc2_results.tsv", "structural_zeros.tsv", "session_info.txt", "warnings.txt", "*.png",
                      "*.pdf")


def run_ancombc2(args, inputs=()) -> Path:
    """Run the R script into a private staging folder inside the output folder; only after it succeeded are its
    files checked against the inputs and moved into place, and the output folder is sealed with a manifest."""
    from . import manifest
    from .io_utils import make_staging_dir, refuse_overlap, remove_staging_dir
    expected = verify_inputs(args.counts, args.features, getattr(args, "unsealed_input", False))
    rscript = find_rscript(args.rscript)
    if not has_ancombc(rscript):
        raise StepError(f"{rscript} has no ANCOMBC package; use the Rscript of the environment-r.yml environment")
    outdir = Path(args.outdir)
    try:
        previous = set((json.loads((outdir / manifest.NAME).read_text()).get("artifacts") or {}))
    except (OSError, ValueError):
        previous = set()
    manifest.begin(outdir, "da")
    stage = make_staging_dir(outdir, "da")
    try:
        ledger = Path(args.counts).parent / "ledger.tsv"
        wanted = {**expected, str(args.metadata): None}
        if ledger.exists():
            wanted[str(ledger)] = (manifest.verify_file(ledger)["artifacts"].get("ledger.tsv")
                                   if expected[str(args.counts)] is not None else None)
        cp = manifest.private_copies(wanted, stage / "inputs")
        cmd = [rscript, str(SCRIPT), "--counts", cp[str(args.counts)], "--features", cp[str(args.features)],
               "--metadata", cp[str(args.metadata)], "--group", args.group, "--outdir", str(stage),
               "--prv-cut", str(args.prv_cut),
               "--lib-cut", str(args.lib_cut), "--alpha", str(args.alpha), "--p-adjust", args.p_adjust,
               "--threads", str(args.threads)]
        if ledger.exists():
            cmd += ["--ledger", str(cp[str(ledger)])]
        if args.formula:
            cmd += ["--formula", args.formula]
        if args.reference_level:
            cmd += ["--reference-level", args.reference_level]
        if args.struc_zero:
            cmd += ["--struc-zero"]
        if getattr(args, "neg_lb", False):
            cmd += ["--neg-lb"]
        if getattr(args, "rand_formula", None):
            cmd += ["--rand-formula", args.rand_formula]
        if getattr(args, "lme_relax_checks", False):
            cmd += ["--lme-relax-checks"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.stderr.strip():
            # R writes into the staging folder, which is removed below; name the published files instead
            print(p.stderr.strip().replace(str(stage), str(outdir)))
        if p.returncode != 0:
            raise StepError(f"ANCOM-BC2 failed (exit {p.returncode}); see the messages above")
        made = sorted(f.name for f in stage.iterdir() if f.is_file() and f.name != ".owned_by_unclassified_meta")
        refuse_overlap([outdir / f for f in made], [*inputs, args.counts, args.features, args.metadata], "da")
        for f in made:
            os.replace(stage / f, outdir / f)
        # outputs of an earlier run that this run did not write again (e.g. structural_zeros.tsv) are removed
        stale = previous | {f.name for pat in DA_OUTPUT_PATTERNS for f in outdir.glob(pat)}
        for name in sorted(stale - set(made)):
            f = outdir / name
            if f.parent == outdir and f.is_file() and not f.is_symlink():
                f.unlink()
        manifest.seal(outdir, made, {"step": "da", "counts": str(Path(args.counts).resolve()),
                                     "inputs_sha256": {k: manifest.sha256(v) for k, v in cp.items()}})
    finally:
        remove_staging_dir(stage)
    return outdir / "ancombc2_results.tsv"
