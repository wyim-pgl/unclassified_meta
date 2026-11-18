"""Run external commands with logging, and mark finished steps so a rerun resumes."""
from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("unclassified_meta")


class StepError(RuntimeError):
    pass


def setup_logging(logfile: Path) -> None:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    root = logging.getLogger("unclassified_meta")
    root.setLevel(logging.INFO)
    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(logfile.resolve())
               for h in root.handlers):
        fh = logging.FileHandler(logfile)
        fh.setFormatter(fmt)
        root.addHandler(fh)
    if not any(type(h) is logging.StreamHandler for h in root.handlers):
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        root.addHandler(sh)


def require(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise StepError(f"required tool(s) not on PATH: {', '.join(missing)} "
                        "(activate the micromamba environment from environment.yml)")


def sh(cmd: str | list[str], stdout: Path | None = None, cwd: Path | None = None) -> None:
    """Run a command (list = exec, str = bash -o pipefail). Raise StepError with the stderr tail on failure."""
    text = cmd if isinstance(cmd, str) else " ".join(shlex.quote(c) for c in cmd)
    log.info("RUN %s%s", text, f" > {stdout}" if stdout else "")
    args = ["bash", "-o", "pipefail", "-c", cmd] if isinstance(cmd, str) else cmd
    out = open(stdout, "w") if stdout else subprocess.PIPE
    try:
        proc = subprocess.run(args, stdout=out, stderr=subprocess.PIPE, text=True, cwd=cwd)
    finally:
        if stdout:
            out.close()
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or "").strip().splitlines()[-20:])
        raise StepError(f"command failed (exit {proc.returncode}): {text}\n{tail}")


class Steps:
    """Done-markers under <workdir>/.done so finished steps are skipped on a rerun."""

    def __init__(self, workdir: Path):
        self.dir = workdir / ".done"
        self.dir.mkdir(parents=True, exist_ok=True)

    def done(self, name: str) -> bool:
        return (self.dir / name).exists()

    def mark(self, name: str) -> None:
        (self.dir / name).touch()
