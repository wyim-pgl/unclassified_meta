"""De novo assembly of the pooled unclassified reads, with a peak-memory estimate before the run."""
from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

from .io_utils import read_fasta, write_fasta
from .runner import StepError, require

log = logging.getLogger("unclassified_meta")

# Rough peak resident memory: BASE_GB + PER_MILLION_GB * (million read pairs or single reads in the pool).
# Measured with /usr/bin/time -v on synthetic and real (Opuntia RNA-seq) inputs; see README "Memory and time".
MEMORY_MODEL = {
    "rnaspades": {"base_gb": 1.5, "per_million_gb": 1.5},
    "metaspades": {"base_gb": 2.0, "per_million_gb": 2.5},
    "megahit": {"base_gb": 0.5, "per_million_gb": 0.5},
}


def choose_assembler(cfg: dict) -> str:
    a = cfg["assembler"]
    if a == "auto":
        return "rnaspades" if cfg["mode"] == "rna" else "megahit"
    if a not in MEMORY_MODEL:
        raise ValueError(f"assembler must be auto, rnaspades, metaspades or megahit, not {a!r}")
    return a


def estimate_memory_gb(assembler: str, n_reads: int) -> float:
    m = MEMORY_MODEL[assembler]
    return round(m["base_gb"] + m["per_million_gb"] * n_reads / 1e6, 1)


def available_memory_gb() -> float | None:
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return None


def memory_warning(assembler: str, n_reads: int, memory_gb: float) -> str | None:
    """Return a warning text when the estimate exceeds the configured or available memory."""
    est = estimate_memory_gb(assembler, n_reads)
    avail = available_memory_gb()
    limit = min(memory_gb, avail) if avail else memory_gb
    if est <= limit:
        return None
    per_m = MEMORY_MODEL[assembler]["per_million_gb"]
    fit = max(int((limit - MEMORY_MODEL[assembler]["base_gb"]) / per_m * 1e6), 0)
    return (f"WARNING: {assembler} on {n_reads:,} pooled reads is estimated to need about {est} GB, "
            f"more than the {limit:.0f} GB available/configured (--memory). Lower the input with "
            f"--subsample (about {fit:,} reads in total fit; divide by the number of samples), "
            f"or raise --memory on a larger machine.")


def timed(cmd: list[str]) -> list[str]:
    """Prefix /usr/bin/time -v when present so the log records peak memory."""
    return (["/usr/bin/time", "-v"] + cmd) if os.path.exists("/usr/bin/time") else cmd


def run_timed(cmd: list[str], logfile: Path) -> dict:
    """Run `cmd`, keep its stdout/stderr in `logfile`, and return {'max_rss_gb', 'wall'} when measurable."""
    log.info("RUN %s", " ".join(cmd))
    with open(logfile, "w") as fo:
        proc = subprocess.run(timed(cmd), stdout=fo, stderr=subprocess.STDOUT, text=True)
    text = logfile.read_text(errors="replace")
    if proc.returncode != 0:
        tail = "\n".join(text.strip().splitlines()[-25:])
        raise StepError(f"assembly failed (exit {proc.returncode}); see {logfile}\n{tail}")
    out = {}
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["max_rss_gb"] = round(int(m.group(1)) / 1024 / 1024, 3)
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", text)
    if m:
        out["wall"] = m.group(1)
    return out


def assemble(assembler: str, pool: dict, outdir: Path, threads: int, memory_gb: float) -> tuple[Path, dict]:
    """pool = {'r1': path|None, 'r2': path|None, 's': path|None}. Returns (contig FASTA, resource dict)."""
    outdir.mkdir(parents=True, exist_ok=True)
    logfile = outdir / f"{assembler}.log"
    if assembler in ("rnaspades", "metaspades"):
        require("spades.py")
        target = outdir / "spades"
        if target.exists():
            shutil.rmtree(target)
        cmd = ["spades.py", "--rna" if assembler == "rnaspades" else "--meta",
               "-t", str(threads), "-m", str(max(int(memory_gb), 1)), "-o", str(target)]
        if pool.get("r1"):
            cmd += ["-1", str(pool["r1"]), "-2", str(pool["r2"])]
        if pool.get("s"):
            if assembler == "metaspades":
                raise StepError("metaSPAdes accepts one paired-end library only; use megahit for single-end reads")
            cmd += ["-s", str(pool["s"])]
        res = run_timed(cmd, logfile)
        fasta = target / ("transcripts.fasta" if assembler == "rnaspades" else "contigs.fasta")
    else:
        require("megahit")
        target = outdir / "megahit"
        if target.exists():
            shutil.rmtree(target)
        cmd = ["megahit", "-t", str(threads), "-m", str(int(memory_gb * 1e9)), "-o", str(target)]
        if pool.get("r1"):
            cmd += ["-1", str(pool["r1"]), "-2", str(pool["r2"])]
        if pool.get("s"):
            cmd += ["-r", str(pool["s"])]
        res = run_timed(cmd, logfile)
        fasta = target / "final.contigs.fa"
    if not fasta.exists():
        raise StepError(f"{assembler} finished but {fasta} is missing; see {logfile}")
    return fasta, res


def filter_and_rename(src: Path, dst: Path, map_tsv: Path, min_len: int) -> int:
    """Keep contigs >= min_len, rename to contig_000001 ... and write the old-to-new name map."""
    kept = []
    for name, seq in read_fasta(src):
        if len(seq) >= min_len:
            kept.append((name, seq))
    kept.sort(key=lambda x: -len(x[1]))
    with open(map_tsv, "w") as fo:
        fo.write("contig\tassembler_name\tlength\n")
        recs = []
        for i, (name, seq) in enumerate(kept, 1):
            new = f"contig_{i:06d}"
            fo.write(f"{new}\t{name}\t{len(seq)}\n")
            recs.append((new, seq))
    return write_fasta(recs, dst)
