"""Reversed-read negative control for read classifiers (Centrifuge, Centrifuger; Kraken2 in part).

A read whose sequence is reversed (NOT reverse-complemented) keeps its length and base composition but matches
no real genome, so every classification of such a read is a chance hit. Classifying control reads with the same
index and options as the real reads shows which --min-score / --min-hitlen keeps chance hits rare.

  control-reads    writes the control FASTQ files (IDs prefixed CTRL_)
  (the user runs the classifier on them, alone or together with the real reads)
  control-summary  reads the per-read output and reports the chance-hit rate over a grid of thresholds

Origin: the reversed-read control and the score calibration are ideas contributed by Eli Holschbach.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from bisect import bisect_right
from pathlib import Path

import pandas as pd

from .config import Sample
from .io_utils import normalise_read_id, open_text, refuse_overlap, staged
from .runner import StepError

log = logging.getLogger("unclassified_meta")
PREFIX = "CTRL_"
SCORE_GRID = (0, 16, 25, 36, 49, 64, 81, 100, 150, 200, 300)
HITLEN_GRID = (0,)


def _records(stream):
    while True:
        head = stream.readline()
        if not head:
            return
        seq, plus, qual = stream.readline(), stream.readline(), stream.readline()
        if not qual:
            raise ValueError(f"truncated FASTQ record at {head.strip()!r}")
        yield head.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n")


def reverse_fastq(src: str | Path, dst: str | Path, n_reads: int | None, seed: int, prefix: str = PREFIX) -> int:
    """Write reversed copies of `n_reads` reads sampled from `src` (all reads when n_reads is None). The sampler
    is waited for and checked INSIDE the staging context: a failed or truncated sampling raises before the
    rename, so the previous output (if any) stays and no partial file is published."""
    refuse_overlap([dst], [src], "control-reads")
    if n_reads and shutil.which("seqtk") is None:
        raise StepError("seqtk not on PATH (activate the environment), or use --n-reads all")
    n = 0
    with staged(dst) as tmp, open_text(tmp, "wt") as fo:
        if n_reads:
            proc = subprocess.Popen(["seqtk", "sample", "-s", str(seed), str(src), str(n_reads)],
                                    stdout=subprocess.PIPE, text=True)
            stream = proc.stdout
        else:
            proc, stream = None, open_text(src)
        try:
            for head, seq, qual in _records(stream):
                fo.write(f"@{prefix}{head[1:]}\n{seq[::-1]}\n+\n{qual[::-1]}\n")
                n += 1
        finally:
            if proc is None:
                stream.close()
        if proc is not None and proc.wait() != 0:
            raise StepError(f"seqtk sample failed on {src}")
    return n


def make_controls(samples: list[Sample], outdir: str | Path, n_reads: int | None = 100_000, seed: int = 11,
                  prefix: str = PREFIX) -> pd.DataFrame:
    """Control FASTQ(s) per sample in <outdir>/control_reads/ and a sample sheet for the classifier run."""
    outdir = Path(outdir)
    d = outdir / "control_reads"
    # before anything is written: no control file may be (or link to) any sample's input
    outs = [outdir / "control_samples.tsv"]
    for s in samples:
        outs += [d / f"{s.name}_ctrl_R1.fastq.gz", d / f"{s.name}_ctrl_R2.fastq.gz"] if s.paired else [d / f"{s.name}_ctrl.fastq.gz"]
    refuse_overlap(outs, [f for s in samples for f in (s.r1, s.r2)], "control-reads")
    from . import manifest
    manifest.verify_consumed([f for s in samples for f in (s.r1, s.r2)], "reads")
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in samples:
        if s.paired:
            o1, o2 = d / f"{s.name}_ctrl_R1.fastq.gz", d / f"{s.name}_ctrl_R2.fastq.gz"
            n1 = reverse_fastq(s.r1, o1, n_reads, seed, prefix)
            n2 = reverse_fastq(s.r2, o2, n_reads, seed, prefix)     # same seed keeps the mates paired
            if n1 != n2:
                raise ValueError(f"{s.name}: {n1} mate-1 and {n2} mate-2 control reads; are R1/R2 in the same order?")
        else:
            o1, o2 = d / f"{s.name}_ctrl.fastq.gz", None
            n1 = reverse_fastq(s.r1, o1, n_reads, seed, prefix)
        log.info("control reads: %s, %d %s", s.name, n1, "pairs" if s.paired else "reads")
        rows.append({"sample": f"{s.name}_ctrl", "r1": str(o1.resolve()), "r2": str(o2.resolve()) if o2 else ".",
                     "control_reads": n1})
    df = pd.DataFrame(rows)
    with staged(outdir / "control_samples.tsv") as tmp:
        df.to_csv(tmp, sep="\t", index=False)
    return df


class _Grid:
    """Counts of reads that stay classified at every (min_hitlen, min_score) cell."""

    def __init__(self, scores, hitlens):
        self.scores, self.hitlens = sorted(set(scores)), sorted(set(hitlens))
        self.cell = [[0] * (len(self.scores) + 1) for _ in self.hitlens]
        self.n = 0

    def add(self, hits: list) -> None:
        """hits = [(score, hitlen), ...] of one read (empty = no assignment)."""
        self.n += 1
        for i, h in enumerate(self.hitlens):
            best = max((s for s, hl in hits if hl >= h), default=None)
            if best is not None:
                self.cell[i][bisect_right(self.scores, best)] += 1      # classified for every score <= best

    def classified(self, i: int, j: int) -> int:
        return sum(self.cell[i][j + 1:])


def _read_hits(path: str | Path, classifier: str):
    """Yield (read ID, [(score, hitLength), ...]) per read; consecutive lines of one read are grouped."""
    current, hits = None, []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if classifier == "kraken2":
                if len(f) < 3 or f[0] not in ("C", "U"):
                    raise ValueError(f"{path}: not Kraken2 per-read output (line: {line[:80]!r})")
                yield normalise_read_id(f[1]), ([(0, 0)] if f[0] == "C" else [])
                continue
            if f[0] == "readID":
                continue
            if len(f) < 6:
                raise ValueError(f"{path}: Centrifuge(r) per-read output needs the score and hitLength columns "
                                 f"(line: {line[:80]!r})")
            rid = normalise_read_id(f[0])
            if rid != current:
                if current is not None:
                    yield current, hits
                current, hits = rid, []
            if f[1] != "unclassified" and f[2].strip() not in ("0", ""):
                hits.append((int(float(f[3])), int(float(f[5]))))
    if current is not None:
        yield current, hits


def summarise(per_read: list, outdir: str | Path, classifier: str = "centrifuger", scores=SCORE_GRID,
              hitlens=HITLEN_GRID, target_pct: float = 1.0, prefix: str = PREFIX, all_control: bool = False) -> pd.DataFrame:
    """Chance-hit rate of the control reads (and the classified share of real reads, if present) per threshold."""
    if classifier not in ("kraken2", "centrifuge", "centrifuger"):
        raise ValueError(f"unknown classifier {classifier!r}")
    if classifier == "kraken2":
        scores, hitlens = (0,), (0,)        # Kraken2 output has no score or hit length: one cell
    if any(v < 0 for v in (*scores, *hitlens)) or not 0 <= target_pct <= 100:
        raise ValueError("thresholds must be 0 or larger and --target-pct between 0 and 100")
    ctrl, real = _Grid(scores, hitlens), _Grid(scores, hitlens)
    seen: dict[str, set] = {"ctrl": set(), "real": set()}
    for path in per_read:
        for rid, hits in _read_hits(path, classifier):
            is_ctrl = all_control or rid.startswith(prefix)
            key = "ctrl" if is_ctrl else "real"
            if rid in seen[key]:
                raise ValueError(f"{path}: read {rid} occurs in separated blocks or in two files; the per-read "
                                 "output must keep the lines of a read together and each read once")
            seen[key].add(rid)
            (ctrl if is_ctrl else real).add(hits)
    if ctrl.n == 0:
        raise ValueError(f"no control reads found (read IDs starting with {prefix!r}). Classify the files made by "
                         "`control-reads`, or pass --all-control when every read in the file is a control read.")
    rows = []
    for i, h in enumerate(ctrl.hitlens):
        for j, sc in enumerate(ctrl.scores):
            c = 100 * ctrl.classified(i, j) / ctrl.n
            row = {"min_hitlen": h, "min_score": sc, "control_reads": ctrl.n,
                   "control_classified": ctrl.classified(i, j), "control_classified_pct": round(c, 4)}
            if real.n:
                r = 100 * real.classified(i, j) / real.n
                row.update({"real_reads": real.n, "real_classified": real.classified(i, j),
                            "real_classified_pct": round(r, 4),
                            # share of real reads classified beyond chance: (R - C) / (1 - C)
                            "real_beyond_chance_pct": round(100 * (r - c) / (100 - c), 4) if c < 100 else 0.0})
            row["meets_target"] = c <= target_pct
            rows.append(row)
    df = pd.DataFrame(rows)
    ok = df[df["meets_target"]]
    df["recommended"] = False
    if len(ok):
        # the least strict cell that meets the target: lowest score, then shortest hit length
        best = ok.sort_values(["min_score", "min_hitlen"]).index[0]
        df.loc[best, "recommended"] = True
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / "control_summary.tsv", sep="\t", index=False)
    _plot(df, outdir, target_pct, classifier)
    return df


def _plot(df: pd.DataFrame, outdir: Path, target_pct: float, classifier: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_real = "real_classified_pct" in df.columns
    fig, axes = plt.subplots(1, 2 if has_real else 1, figsize=(9 if has_real else 4.8, 3.6), squeeze=False)
    floor = 100 / max(int(df["control_reads"].iloc[0]), 1) / 2          # half a read: what 0 % is drawn as
    for h, g in df.groupby("min_hitlen"):
        g = g.sort_values("min_score")
        axes[0][0].plot(g["min_score"], g["control_classified_pct"].clip(lower=floor), marker="o", ms=3,
                        label=f"min hit length {h}")
        if has_real:
            axes[0][1].plot(g["min_score"], g["real_classified_pct"], marker="o", ms=3, label=f"min hit length {h}")
    axes[0][0].axhline(target_pct, color="grey", ls="--", lw=0.8)
    axes[0][0].set_yscale("log")
    axes[0][0].set_ylabel("control reads classified (%)\n= chance-hit rate")
    axes[0][0].set_title("Reversed control reads", fontsize=9)
    if has_real:
        axes[0][1].set_ylabel("real reads classified (%)")
        axes[0][1].set_title("Real reads", fontsize=9)
    for ax in axes[0]:
        ax.set_xlabel("minimum score" if classifier != "kraken2" else "Kraken2 (no score threshold)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"control_summary.{ext}", dpi=200)
    plt.close(fig)
