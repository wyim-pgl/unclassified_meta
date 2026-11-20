"""Repeat the four CLI commands on the same input files using an immutable local source export."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "pinned/src"


def retry(name):
    folder = ROOT / "runs" / name
    rec = json.loads((folder / "run.json").read_text())
    out = folder / "pinned_fit"
    out.mkdir()
    cli = [sys.executable, "-m", "unclassified_meta.cli"]
    env = dict(os.environ, PYTHONPATH=str(SOURCE))
    tables = ["taxa-table", "--reports", "reports.tsv", "--per-read", "--nodes", "../../nodes.dmp",
              "--names", "../../names.dmp", "--taxonomy-name", "sim", "--rank", "S"]
    commands = [("table_lca", [*tables, "--assign", "lca", "-o", "pinned_fit/lca"]),
                ("table_em", [*tables, "--assign", "em", "--imputations", "5", "--seed", str(rec["seed"]), "-o", "pinned_fit/em"])]
    for mode in ["lca", "em"]:
        commands.append(("da_" + mode, ["da", "--counts", f"pinned_fit/{mode}/counts.tsv", "--features",
                         f"pinned_fit/{mode}/features.tsv", "--metadata", "meta.tsv", "--group", "group",
                         "--formula", "group", "--reference-level", "control", "--p-adjust", "BH", "--alpha",
                         "0.05", "--threads", "1", "-o", f"pinned_fit/da_{mode}", *(["--imputations"] if mode == "em" else [])]))
    start = time.monotonic()
    record = dict(started_utc=dt.datetime.now(dt.timezone.utc).isoformat(), commands=[], status="ok")
    for label, args in commands:
        t = time.monotonic()
        with (out / f"{label}.log").open("w") as f:
            p = subprocess.run(cli+args, cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT)
        record["commands"].append(dict(label=label, command=shlex.join(cli+args), seconds=time.monotonic()-t, exit=p.returncode))
        if p.returncode:
            record["status"] = "failed"
            break
    record["elapsed_s"] = time.monotonic()-start
    (out / "run.json").write_text(json.dumps(record, indent=2)+"\n")
    return name, record["status"]


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("names", nargs="+")
    parser.add_argument("--workers", type=int, choices=range(1,9), default=2)
    args=parser.parse_args()
    with ThreadPoolExecutor(args.workers) as workers:
        for result in workers.map(retry, args.names):
            print(*result, flush=True)
