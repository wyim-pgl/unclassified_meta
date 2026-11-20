"""Generate truth and hit files, then invoke the production pipeline on each dataset."""
import argparse
import collections
import concurrent.futures
import datetime
import itertools
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DEPTH = 200_000
PROBS = [0.30, 0.40, 0.60, 0.70, 0.50, 0.50, 0, 0] * 5
UP, DOWN = {0, 10, 20}, {5, 25, 35}
SAMPLES = [f"s{s:02d}" for s in range(12)]
GROUPS = ["control"] * 6 + ["treated"] * 6
HEADER = "readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n"


def taxonomy():
    rows = [(1, 1, "no rank", "root"), (2, 1, "superkingdom", "Synthetic domain")]
    for g in range(5):
        rows.append((10 + g, 2, "genus", f"Genus{g}"))
        rows.extend((100 + g * 8 + j, 10 + g, "species", f"Genus{g} species{j}") for j in range(8))
    (ROOT / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n" for t, p, r, _ in rows))
    (ROOT / "names.dmp").write_text("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in rows))


def candidates(i):
    """Conditional on ambiguity: uniform number of sisters, then uniform subset."""
    j = i % 8
    members = range(i - j, i - j + 4) if j < 4 else range(i - j + 4, i - j + 6)
    sisters = [s for s in members if s != i]
    sets, weights = [], []
    for k in range(1, len(sisters) + 1):
        for others in itertools.combinations(sisters, k):
            sets.append(tuple(sorted((i, *others))))
            weights.append(1 / (len(sisters) * math.comb(len(sisters), k)))
    return sets, weights


def generate(folder, scenario, seed):
    rng = np.random.default_rng(seed)
    mu = rng.normal(0, 0.85, 40)
    beta = np.array([math.log(4) * ((i in UP) - (i in DOWN)) for i in range(40)])
    if scenario == "null":
        beta[:] = 0
    ids = [f"sim|{100+i}|S" for i in range(40)]
    truth = pd.DataFrame(dict(feature_id=ids, baseline_log=mu, true_lfc=beta, prob_ambig=PROBS))
    truth.to_csv(folder / "truth.tsv", sep="\t", index=False)
    pd.DataFrame(dict(sample=SAMPLES, group=GROUPS)).to_csv(folder / "meta.tsv", sep="\t", index=False)
    counts, singles = {}, {}
    for sample, group in zip(SAMPLES, GROUPS):
        abundance = np.exp(mu + (group == "treated") * beta + rng.normal(0, 0.65, 40))
        counts[sample] = rng.multinomial(DEPTH, abundance / abundance.sum())
        amb = rng.binomial(counts[sample], PROBS)
        singles[sample] = counts[sample] - amb
        classes = collections.Counter()
        for i, n in enumerate(amb):
            classes[(i,)] += singles[sample][i]
            if n:
                sets, weights = candidates(i)
                classes.update(dict(zip(sets, rng.multinomial(n, weights))))
        index = 0
        with (folder / f"{sample}.cfr.tsv").open("w") as f:
            f.write(HEADER)
            for hits, n in sorted(classes.items()):
                pattern = "".join(f"r{{}}\tref{100+j}\t{100+j}\t900\t0\t100\t100\t{len(hits)}\n" for j in hits)
                for rid in range(index, index + n):
                    f.write(pattern.format(*([rid] * len(hits))))
                index += n
        assert index == DEPTH
    for name, table in [("origin", counts), ("singleton", singles)]:
        pd.DataFrame(table, index=ids).rename_axis("feature_id").to_csv(folder / f"{name}.tsv", sep="\t")


def run_dataset(job):
    scenario, rep, reuse = job
    seed = 202611190 + 2 * rep + int(scenario == "alt")
    folder = ROOT / "runs" / f"{scenario}_{rep:03d}"
    folder.mkdir(parents=True, exist_ok=reuse)
    started = time.monotonic()
    record = {"scenario": scenario, "rep": rep, "seed": seed, "pythonpath": os.environ["PYTHONPATH"], "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "commands": []}
    def run(label, args):
        cmd = [sys.executable, "-m", "unclassified_meta.cli", *map(str, args)]
        t = time.monotonic()
        with (folder / f"{label}.log").open("w") as f:
            p = subprocess.run(cmd, cwd=folder, stdout=f, stderr=subprocess.STDOUT)
        record["commands"].append(dict(label=label, command=shlex.join(cmd), seconds=time.monotonic()-t, exit=p.returncode))
        p.check_returncode()

    try:
        if reuse:
            old = json.loads((folder / "run.json").read_text())
            assert old["status"] == "failed" and old["seed"] == seed
            (folder / "run.json").rename(folder / "initial_run.json")
            (folder / "table_lca.log").rename(folder / "initial_table_lca.log")
            record["generation_s"] = old["generation_s"]
        else:
            generate(folder, scenario, seed)
            record["generation_s"] = time.monotonic() - started
        pd.DataFrame({"sample": SAMPLES, "report": [s + ".cfr.tsv" for s in SAMPLES]}).to_csv(
            folder / "reports.tsv", sep="\t", index=False)
        common = ["taxa-table", "--reports", "reports.tsv", "--per-read", "--nodes", "../../nodes.dmp",
                  "--names", "../../names.dmp", "--taxonomy-name", "sim", "--rank", "S"]
        run("table_lca", [*common, "--assign", "lca", "-o", "lca"])
        run("table_em", [*common, "--assign", "em", "--imputations", 5, "--seed", seed, "-o", "em"])
        for method in ["lca", "em"]:
            args = ["da", "--counts", f"{method}/counts.tsv", "--features", f"{method}/features.tsv",
                    "--metadata", "meta.tsv", "--group", "group", "--formula", "group", "--reference-level", "control",
                    "--p-adjust", "BH", "--alpha", "0.05", "--threads", 1, "-o", f"da_{method}"]
            run("da_" + method, args + (["--imputations"] if method == "em" else []))
        record["status"] = "ok"
    except Exception as e:
        record.update(status="failed", error=repr(e))
    finally:
        record["elapsed_s"] = time.monotonic() - started
        (folder / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    return folder.name, record["status"], round(record["elapsed_s"], 1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--alt", nargs=2, type=int, metavar=("FIRST", "LAST"), default=[1, 0])
    p.add_argument("--null", nargs=2, type=int, metavar=("FIRST", "LAST"), default=[1, 0])
    p.add_argument("--workers", type=int, choices=range(1, 9), default=8)
    p.add_argument("--reuse-inputs", action="store_true", help="retry original input-directory overlap failures")
    opt = p.parse_args()
    if not (ROOT / "nodes.dmp").exists():
        taxonomy()
    jobs = sorted(((scenario, r) for scenario in ["alt", "null"] for r in range(getattr(opt, scenario)[0], getattr(opt, scenario)[1]+1)), key=lambda j: j[1])
    with concurrent.futures.ThreadPoolExecutor(opt.workers) as ex:
        futures = [ex.submit(run_dataset, (*job, opt.reuse_inputs)) for job in jobs]
        failed = False
        for f in concurrent.futures.as_completed(futures):
            result = f.result()
            print(*result, flush=True)
            failed |= result[1] != "ok"
    raise SystemExit(1 if failed else 0)
