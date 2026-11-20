"""Independent arithmetic check of the report's primary and extended endpoints."""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

root = Path(__file__).resolve().parent

def rows(name):
    with (root/name).open() as f:
        return list(csv.DictReader(f,delimiter="\t"))

data = rows("results.tsv")
runtime = rows("runtimes.tsv")
accepted = {r["dataset"] for r in runtime if r["status"] == "ok"}
assert accepted == {r["dataset"] for r in data}
assert len({(r["dataset"],r["method"],r["feature_id"]) for r in data}) == len(data) == 120*len(accepted)
proofs = {r["path"] for r in rows("repool_audit.tsv")}
for row in runtime:
    if row["status"] != "ok":
        continue
    path = "runs/"+row["dataset"]+("/pinned_fit" if row["source_cohort"]=="pinned_replay" else "")
    assert path in proofs
for strict, filename in [(True,"strict_metrics.tsv"),(False,"metrics.tsv")]:
    summary = {(r["scenario"],r["method"],r["subset"]):r for r in rows(filename)}
    for method in ["LCA","EM","MI"]:
        z = [r for r in data if r["method"]==method and (not strict or r["strict_revision"]=="True")]
        alt = [r for r in z if r["scenario"]=="alt" and float(r["true_lfc"])!=0]
        null = [r for r in z if r["scenario"]=="global_null"]
        pooled = summary["alt",method,"changed"]
        calls = sum(float(r["q"])<.05 for r in alt)
        raw = sum(float(r["p"])<.05 for r in null)
        cover = sum(abs(float(r["lfc_centered"])-float(r["true_lfc"])) <= float(r["halfwidth"]) for r in alt)
        assert all(math.isclose(x,float(pooled[k]),abs_tol=1e-12) for k,x in [
            ("bh",calls/len(alt)), ("coverage",cover/len(alt)),
            ("bias",mean(float(r["error"]) for r in alt)), ("rmse",math.sqrt(mean(float(r["error"])**2 for r in alt)))])
        assert math.isclose(raw/len(null),float(summary["global_null",method,"all"]["p_lt05"]),abs_tol=1e-12)
        print("strict" if strict else "all",method,"TP",calls,len(alt),"covered",cover,len(alt),"raw null positives",raw,len(null))
print(len(runtime),"seeds attempted;",len(accepted),"accepted;",len(data),"unique species-method rows")
