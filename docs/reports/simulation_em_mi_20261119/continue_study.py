"""Fill eight dataset slots using only run durations for the budget decision."""
import concurrent.futures as cf
import datetime as dt
import json
from pathlib import Path
import time

from simulate import ROOT, run_dataset


def main():
    runs = ROOT / "runs"
    initial = sorted(runs.iterdir())
    if len(initial) != 8:
        raise SystemExit("Start after launching the eight pilot-wave datasets in a fresh study.")
    start = min(dt.datetime.fromisoformat(json.loads(p.read_text())["started_utc"]).timestamp()
                for p in runs.glob("*/initial_run.json"))
    # Reserve 12 min for data audit, final summaries and wall-time variation.
    deadline = start + (4 * 60 - 12) * 60
    planned = {"alt": 5, "null": 3}
    durations = []
    # First target 12 mixed and 8 null datasets. Then balance the two scenarios.
    def next_job():
        if planned["alt"] < 12 and planned["null"] >= 8:
            arm = "alt"
        elif planned["alt"] < 12 and planned["alt"] / 12 < planned["null"] / 8:
            arm = "alt"
        else:
            arm = min(planned, key=planned.get)
        planned[arm] += 1
        return arm, planned[arm], False

    with cf.ThreadPoolExecutor(8) as pool:
        active = set()
        pending_pilots = set(initial)
        while True:
            for d in list(pending_pilots):
                p = d / "run.json"
                if not p.exists():
                    continue
                rec = json.loads(p.read_text())
                if rec["status"] != "ok":
                    raise RuntimeError(f"pilot {d.name} failed, inspect its log")
                durations.append(rec["elapsed_s"] + rec["generation_s"])
                pending_pilots.remove(d)
            for f in list(active):
                if f.done():
                    name, status, elapsed = f.result()
                    print(name, status, elapsed, flush=True)
                    if status != "ok":
                        raise RuntimeError(f"{name} failed, inspect its log")
                    durations.append(elapsed)
                    active.remove(f)
            if durations:
                reserve = 1.2 * max(durations[-24:])
                free = 8 - len(pending_pilots) - len(active)
                while free > 0 and time.time() + reserve < deadline:
                    job = next_job()
                    print("start", *job[:2], flush=True)
                    active.add(pool.submit(run_dataset, job))
                    free -= 1
                if not pending_pilots and not active:
                    break
            time.sleep(5)
    print("planned", planned, flush=True)


if __name__ == "__main__":
    main()
