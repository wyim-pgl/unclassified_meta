"""Output compatibility of two taxa-table runs (baseline, optimised): compat.py <base out> <opt out>
Checks what downstream steps (merge, da --imputations) rely on, for the table folder and merged/ if present."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from unclassified_meta import manifest
from unclassified_meta.taxa import read_draw


def folder(p: Path):
    name = "provenance.json" if (p / "provenance.json").exists() else manifest.NAME
    rec = manifest.verify(p, name, what="count tables")          # raises when a sealed file does not match
    return rec


def check(a: Path, b: Path):
    ra, rb = folder(a), folder(b)
    print(f"[{a.name} vs {b.name}] manifests verify; files identical set: {set(ra['artifacts']) == set(rb['artifacts'])}"
          f"; imputations {ra.get('imputations')} / {rb.get('imputations')}; seed {ra.get('seed')} / {rb.get('seed')}")
    ca, cb = pd.read_csv(a / "counts.tsv", sep="\t"), pd.read_csv(b / "counts.tsv", sep="\t")
    print(f"  columns equal {list(ca.columns) == list(cb.columns)}; features {len(ca)} / {len(cb)}, "
          f"only in base {len(set(ca.feature_id) - set(cb.feature_id))}, only in opt "
          f"{len(set(cb.feature_id) - set(ca.feature_id))}")
    samples = [c for c in ca.columns if c != "feature_id" and pd.api.types.is_numeric_dtype(ca[c])]
    m = ca.set_index("feature_id")[samples].join(cb.set_index("feature_id")[samples], rsuffix="_o", how="outer").fillna(0)
    d = np.abs(m[samples].values - m[[s + "_o" for s in samples]].values)
    print(f"  expected counts: column sums equal to 1e-5: "
          f"{np.allclose(m[samples].sum().values, m[[s + '_o' for s in samples]].sum().values, atol=1e-5)}; "
          f"max |diff| {d.max():.3g} reads; cells |diff| > 0.5: {(d > 0.5).sum()} of {d.size}")
    n = int(rb.get("imputations") or 0)
    rows = list(cb["feature_id"])
    for k in range(n):
        f = b / f"imputations/counts_{k + 1:03d}.tsv"
        t = read_draw(f)                          # whole reads only (raises otherwise)
        fa = read_draw(a / f"imputations/counts_{k + 1:03d}.tsv")
        assert list(t.index) == rows, f"{f}: rows differ from counts.tsv (da --imputations would refuse it)"
        assert list(t.columns) == list(fa.columns), f"{f}: sample columns differ"
        if (b / f"imputations/ledger_{k + 1:03d}.tsv").exists():
            for side in (a, b):
                led = pd.read_csv(side / f"imputations/ledger_{k + 1:03d}.tsv", sep="\t").set_index("sample")
                parts = led[["unclassified", "below_threshold", "dropped_lineages", "kept_at_rank",
                             "unresolved_above_rank"]].sum(axis=1)
                assert (parts == led["input_fragments"]).all(), f"{side}: imputation {k + 1} ledger does not balance"
                assert (led["kept_at_rank"] == read_draw(side / f"imputations/counts_{k + 1:03d}.tsv")
                        .sum().reindex(led.index)).all(), f"{side}: imputation {k + 1} kept reads != table"
    print(f"  {n} imputation tables: whole reads, rows = counts.tsv rows, same samples"
          + (", every imputation ledger balances (input = unclassified + below + dropped + kept + unresolved) and "
             "matches its table" if (b / "imputations/ledger_001.tsv").exists() else ""))
    only = set(ca.feature_id) ^ set(cb.feature_id)
    both = pd.concat([ca, cb]).set_index("feature_id")[samples]
    print(f"  features in one run only: {len(only)}; their largest expected count in any sample "
          f"{both.loc[both.index.isin(only)].max().max():.3g} reads (rows kept because an imputation drew them)")
    if (a / "ledger.tsv").exists():                  # merged/ has no ledger of its own
        la, lb = pd.read_csv(a / "ledger.tsv", sep="\t"), pd.read_csv(b / "ledger.tsv", sep="\t")
        print(f"  ledger columns equal {list(la.columns) == list(lb.columns)}")


A, B = Path(sys.argv[1]), Path(sys.argv[2])
check(A, B)
if (A / "merged").exists():
    check(A / "merged", B / "merged")
