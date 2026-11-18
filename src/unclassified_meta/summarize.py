"""Tables: input summary, per-sample read categories, per-reference read shares."""
from __future__ import annotations

from collections import Counter, defaultdict

import pandas as pd

CATEGORY_ORDER = ["virus", "bacteria", "archaea", "fungi", "human", "animal", "plant", "other_eukaryote",
                  "protein", "rRNA", "organelle", "other", "kraken2_hit", "unclassified", "not_assembled"]


def ordered_categories(cats) -> list[str]:
    cats = set(cats)
    hosts = sorted(c for c in cats if c.startswith("host:"))
    rest = [c for c in CATEGORY_ORDER if c in cats]
    extra = sorted(cats - set(hosts) - set(rest))
    # host first (expected), then informative classes, then leftovers
    return hosts + [c for c in rest if c not in ("unclassified", "not_assembled")] + extra + \
        [c for c in ("unclassified", "not_assembled") if c in cats]


def read_category_table(inputs: pd.DataFrame, per_sample: dict[str, dict[str, int]]) -> pd.DataFrame:
    """One row per sample: counts n_<cat>, fractions frac_<cat> (of reads used) and est_<cat> (scaled to all
    unclassified reads of the sample)."""
    cats = ordered_categories({c for d in per_sample.values() for c in d} | {"not_assembled"})
    rows = []
    for _, r in inputs.iterrows():
        counts = per_sample.get(r["sample"], {})
        used = int(r["reads_used"])
        row = {"sample": r["sample"], "reads_used": used, "unclassified_reads": int(r["unclassified_reads"])}
        for c in cats:
            n = int(counts.get(c, 0))
            row[f"n_{c}"] = n
            row[f"frac_{c}"] = n / used if used else 0.0
            row[f"est_{c}"] = round(row[f"frac_{c}"] * int(r["unclassified_reads"]))
        rows.append(row)
    return pd.DataFrame(rows)


def reference_table(contigs: pd.DataFrame, contig_counts: dict[str, Counter], reads_used: dict[str, int]) -> pd.DataFrame:
    """Read share per (category, best hit), pooled and per sample, with read-weighted mean identity."""
    info = contigs.set_index("contig")
    agg: dict[tuple, dict] = defaultdict(lambda: {"n_contigs": 0, "reads": Counter(), "id_w": 0.0, "id_n": 0})
    for c, row in info.iterrows():
        key = (row["category"], row["best_hit"] if isinstance(row["best_hit"], str) and row["best_hit"] else "-",
               row["best_hit_title"] if isinstance(row["best_hit_title"], str) else "")
        agg[key]["n_contigs"] += 1
        for s, cnt in contig_counts.items():
            n = cnt.get(c, 0)
            if n:
                agg[key]["reads"][s] += n
                if pd.notna(row.get("identity")) and row.get("identity") != "":
                    agg[key]["id_w"] += float(row["identity"]) * n
                    agg[key]["id_n"] += n
    total_used = sum(reads_used.values()) or 1
    rows = []
    for (cat, hit, title), a in agg.items():
        r = {"category": cat, "best_hit": hit, "best_hit_title": title, "n_contigs": a["n_contigs"],
             "reads_pooled": sum(a["reads"].values()),
             "pct_of_reads_used_pooled": round(100 * sum(a["reads"].values()) / total_used, 3),
             "read_weighted_identity": round(a["id_w"] / a["id_n"], 2) if a["id_n"] else None}
        for s, used in reads_used.items():
            r[f"pct_{s}"] = round(100 * a["reads"][s] / used, 3) if used else 0.0
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["category", "reads_pooled"], ascending=[True, False]).reset_index(drop=True)
