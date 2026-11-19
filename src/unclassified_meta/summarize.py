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


IDENTITY_METHODS = ("blastn", "tblastx", "diamond", "minimap2")


def identity_method(method) -> str | None:
    """Which identity a classification method reports: BLASTN and minimap2 nucleotide identity, TBLASTX and
    DIAMOND amino-acid identity. They are not comparable, so they are kept in separate columns."""
    m = str(method or "")
    if "tblastx" in m:
        return "tblastx"
    if "diamond" in m:
        return "diamond"
    if "blastn" in m:
        return "blastn"
    if "minimap2" in m:
        return "minimap2"
    return None


def reference_table(contigs: pd.DataFrame, contig_counts: dict[str, Counter], reads_used: dict[str, int]) -> pd.DataFrame:
    """Read share per (category, best hit), pooled and per sample, with the read-weighted identity per method
    (blastn_identity, tblastx_identity, diamond_identity, minimap2_identity)."""
    from .config import REFERENCE_TABLE_COLUMNS
    clash = sorted(s for s in reads_used if f"pct_{s}" in REFERENCE_TABLE_COLUMNS)
    if clash:
        raise ValueError(f"sample name(s) {', '.join(clash)} would overwrite a column of the reference table")
    info = contigs.set_index("contig")
    agg: dict[tuple, dict] = defaultdict(lambda: {"n_contigs": 0, "reads": Counter(),
                                                  "id": {m: [0.0, 0] for m in IDENTITY_METHODS}})
    for c, row in info.iterrows():
        key = (row["category"], row["best_hit"] if isinstance(row["best_hit"], str) and row["best_hit"] else "-",
               row["best_hit_title"] if isinstance(row["best_hit_title"], str) else "")
        agg[key]["n_contigs"] += 1
        for s, cnt in contig_counts.items():
            n = cnt.get(c, 0)
            if n:
                agg[key]["reads"][s] += n
                m = identity_method(row.get("method"))
                if m and pd.notna(row.get("identity")) and row.get("identity") != "":
                    agg[key]["id"][m][0] += float(row["identity"]) * n
                    agg[key]["id"][m][1] += n
    total_used = sum(reads_used.values()) or 1
    rows = []
    for (cat, hit, title), a in agg.items():
        r = {"category": cat, "best_hit": hit, "best_hit_title": title, "n_contigs": a["n_contigs"],
             "reads_pooled": sum(a["reads"].values()),
             "pct_of_reads_used_pooled": 100 * sum(a["reads"].values()) / total_used}
        for m in IDENTITY_METHODS:
            w, n = a["id"][m]
            r[f"{m}_identity"] = round(w / n, 2) if n else None
        for s, used in reads_used.items():
            r[f"pct_{s}"] = 100 * a["reads"][s] / used if used else 0.0   # full precision; rounded only for display
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["category", "reads_pooled"], ascending=[True, False]).reset_index(drop=True)


def contig_count_matrix(contigs: pd.DataFrame, contig_counts: dict[str, Counter]) -> pd.DataFrame:
    """Raw fragment counts, contig x sample (every contig of contig_category.tsv, zeros included). The counts
    come from the reads used for read-back (see `reads_used` in input_summary.tsv), not scaled or rounded."""
    df = pd.DataFrame({"contig": list(contigs["contig"])})
    for s, cnt in contig_counts.items():
        df[s] = [int(cnt.get(c, 0)) for c in df["contig"]]
    return df
