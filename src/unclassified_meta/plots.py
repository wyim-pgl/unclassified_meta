"""Figures: stacked bars per sample (matplotlib) and Sankey diagrams (plotly HTML)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

PALETTE = {"virus": "#E64B35", "bacteria": "#4DBBD5", "archaea": "#3C5488", "fungi": "#F39B7F",
           "human": "#8491B4", "animal": "#91D1C2", "plant": "#00A087", "other_eukaryote": "#B09C85",
           "protein": "#7E6148", "rRNA": "#DC0000", "organelle": "#9FD39F", "other": "#C9B26B",
           "kraken2_hit": "#6A3D9A", "unclassified": "#B3B3B3", "not_assembled": "#E0E0E0",
           "classified_by_classifier": "#404040"}
HOST_COLORS = ["#2E7D32", "#66BB6A", "#A5D6A7", "#1B5E20"]


def color(cat: str, hosts: list[str]) -> str:
    if cat.startswith("host:"):
        return HOST_COLORS[hosts.index(cat) % len(HOST_COLORS)]
    return PALETTE.get(cat, "#999999")


def label(cat: str) -> str:
    if cat.startswith("host:"):
        return f"host genome ({cat[5:]})"
    return {"not_assembled": "not assembled", "protein": "plant protein match", "rRNA": "rRNA",
            "classified_by_classifier": "classified by the read classifier"}.get(cat, cat.replace("_", " "))


def categories(table: pd.DataFrame) -> list[str]:
    return [c[5:] for c in table.columns if c.startswith("frac_")]


def stacked_bars(table: pd.DataFrame, inputs: pd.DataFrame, out_prefix: Path) -> list[Path]:
    """Panel A: categories as % of unclassified reads; panel B (when classifier totals exist): % of all reads."""
    cats = categories(table)
    hosts = [c for c in cats if c.startswith("host:")]
    has_total = bool((inputs["classified_reads"] > 0).any())
    n = len(table)
    fig, axes = plt.subplots(1 + has_total, 1, figsize=(max(6, 0.45 * n + 3), 4.2 * (1 + has_total)),
                             squeeze=False)
    x = range(n)
    ax = axes[0][0]
    bottom = [0.0] * n
    for c in cats:
        v = (100 * table[f"frac_{c}"]).tolist()
        ax.bar(x, v, bottom=bottom, color=color(c, hosts), label=label(c), width=0.8)
        bottom = [b + y for b, y in zip(bottom, v)]
    ax.set_ylabel("% of unclassified reads")
    ax.set_ylim(0, 100)
    ax.set_xticks(list(x))
    ax.set_xticklabels(table["sample"], rotation=60, ha="right", fontsize=8)
    ax.set_title("A  What the unclassified reads are", loc="left", fontsize=10)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    if has_total:
        ax = axes[1][0]
        tot = inputs.set_index("sample").loc[table["sample"], "total_reads"].to_numpy()
        clf = inputs.set_index("sample").loc[table["sample"], "classified_reads"].to_numpy()
        v = 100 * clf / tot
        ax.bar(x, v, color=PALETTE["classified_by_classifier"], label=label("classified_by_classifier"), width=0.8)
        bottom = list(v)
        for c in cats:
            w = (100 * table[f"est_{c}"].to_numpy() / tot).tolist()
            ax.bar(x, w, bottom=bottom, color=color(c, hosts), label=label(c), width=0.8)
            bottom = [b + y for b, y in zip(bottom, w)]
        ax.set_ylabel("% of all reads")
        ax.set_ylim(0, 100)
        ax.set_xticks(list(x))
        ax.set_xticklabels(table["sample"], rotation=60, ha="right", fontsize=8)
        ax.set_title("B  All reads, including those the classifier assigned", loc="left", fontsize=10)
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    for a in axes.flat:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    paths = []
    for ext in ("png", "pdf"):
        p = out_prefix.with_suffix(f".{ext}")
        fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def sankey(inputs: pd.DataFrame, table: pd.DataFrame, refs: pd.DataFrame, out_html: Path,
           sample: str | None = None, top_refs: int = 8, min_ref_share: float = 0.01) -> Path:
    """Total reads -> classified / unclassified -> categories -> (virus) top reference genomes."""
    import plotly.graph_objects as go

    cats = categories(table)
    hosts = [c for c in cats if c.startswith("host:")]
    t = table if sample is None else table[table["sample"] == sample]
    i = inputs if sample is None else inputs[inputs["sample"] == sample]
    total, classified, unclassified = int(i["total_reads"].sum()), int(i["classified_reads"].sum()), \
        int(i["unclassified_reads"].sum())
    nodes, colors, src, dst, val, lcol = [], [], [], [], [], []

    def node(name, col):
        nodes.append(name)
        colors.append(col)
        return len(nodes) - 1

    def link(a, b, v, col):
        if v > 0:
            src.append(a)
            dst.append(b)
            val.append(v)
            lcol.append(col)

    start = node(f"All reads<br>{total:,}", "#404040")
    unc = node(f"Unclassified<br>{unclassified:,}", "#B3B3B3")
    if classified:
        link(start, node(f"Classified by read classifier<br>{classified:,}", PALETTE["classified_by_classifier"]),
             classified, "rgba(64,64,64,0.3)")
    link(start, unc, unclassified, "rgba(150,150,150,0.35)")
    for c in cats:
        v = int(t[f"est_{c}"].sum())
        if v <= 0:
            continue
        col = color(c, hosts)
        n = node(f"{label(c)}<br>{v:,} ({100 * v / max(unclassified, 1):.1f} %)", col)
        link(unc, n, v, _rgba(col, 0.45))
        if c == "virus" and refs is not None and not refs.empty:
            vr = refs[refs["category"] == "virus"].copy()
            if sample is not None:
                vr["w"] = vr.get(f"pct_{sample}", 0)
            else:
                vr["w"] = vr["reads_pooled"]
            vr = vr.sort_values("w", ascending=False)
            share_total = vr["w"].sum()
            if share_total > 0:
                shown = 0
                vr = vr[vr["w"] >= min_ref_share * share_total]  # tiny hits go to "other viruses"
                for _, r in vr.head(top_refs).iterrows():
                    part = int(round(v * r["w"] / share_total))
                    name = (r["best_hit_title"] or r["best_hit"])[:60]
                    link(n, node(name, col), part, _rgba(col, 0.3))
                    shown += part
                if v - shown > 0:
                    link(n, node("other viruses", col), v - shown, _rgba(col, 0.3))
    fig = go.Figure(go.Sankey(arrangement="snap",
                              node=dict(label=nodes, color=colors, pad=14, thickness=16),
                              link=dict(source=src, target=dst, value=val, color=lcol)))
    fig.update_layout(title=f"Read fate{'' if sample is None else ' — ' + sample} "
                            "(counts after subsampling are scaled to all unclassified reads)",
                      font_size=11, height=650)
    fig.write_html(out_html, include_plotlyjs="cdn")
    _static_png(fig, out_html.with_suffix(".png"))
    return out_html


def _static_png(fig, out_png: Path) -> None:
    """Also save a PNG of a plotly figure when the optional `kaleido` package is installed."""
    try:
        import kaleido  # noqa: F401
    except ImportError:
        print(f"[unclassified-meta] note: kaleido not installed, {out_png.name} not written "
              "(pip install kaleido==0.2.1 to get static Sankey PNGs)")
        return
    fig.write_image(out_png, width=1100, height=650, scale=2)


def _rgba(hexcol: str, alpha: float) -> str:
    h = hexcol.lstrip("#")
    r, g, b = (int(h[k:k + 2], 16) for k in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
