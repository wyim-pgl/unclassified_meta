"""One self-contained HTML summary page (tables + bar chart, links to the Sankey pages)."""
from __future__ import annotations

import base64
import html
from pathlib import Path

import pandas as pd

from . import __version__


def _table(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df is None or df.empty:
        return "<p><i>none</i></p>"
    return df.head(max_rows).to_html(index=False, border=0, classes="t", float_format=lambda v: f"{v:.4g}")


def write_report(outdir: Path, cfg: dict, inputs: pd.DataFrame, table: pd.DataFrame, refs: pd.DataFrame,
                 contigs: pd.DataFrame, bar_png: Path, sankeys: list[Path], resources: pd.DataFrame) -> Path:
    img = base64.b64encode(bar_png.read_bytes()).decode() if bar_png.exists() else ""
    frac = table[["sample", "reads_used"] + [c for c in table.columns if c.startswith("frac_")]].copy()
    for c in frac.columns:
        if c.startswith("frac_"):
            frac[c] = (100 * frac[c]).round(2)
    frac.columns = [c.replace("frac_", "% ") for c in frac.columns]
    cat_counts = contigs.groupby("category").agg(contigs=("contig", "size"), total_bp=("length", "sum")) \
        .reset_index().sort_values("contigs", ascending=False)
    links = "".join(f'<li><a href="{html.escape(str(p.relative_to(outdir)))}">{html.escape(p.stem)}</a></li>'
                    for p in sankeys)
    css = ("body{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}"
           "table.t{border-collapse:collapse;font-size:12px;margin:8px 0 20px}"
           "table.t th,table.t td{border-bottom:1px solid #ddd;padding:3px 8px;text-align:right}"
           "table.t td:first-child,table.t th:first-child{text-align:left}h2{margin-top:32px}"
           "img{max-width:100%}")
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>unclassified-meta report</title>
<style>{css}</style></head><body>
<h1>unclassified-meta report</h1>
<p>Mode <b>{cfg['mode']}</b>, assembler <b>{html.escape(str(cfg.get('_assembler', '')))}</b>,
subsample {cfg['subsample'] if cfg['subsample'] else 'all reads'} per sample, contigs &ge; {cfg['min_contig_len']} bp,
layer order {', '.join(cfg['priority'])}. unclassified-meta {__version__}.</p>
<h2>1. Input reads</h2>{_table(inputs)}
<h2>2. What the unclassified reads are</h2>
<img src="data:image/png;base64,{img}" alt="stacked bar chart">
<p>Sankey diagrams (interactive):</p><ul>{links}</ul>
<h2>3. Read categories per sample (% of reads used)</h2>{_table(frac)}
<h2>4. Top reference hits (pooled)</h2>
{_table(refs.sort_values('reads_pooled', ascending=False)[['category', 'best_hit_title', 'n_contigs', 'pct_of_reads_used_pooled', 'blastn_identity', 'tblastx_identity', 'diamond_identity']] if not refs.empty else refs)}
<h2>5. Contigs per category</h2>{_table(cat_counts)}
<h2>6. Resources of the heavy steps</h2>{_table(resources)}
<p>Full tables are in <code>tables/</code>; figures in <code>figures/</code>; the log is <code>logs/pipeline.log</code>.</p>
</body></html>"""
    out = outdir / "report.html"
    out.write_text(page)
    return out
