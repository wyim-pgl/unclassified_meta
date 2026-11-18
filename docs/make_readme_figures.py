#!/usr/bin/env python3
"""Rebuild the README figures from a finished run: stacked bars + all-sample Sankey PNG.

Usage: python docs/make_readme_figures.py <run_output_dir> <name_prefix>
Needs `kaleido==0.2.1` for the Sankey PNG.
"""
import shutil
import sys
from pathlib import Path

import pandas as pd

from unclassified_meta import plots

run, prefix = Path(sys.argv[1]), sys.argv[2]
out = Path(__file__).resolve().parent / "img"
out.mkdir(exist_ok=True)
tab = run / "tables"
inputs = pd.read_csv(tab / "input_summary.tsv", sep="\t")
table = pd.read_csv(tab / "read_category_by_sample.tsv", sep="\t")
refs = pd.read_csv(tab / "reference_by_sample.tsv", sep="\t").fillna("")
shutil.copy(run / "figures" / "read_categories.png", out / f"{prefix}_read_categories.png")
html = out / f"{prefix}_sankey.html"
plots.sankey(inputs, table, refs, html)
html.unlink()
print("written:", sorted(p.name for p in out.glob(f"{prefix}_*")))
