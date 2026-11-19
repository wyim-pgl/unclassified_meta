"""taxa-table --merge-panel with expected counts, imputations and --panel-taxonomy."""
import pandas as pd
import pytest

from unclassified_meta import manifest, panel, taxa

from test_em_taxa import em_file
from test_ported_units import taxonomy_files


def quant(tmp_path, counts, in_frag, imputations=()):
    """A sealed panel-quant folder counted on the unclassified reads of sample s."""
    q = tmp_path / "quant"
    q.mkdir()
    ids = ["ns|V1|like", "ns|V2|same", "ns|V3|like"]
    pd.DataFrame({"feature_id": ids, "s": counts}).to_csv(q / "counts.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ids, "label": ["V1-like", "V2", "V3-like"],
                  "identity_class": ["like", "same", "like"]}).to_csv(q / "features.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": ["s"], "input_fragments": [in_frag]}).to_csv(q / "ledger.tsv", sep="\t", index=False)
    extra = []
    if imputations:
        (q / "imputations").mkdir()
    for m, c in enumerate(imputations):
        name = f"imputations/counts_{m + 1:03d}.tsv"
        pd.DataFrame({"feature_id": ids, "s": c}).to_csv(q / name, sep="\t", index=False)
        led = f"imputations/ledger_{m + 1:03d}.tsv"
        pd.DataFrame({"sample": ["s"], "imputation": [m + 1], "input_fragments": [in_frag],
                      "assigned_to_features": [sum(c)]}).to_csv(q / led, sep="\t", index=False)
        extra += [name, led]
    panel.write_quant_provenance(q, {"reads_scope": "unclassified", "imputations": len(imputations)}, extra)
    return q / "counts.tsv"


def test_expected_counts_reconcile_and_imputations_are_merged(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=3, seed=5)
    pc = quant(tmp_path, [1, 0, 0], 1, imputations=[[1, 0, 0], [0, 0, 1], [1, 0, 0]])
    taxa.merge_with_panel(out, pc, None, out / "merged")
    rec = pd.read_csv(out / "merged/reconciliation.tsv", sep="\t").iloc[0]
    assert rec.status == "ok" and abs(rec.difference) < 1e-6 and rec.total_input_reads == 27
    seal = manifest.verify(out / "merged")
    assert seal["imputations"] == 3
    merged = pd.read_csv(out / "merged/counts.tsv", sep="\t").set_index("feature_id")
    for m in range(1, 4):
        mc = pd.read_csv(out / f"merged/imputations/counts_{m:03d}.tsv", sep="\t").set_index("feature_id")
        assert list(mc.index) == list(merged.index) and mc["s"].dtype.kind == "i"
        tl = pd.read_csv(out / f"imputations/ledger_{m:03d}.tsv", sep="\t").iloc[0]
        assert mc["s"].sum() == tl.kept_at_rank + 1               # taxa kept + the one panel fragment


def test_unequal_imputation_numbers_are_refused(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=2)
    pc = quant(tmp_path, [1, 0, 0], 1, imputations=[[1, 0, 0], [1, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError, match="2 imputations and the panel table 3"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_panel_features_join_the_row_of_their_taxon(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=2)
    pc = quant(tmp_path, [1, 0, 0], 1, imputations=[[1, 0, 0], [0, 1, 0]])
    tx = tmp_path / "panel_taxonomy.tsv"
    pd.DataFrame({"feature_id": ["ns|V1|like", "ns|V2|same", "ns|V3|like"],
                  "taxid": ["83333", "561", "28901"]}).to_csv(tx, sep="\t", index=False)
    taxa.merge_with_panel(out, pc, None, out / "merged", tx, nodes, names, ("Viruses",))
    c = pd.read_csv(out / "merged/counts.tsv", sep="\t").set_index("feature_id")["s"]
    assert "ns|V1|like" not in c and c["t|562|S"] == pytest.approx(19, abs=1e-4)   # 18 classified + 1 panel (strain -> species)
    assert "ns|V2|same" in c                                         # genus: above the species rank, stays
    assert "ns|V3|like" not in c and c["t|28901|S"] == 0             # new taxon row (Salmonella had ~0 classified)
    plan = pd.read_csv(out / "merged/collapsed.tsv", sep="\t").set_index("feature_id")
    assert plan.loc["ns|V1|like", "status"] == "collapsed" and plan.loc["ns|V2|same", "status"] == "above_rank"
    f = pd.read_csv(out / "merged/features.tsv", sep="\t", keep_default_na=False).set_index("feature_id")
    assert f.loc["t|562|S", "collapsed_panel_features"] == "ns|V1|like"
    assert f.loc["t|28901|S", "label"] == "Salmonella enterica" and list(f.index) == list(c.index)
    for m, want in ((1, 1), (2, 0)):
        mc = pd.read_csv(out / f"merged/imputations/counts_{m:03d}.tsv", sep="\t").set_index("feature_id")["s"]
        assert list(mc.index) == list(c.index)
        assert mc["ns|V2|same"] == (0 if m == 1 else 1)
    rec = pd.read_csv(out / "merged/reconciliation.tsv", sep="\t").iloc[0]
    assert rec.status == "ok"
    assert manifest.verify(out / "merged")["panel_taxonomy_sha256"]
