"""Findings of the second review of the EM / imputation change."""
import json

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import da, manifest, taxa

from test_em_da import result
from test_em_merge import quant
from test_em_taxa import em_file
from test_ported_units import taxonomy_files


def test_a_hash_inside_a_label_is_data(tmp_path):
    f = tmp_path / "r.tsv"
    f.write_text("# note\nfeature_id\tlabel\tlfc_groupB\nf1\tstrain #2\t1.5\n")
    r = da.read_results(f)
    assert r.loc[0, "label"] == "strain #2" and r.loc[0, "lfc_groupB"] == 1.5


def test_a_fit_whose_own_test_failed_does_not_count():
    fits = [result([2.0, 0.1], [0.5, 0.5]).assign(**{"W_groupB": [4.0, 0.2], "p_groupB": [1e-4, 0.8]}) for _ in range(3)]
    fits[1].loc[0, "p_groupB"] = np.nan
    out, _ = da.combine_imputations(fits, 0.05, "holm")
    assert out.loc[0, "n_fits_groupB"] == 2 and out.loc[0, "p_groupB"] == 1.0 and not bool(out.loc[0, "diff_groupB"])


def test_stale_cleanup_never_leaves_the_folder(tmp_path):
    out = tmp_path / "out"
    (out / "imputations").mkdir(parents=True)
    victim = tmp_path / "victim.tsv"
    victim.write_text("x\n")
    (out / "imputations" / "counts_001.tsv").write_text("y\n")
    (out / "imputations" / "counts_002.tsv").write_text("replaced by the user\n")
    manifest.remove_stale(out, {"imputations/../../victim.tsv": manifest.sha256(victim),
                                "imputations/counts_001.tsv": manifest.sha256(out / "imputations/counts_001.tsv"),
                                "imputations/counts_002.tsv": "0" * 64}, [])
    assert victim.exists() and not (out / "imputations" / "counts_001.tsv").exists()
    assert (out / "imputations" / "counts_002.tsv").exists()            # not the sealed file any more: kept


def test_imputation_rows_must_be_the_merged_rows(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", assign="em",
                                   imputations=2)
    pc = quant(tmp_path, [1, 0, 0], 1, imputations=[[1, 0, 0], [0, 1, 0]])
    bad = pc.parent / "imputations/counts_002.tsv"
    pd.DataFrame({"feature_id": ["ns|V1|like", "ns|V2|same", "ns|EXTRA|like"], "s": [0, 1, 0]}).to_csv(
        bad, sep="\t", index=False)
    rec = json.loads((pc.parent / "provenance.json").read_text())
    rec["artifacts"]["imputations/counts_002.tsv"] = manifest.sha256(bad)
    (pc.parent / "provenance.json").write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="imputation 2: its rows are not the rows"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_a_panel_taxonomy_of_another_panel_is_refused(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t")
    q = tmp_path / "quant"
    pc = quant(tmp_path, [1, 0, 0], 1)
    rec = json.loads((q / "provenance.json").read_text())
    rec["panel"] = {"panel_id": "aaaa"}
    (q / "provenance.json").write_text(json.dumps(rec))
    d = tmp_path / "ptx"
    d.mkdir()
    pd.DataFrame({"feature_id": ["ns|V1|like"], "taxid": ["562"]}).to_csv(d / "panel_taxonomy.tsv", sep="\t", index=False)
    manifest.begin(d, "panel-taxonomy")
    manifest.seal(d, ["panel_taxonomy.tsv"], {"step": "panel-taxonomy", "panel_id": "bbbb"})
    with pytest.raises(ValueError, match="made for another panel"):
        taxa.merge_with_panel(out, pc, None, out / "m", d / "panel_taxonomy.tsv", nodes, names, ("Viruses",))


def test_reads_that_may_be_the_panels_organism_leave_the_background_whole(tmp_path):
    from test_ported_units import CFR_HEADER, cfr
    nodes, names = taxonomy_files(tmp_path)
    f = tmp_path / "s.cfr.tsv"
    f.write_text(CFR_HEADER + cfr("tie", 70001, 400, 60) + cfr("tie", 562, 400, 60) + cfr("ec", 562, 400, 60))
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", (), "t", panel_lineage=("Viruses",),
                                   assign="em", imputations=2)
    c = pd.read_csv(out / "counts.tsv", sep="\t").set_index("feature_id")["s"]
    led = pd.read_csv(out / "ledger.tsv", sep="\t").iloc[0]
    assert c["t|562|S"] == 1 and led.dropped_lineages == 1         # not 1.5 E. coli: the tie leaves whole
    il = pd.read_csv(out / "imputations/ledger_001.tsv", sep="\t").iloc[0]
    assert il.dropped_lineages == 1 and il.kept_at_rank == 1


def test_a_merged_draw_must_match_its_own_ledger(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", assign="em",
                                   imputations=2)
    led = out / "imputations/ledger_002.tsv"
    df = pd.read_csv(led, sep="\t")
    df.loc[0, "kept_at_rank"] += 1
    df.to_csv(led, sep="\t", index=False)
    rec = json.loads((out / "MANIFEST.json").read_text())
    rec["artifacts"]["imputations/ledger_002.tsv"] = manifest.sha256(led)
    (out / "MANIFEST.json").write_text(json.dumps(rec))
    pc = quant(tmp_path, [1, 0, 0], 1)
    with pytest.raises(ValueError, match="imputation 2, s: the taxa draw holds"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_collapse_is_refused_for_a_panel_counted_on_all_reads(tmp_path):
    from unclassified_meta import panel
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t")
    q = tmp_path / "quant"
    pc = quant(tmp_path, [1, 0, 0], 1)
    rec = json.loads((q / "provenance.json").read_text())
    panel.write_quant_provenance(q, {**{k: v for k, v in rec.items() if k not in ("artifacts", "sealed")},
                                     "reads_scope": "all"})
    tx = tmp_path / "ptx.tsv"
    pd.DataFrame({"feature_id": ["ns|V1|like"], "taxid": ["562"]}).to_csv(tx, sep="\t", index=False)
    with pytest.raises(ValueError, match="needs a panel counted on the unclassified reads"):
        taxa.merge_with_panel(out, pc, None, out / "m", tx, nodes, names, ("Viruses",))


def test_a_sample_named_NA_survives_the_imputation_ledgers(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"NA": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t",
                                   assign="em", imputations=2)
    q = tmp_path / "quant"
    q.mkdir()
    ids = ["ns|V1|like"]
    pd.DataFrame({"feature_id": ids, "NA": [1]}).to_csv(q / "counts.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ids, "label": ["V1"], "identity_class": ["like"]}).to_csv(q / "features.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": ["NA"], "input_fragments": [1]}).to_csv(q / "ledger.tsv", sep="\t", index=False)
    (q / "imputations").mkdir()
    extra = []
    for m in (1, 2):
        pd.DataFrame({"feature_id": ids, "NA": [1]}).to_csv(q / f"imputations/counts_{m:03d}.tsv", sep="\t", index=False)
        pd.DataFrame({"sample": ["NA"], "assigned_to_features": [1]}).to_csv(q / f"imputations/ledger_{m:03d}.tsv",
                                                                             sep="\t", index=False)
        extra += [f"imputations/counts_{m:03d}.tsv", f"imputations/ledger_{m:03d}.tsv"]
    from unclassified_meta import panel
    panel.write_quant_provenance(q, {"reads_scope": "unclassified", "imputations": 2}, extra)
    taxa.merge_with_panel(out, q / "counts.tsv", None, out / "merged")
    assert (out / "merged/imputations/counts_002.tsv").exists()


def test_an_unsealed_mi_results_is_not_deleted(tmp_path):
    (tmp_path / "mi_results.tsv").write_text("mine\n")
    manifest.remove_stale(tmp_path, {}, [])
    assert (tmp_path / "mi_results.tsv").exists()


def test_quote_characters_in_labels_are_data(tmp_path):
    f = tmp_path / "r.tsv"
    f.write_text('feature_id\tlabel\tlfc_groupB\nf1\t"candidate one\t1.0\nf2\tspecies"\t2.0\n')
    r = da.read_results(f)
    assert list(r["feature_id"]) == ["f1", "f2"] and list(r["lfc_groupB"]) == [1.0, 2.0]


def test_huge_whole_counts_are_compared_exactly():
    assert not taxa._same_whole(2 ** 53, 2 ** 53 + 1) and taxa._same_whole(np.int64(7), "7")
