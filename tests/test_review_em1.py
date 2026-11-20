"""Findings of the first review of the EM / imputation / panel-taxonomy change (each failed before its fix)."""
import numpy as np
import pandas as pd
import pytest

from unclassified_meta import manifest, panel, panel_taxonomy as ptx, taxa

from test_em_merge import quant
from test_em_taxa import em_file
from test_ported_units import CFR_HEADER, cfr, taxonomy_files
from test_review_fixes import HAVE_PANEL_TOOLS, tiny_run


def test_distinct_dropped_taxa_are_distinct_candidates(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    f = tmp_path / "s.cfr.tsv"
    f.write_text(CFR_HEADER + "".join(cfr(f"x{i}", 70001, 400, 60) for i in range(100))     # phage, dropped
                 + "".join(cfr(f"a{i}", 564, 400, 60) for i in range(10))
                 + "".join(cfr(f"y{i}", 564, 400, 60) + cfr(f"y{i}", 28901, 400, 60) for i in range(10)))
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses", "Salmonella"), "t", assign="em")
    c = pd.read_csv(out / "counts.tsv", sep="\t").set_index("feature_id")["s"]
    assert c["t|564|S"] == pytest.approx(20, abs=1e-3)       # Salmonella has no read of its own, the phage is unrelated


def test_samples_with_the_same_ambiguous_reads_get_independent_imputations(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(); b.mkdir()
    reports = {"s1": em_file(a), "s2": em_file(b)}
    out = tmp_path / "t"
    taxa.build_table_from_per_read(reports, out, nodes, names, "S", ("Viruses",), "t", assign="em", imputations=6)
    tabs = [pd.read_csv(out / f"imputations/counts_{m:03d}.tsv", sep="\t").set_index("feature_id")
            for m in range(1, 7)]
    assert any((t["s1"] != t["s2"]).any() for t in tabs)


def test_merge_refuses_a_feature_table_that_misses_a_counted_feature(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t")
    pc = quant(tmp_path, [1, 0, 0], 1)
    feats = pc.parent / "features.tsv"
    pd.read_csv(feats, sep="\t").iloc[:2].to_csv(feats, sep="\t", index=False)
    panel.write_quant_provenance(pc.parent, {"reads_scope": "unclassified"})          # sealed as it is
    with pytest.raises(ValueError, match="counted but not described: ns\\|V3\\|like"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_a_side_with_em_but_no_imputations_is_refused_when_the_other_has_them(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", assign="em")
    pc = quant(tmp_path, [1, 0, 0], 1, imputations=[[1, 0, 0], [0, 1, 0]])
    with pytest.raises(ValueError, match="rerun taxa-table with --imputations 2"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_a_panel_taxonomy_edited_after_its_seal_is_refused(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t")
    pc = quant(tmp_path, [1, 0, 0], 1)
    d = tmp_path / "ptx"
    d.mkdir()
    pd.DataFrame({"feature_id": ["ns|V1|like"], "taxid": ["562"]}).to_csv(d / "panel_taxonomy.tsv", sep="\t", index=False)
    manifest.begin(d, "panel-taxonomy")
    manifest.seal(d, ["panel_taxonomy.tsv"], {"step": "panel-taxonomy"})
    taxa.merge_with_panel(out, pc, None, out / "m1", d / "panel_taxonomy.tsv", nodes, names, ("Viruses",))
    pd.DataFrame({"feature_id": ["ns|V1|like"], "taxid": ["28901"]}).to_csv(d / "panel_taxonomy.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="changed"):
        taxa.merge_with_panel(out, pc, None, out / "m2", d / "panel_taxonomy.tsv", nodes, names, ("Viruses",))


def test_rerun_removes_only_the_imputations_it_sealed_before(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "t"
    f = em_file(tmp_path)
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses",), "t", assign="em", imputations=2)
    mine = out / "imputations" / "my_notes.tsv"
    mine.write_text("kept\n")
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses",), "t", assign="em")
    assert mine.exists() and not (out / "imputations/counts_001.tsv").exists()


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_write_fasta_refuses_a_path_inside_the_panel(tmp_path):
    run, _ = tiny_run(tmp_path)
    p = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    with pytest.raises(ValueError, match="inside the panel folder"):
        ptx.write_fasta(p, p / "seq2feature.tsv")
