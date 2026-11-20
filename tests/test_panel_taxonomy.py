"""panel-taxonomy: consensus taxID per panel feature."""
import pandas as pd
import pytest

from unclassified_meta import manifest, panel, panel_taxonomy as ptx, taxa

from test_ported_units import CFR_HEADER, cfr, taxonomy_files
from test_review_fixes import HAVE_PANEL_TOOLS, tiny_run


def test_consensus_takes_the_deepest_taxon_holding_most_of_the_classified_length(tmp_path):
    nodes, _ = taxonomy_files(tmp_path)
    tax = taxa.Taxonomy(nodes)
    lengths = {"c1": 700, "c2": 150, "c3": 150}
    t, sup, frac = ptx.consensus(lengths, {"c1": "83333", "c2": "562", "c3": "564"}, tax, 0.8, 0.5)
    assert t == "562" and sup == pytest.approx(0.85) and frac == 1.0          # strain 70 %, E. coli 85 %
    t, sup, _ = ptx.consensus(lengths, {"c1": "562", "c2": "564", "c3": "564"}, tax, 0.9, 0.5)
    assert t == "561"                                                          # two species: the genus holds 100 %
    t, _, frac = ptx.consensus(lengths, {"c1": None, "c2": "562", "c3": None}, tax, 0.8, 0.5)
    assert t is None and frac == pytest.approx(0.15)                           # too little classified


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_panel_taxonomy_from_classified_panel_sequences(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    run, _ = tiny_run(tmp_path)
    p = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    fa = ptx.write_fasta(p, tmp_path / "targets.fa")
    targets = ptx.target_sequences(p)
    ids = [l[1:].strip() for l in fa.read_text().splitlines() if l.startswith(">")]
    tag = panel.verify_panel(p)["panel_id"]
    assert all(i.startswith(tag + ":") for i in ids) and {i.split(":", 1)[1] for i in ids} == set(targets)
    first = sorted(set(targets.values()))[0]
    feat = lambda i: targets[i.split(":", 1)[1]]
    body = CFR_HEADER + "".join(cfr(i, 83333 if feat(i) == first else 0, 500 if feat(i) == first else 0, 60)
                                for i in ids)
    cl = tmp_path / "targets.cfr.tsv"
    cl.write_text(body)
    out = ptx.build(p, cl, nodes, names, tmp_path / "ptx")
    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False).set_index("feature_id")
    assert df.loc[first, "taxid"] == "83333" and df.loc[first, "status"] == "assigned"
    assert df.loc[first, "name"] == "Escherichia coli K-12"
    assert (df.drop(index=first)["status"] == "unclassified").all()
    assert manifest.verify(tmp_path / "ptx")["min_support"] == 0.8
    cl.write_text(CFR_HEADER + cfr(f"{tag}:not_a_panel_sequence", 562, 500, 60))
    with pytest.raises(ValueError, match="not a target sequence"):
        ptx.build(p, cl, nodes, names, tmp_path / "ptx2")
    cl.write_text(body.replace(tag + ":", tag[:12] + "0" * 20 + ":"))      # same 12-character prefix, other panel           # classified FASTA of another panel generation
    with pytest.raises(ValueError, match="not written by --write-fasta for this panel"):
        ptx.build(p, cl, nodes, names, tmp_path / "ptx3")
    assert manifest.verify(tmp_path / "ptx")["nodes_sha256"] == manifest.sha256(nodes)
