"""panel-quant --ambiguous em / --imputations."""
import itertools

import pandas as pd
import pytest

from unclassified_meta import fragments as fr
from unclassified_meta import manifest, panel
from unclassified_meta.config import Sample

from test_panel_units import S2F, sam
from test_ported_units import fastq
from test_review_fixes import HAVE_PANEL_TOOLS, tiny_run


RECORD_SETS = [
    [(99, "contig|a1"), (147, "contig|a1")], [(73, "contig|a1"), (133, "*")], [(65, "contig|a1"), (129, "ref|B")],
    [(0, "contig|a1", -5), (256, "ref|B", -5)], [(0, "contig|a1", 0), (256, "ref|B", -40)],
    [(0, "contig|a1", -6), (256, "decoy|host:h|chr1", -6)], [(0, "decoy|host:h|chr1")],
    [(0, "decoy|host:h|chr1", -3), (256, "decoy|rrna|r1", -3)], [(77, "*"), (141, "*")],
]


@pytest.mark.parametrize("records", RECORD_SETS)
def test_the_class_of_a_fragment_gives_the_same_ledger_outcome_as_before(records):
    lines = [sam("f", *r) for r in records]
    feats, ledger, _ = fr.count_features(lines, S2F)
    classes, n = fr.count_classes(lines, S2F)
    assert n == 1
    assert fr.outcome_of_class(next(iter(classes))) == next(iter(feats or ledger))


def fake_bowtie2(lines):
    def run(prefix, r1, r2, logfile, threads, options, consume):
        return {"result": consume(iter(lines)), "input_reads": None, "wall_s": 0.0, "max_rss_gb": 0.0}
    return run


@pytest.mark.skipif(not HAVE_PANEL_TOOLS, reason="bowtie2 / seqkit not on PATH (activate the environment)")
def test_em_shares_fragments_that_tie_between_features_and_imputations_keep_every_fragment(tmp_path, monkeypatch):
    run, _ = tiny_run(tmp_path)
    p = panel.build_panel(run, tmp_path / "panel", ["virus"], [], None, None, None, threads=1)
    s2f = panel.load_seq2feature(p)
    targets = sorted({(seq, f) for seq, f in s2f.items() if not f.startswith(fr.DECOY_PREFIX)}, key=lambda x: x[1])
    feats = sorted({f for _, f in targets})
    if len(feats) < 2:
        pytest.skip("the tiny panel has fewer than two features")
    seq_a = next(s for s, f in targets if f == feats[0])
    seq_b = next(s for s, f in targets if f == feats[1])
    names = (f"r{i}" for i in itertools.count())
    lines = ([sam(next(names), 0, seq_a) for _ in range(30)] + [sam(next(names), 0, seq_b) for _ in range(10)]
             + [l for _ in range(20) for n in [next(names)] for l in (sam(n, 0, seq_a, -2), sam(n, 256, seq_b, -2))]
             + [sam(next(names), 4, "*") for _ in range(5)])
    monkeypatch.setattr(panel, "run_bowtie2", fake_bowtie2(lines))
    fastq(tmp_path / "r.fastq.gz", "ab")
    samples = [Sample("s1", tmp_path / "r.fastq.gz")]

    q0 = tmp_path / "q0"
    panel.quantify(p, samples, q0, threads=1)
    c0 = pd.read_csv(q0 / "counts.tsv", sep="\t").set_index("feature_id")["s1"]
    l0 = pd.read_csv(q0 / "ledger.tsv", sep="\t").iloc[0]
    assert c0[feats[0]] == 30 and c0[feats[1]] == 10 and l0[fr.AMBIGUOUS] == 20      # the old policy, unchanged

    q = tmp_path / "q"
    panel.quantify(p, samples, q, threads=1, ambiguous="em", imputations=4, seed=2)
    c = pd.read_csv(q / "counts.tsv", sep="\t").set_index("feature_id")["s1"]
    led = pd.read_csv(q / "ledger.tsv", sep="\t").iloc[0]
    assert c[feats[0]] == pytest.approx(45, abs=1e-3) and c[feats[1]] == pytest.approx(15, abs=1e-3)
    assert led[fr.AMBIGUOUS] == 0 and led[fr.UNMAPPED] == 5
    outcome_cols = [c for c in led.index if c not in ("sample", "paired", "input_fragments", "wall_s", "max_rss_gb")]
    assert sum(led[c] for c in outcome_cols) == pytest.approx(led.input_fragments, abs=1e-6)
    rec = manifest.verify(q, "provenance.json")
    assert rec["em"]["s1"]["shared_fragments"] == 20 and rec["em"]["s1"]["converged"]
    assert rec["ambiguous"] == "em" and rec["imputations"] == 4
    for m in range(1, 5):
        assert f"imputations/counts_{m:03d}.tsv" in rec["artifacts"]
        ic = pd.read_csv(q / f"imputations/counts_{m:03d}.tsv", sep="\t").set_index("feature_id")["s1"]
        il = pd.read_csv(q / f"imputations/ledger_{m:03d}.tsv", sep="\t").iloc[0]
        assert ic.sum() == 60 == il.assigned_to_features and il.input_fragments == 65 and il[fr.UNMAPPED] == 5
        assert 30 <= ic[feats[0]] <= 50
    panel.quantify(p, samples, q, threads=1, ambiguous="em")          # rerun without imputations removes them
    assert not list((q / "imputations").glob("*.tsv"))


def test_panel_imputations_need_em(tmp_path):
    with pytest.raises(ValueError, match="--imputations needs --ambiguous em"):
        panel.quantify(tmp_path / "p", [], tmp_path / "q", imputations=2)
