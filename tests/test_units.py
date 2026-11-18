"""Unit tests that need no external tools."""
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from unclassified_meta import assemble, classify, config, extract, readback, summarize
from unclassified_meta.io_utils import normalise_read_id


def test_normalise_read_id_strips_mate_suffix_and_comment():
    assert normalise_read_id("@read1/1 extra") == "read1"
    assert normalise_read_id("read1/2") == "read1"
    assert normalise_read_id("@r/3") == "r/3"


def test_kraken2_parser_counts_classified_and_unclassified(tmp_path):
    f = tmp_path / "k.out"
    f.write_text("C\tr1\t2\t100\t2:66\nU\tr2\t0\t100\t0:66\nU\tr3/1\t0\t100\t0:66\n")
    ids, total, classified = extract.parse_kraken2(f)
    assert ids == {"r2", "r3"} and total == 3 and classified == 1


def test_kraken2_parser_rejects_report_file(tmp_path):
    f = tmp_path / "k.report"
    f.write_text(" 50.00\t10\t10\tU\t0\tunclassified\n")
    with pytest.raises(ValueError, match="per-read"):
        extract.parse_kraken2(f)


def test_centrifuge_parser_read_with_any_hit_is_classified(tmp_path):
    f = tmp_path / "c.out"
    f.write_text("readID\tseqID\ttaxID\tscore\n"
                 "a\tunclassified\t0\t0\n"
                 "b\tNC_1\t562\t900\n"
                 "c\tunclassified\t0\t0\n"
                 "c\tNC_2\t10239\t500\n")
    ids, total, classified = extract.parse_centrifuge(f)
    assert ids == {"a"} and total == 3 and classified == 2


def test_filter_fastq_keeps_listed_reads(tmp_path):
    src = tmp_path / "r.fq"
    src.write_text("@a/1\nACGT\n+\nIIII\n@b/1\nCCCC\n+\nIIII\n")
    dst = tmp_path / "o.fq.gz"
    assert extract.filter_fastq(src, dst, {"b"}) == 1


def test_merged_cover_merges_overlaps_and_reverse_intervals():
    assert classify.merged_cover([(1, 100), (50, 150), (300, 201)]) == 250
    assert classify.merged_cover([]) == 0


def test_parse_blast_cover_and_thresholds(tmp_path):
    f = tmp_path / "b.tsv"
    f.write_text("c1\tref\t95.0\t200\t400\t1\t200\t1e-50\t300\tref title\n"
                 "c1\tref\t90.0\t150\t400\t201\t350\t1e-30\t200\tref title\n"
                 "c2\tref\t70.0\t300\t300\t1\t300\t1e-40\t250\tref title\n")
    hits = classify.parse_blast(f, min_identity=85)
    assert set(hits) == {"c1"}
    assert hits["c1"]["cover"] == pytest.approx(350 / 400)
    assert hits["c1"]["title"] == "ref title"


def test_assign_by_priority_first_layer_wins():
    layers = [("rRNA", {"a": {}}), ("host:x", {"a": {}, "b": {}}), ("virus", {"b": {}, "c": {}})]
    out = classify.assign_by_priority(["a", "b", "c", "d"], layers)
    assert [out[k]["category"] for k in "abcd"] == ["rRNA", "host:x", "virus", "unclassified"]


def test_kraken_category_walks_lineage():
    parent = {"1": "1", "2759": "1", "33090": "2759", "3398": "33090", "10239": "1", "12345": "10239"}
    assert classify.kraken_category("3398", parent) == "plant"
    assert classify.kraken_category("12345", parent) == "virus"
    assert classify.kraken_category("999", parent) == "other"


def test_category_counts_adds_not_assembled():
    out = readback.category_counts(Counter({"c1": 6, "c2": 2}), {"c1": "virus", "c2": "host:h"}, 10)
    assert out == {"virus": 6, "host:h": 2, "not_assembled": 2}
    with pytest.raises(ValueError):
        readback.category_counts(Counter({"c1": 11}), {"c1": "virus"}, 10)


def test_read_category_table_fractions_sum_to_one_and_scale():
    inputs = pd.DataFrame([{"sample": "S", "reads_used": 10, "unclassified_reads": 100}])
    t = summarize.read_category_table(inputs, {"S": {"virus": 6, "host:h": 2, "not_assembled": 2}})
    fr = [c for c in t.columns if c.startswith("frac_")]
    assert t[fr].sum(axis=1).iloc[0] == pytest.approx(1.0)
    assert t["est_virus"].iloc[0] == 60


def test_memory_warning_when_estimate_exceeds_limit():
    assert assemble.memory_warning("rnaspades", 1_000, 64) is None
    w = assemble.memory_warning("rnaspades", 500_000_000, 16)
    assert w and "--subsample" in w


def test_subsample_parsing():
    assert config.parse_subsample("all") is None
    assert config.parse_subsample("2e6") == 2_000_000
    assert config.load_config(None)["subsample"] == 1_000_000


def test_detect_samples_pairs_and_classifier_outputs(tmp_path):
    for n in ("A_R1.fastq.gz", "A_R2.fastq.gz", "B_1.fq.gz", "B_2.fq.gz", "C.fastq.gz"):
        (tmp_path / n).write_bytes(b"")
    (tmp_path / "A.kraken2").write_text("U\tr\t0\t1\t0:1\n")
    (tmp_path / "A.report").write_text(" 1.0\t1\t1\tU\t0\tunclassified\n")
    (tmp_path / "B.centrifuge").write_text("readID\tseqID\ttaxID\tscore\n")
    s = {x.name: x for x in config.detect_samples(tmp_path)}
    assert set(s) == {"A", "B", "C"}
    assert s["A"].paired and s["A"].classifier == "kraken2"
    assert s["B"].classifier == "centrifuge" and s["C"].classifier == "none" and not s["C"].paired


def test_missing_references_fail_early(tmp_path):
    cfg = config.load_config(None, {"references": {"refs_dir": str(tmp_path)}})
    with pytest.raises(ValueError, match="download-refs"):
        config.check_references(cfg)
