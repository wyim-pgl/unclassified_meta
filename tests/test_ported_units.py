"""Unit tests for extract thresholds, taxa-table, merge modes and control reads."""
import gzip
import json

import pandas as pd
import pytest

from unclassified_meta import control, extract, panel, taxa
from unclassified_meta.config import Sample
from unclassified_meta.pipeline import extract_all

CFR_HEADER = "readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n"


def cfr(rid, taxid, score, hitlen, seq="NC_1"):
    return f"{rid}\t{seq if taxid else 'unclassified'}\t{taxid}\t{score}\t0\t{hitlen}\t200\t1\n"


# ------------------------------------------------------------------ extract --min-score / --min-hitlen

def test_read_with_several_lines_is_classified_when_any_line_passes(tmp_path):
    f = tmp_path / "c.tsv"
    f.write_text(CFR_HEADER
                 + cfr("a", 0, 0, 0)                               # no hit
                 + cfr("b", 562, 900, 80)                          # strong
                 + cfr("c", 562, 49, 22) + cfr("c", 1280, 400, 60)  # -k 2: one weak line, one passing line
                 + cfr("d", 562, 49, 22) + cfr("d", 1280, 64, 23)   # -k 2: both below
                 + cfr("e", 562, 900, 20))                         # score fine, hit too short
    plain = extract.parse_centrifuge(f)
    assert plain.unclassified == {"a"} and (plain.total, plain.classified, plain.low_score) == (5, 4, 0)
    thr = extract.parse_centrifuge(f, min_score=100, min_hitlen=30)
    assert thr.unclassified == {"a", "d", "e"} and (thr.total, thr.classified, thr.low_score) == (5, 2, 2)
    assert thr.classified + thr.low_score + (len(thr.unclassified) - thr.low_score) == thr.total
    only_len = extract.parse_centrifuge(f, min_hitlen=30)
    assert only_len.unclassified == {"a", "d", "e"} and only_len.low_score == 2
    edge = extract.parse_centrifuge(f, min_score=49, min_hitlen=22)       # thresholds are inclusive, as in quant
    assert edge.classified == 3 and edge.unclassified == {"a", "e"}       # c and d pass exactly; e has a 20-long hit


def test_thresholds_need_the_score_columns_and_do_not_apply_to_kraken2(tmp_path):
    short = tmp_path / "short.tsv"
    short.write_text("readID\tseqID\ttaxID\na\tNC_1\t562\n")
    assert extract.parse_centrifuge(short).classified == 1
    with pytest.raises(ValueError, match="6 or more columns"):
        extract.parse_centrifuge(short, min_score=10)
    k = tmp_path / "k.out"
    k.write_text("C\tr1\t2\t100\t2:66\n")
    with pytest.raises(ValueError, match="Kraken2"):
        extract.parse_classifier(k, "kraken2", min_score=10)
    assert extract.parse_classifier(k, "kraken2").low_score == 0


def fastq(path, ids):
    with gzip.open(path, "wt") as fh:
        for i in ids:
            fh.write(f"@{i}\nACGTACGT\n+\nIIIIIIII\n")


def test_extract_all_keeps_the_accounting_closed_and_refuses_kraken2_with_thresholds(tmp_path):
    fastq(tmp_path / "s.fastq.gz", "abcde")
    c = tmp_path / "s.cfr"
    c.write_text(CFR_HEADER + cfr("a", 0, 0, 0) + cfr("b", 562, 900, 80) + cfr("c", 562, 49, 22)
                 + cfr("d", 562, 49, 22) + cfr("e", 562, 900, 80))
    s = Sample("s", tmp_path / "s.fastq.gz", None, c, "centrifuger")
    row = extract_all([s], tmp_path / "out", min_score=100).iloc[0]
    assert (row.total_reads, row.classified_reads, row.low_score_reads, row.no_hit_reads, row.unclassified_reads) \
        == (5, 2, 2, 1, 3)
    again = extract_all([s], tmp_path / "out", min_score=0).iloc[0]          # other thresholds: recomputed, not cached
    assert again.low_score_reads == 0 and again.unclassified_reads == 1
    k = tmp_path / "k.out"
    k.write_text("".join(f"U\t{i}\t0\t8\t0:1\n" for i in "abcde"))
    with pytest.raises(ValueError, match="Kraken2 output: k"):
        extract_all([Sample("k", tmp_path / "s.fastq.gz", None, k, "kraken2")], tmp_path / "o2", min_hitlen=20)
    with pytest.raises(ValueError, match="no sample has Centrifuge"):
        extract_all([Sample("n", tmp_path / "s.fastq.gz")], tmp_path / "o3", min_score=5)


# ------------------------------------------------------------------ taxa-table

K2 = """ 20.00\t20\t20\tU\t0\tunclassified
 80.00\t80\t1\tR\t1\troot
 70.00\t70\t4\tD\t2\t  Bacteria
 60.00\t60\t6\tG\t561\t    Escherichia
 40.00\t40\t30\tS\t562\t      Escherichia coli
 10.00\t10\t10\tS1\t83333\t        Escherichia coli K-12
 14.00\t14\t14\tS\t564\t      Escherichia fergusonii
  6.00\t6\t6\tG\t590\t    Salmonella
  9.00\t9\t0\tD\t10239\t  Viruses
  9.00\t9\t9\tS\t70001\t    Test phage
"""


def write(path, text):
    path.write_text(text)
    return path


def test_species_table_uses_clade_reads_and_every_read_has_one_outcome(tmp_path):
    r = taxa.parse_report(write(tmp_path / "a.kreport", K2), "S", report_format="kraken2")
    assert r.counts == {"562": 40, "564": 14, "70001": 9}           # 562 includes its strain (clade reads)
    assert r.ledger == {"input_fragments": 100, "unclassified": 20, "below_threshold": 0, "dropped_lineages": 0,
                        "kept_at_rank": 63, "unresolved_above_rank": 17, "rounding_residual": 0, "report_rows": 10,
                        "report_kind": "kraken2"}
    assert r.lineages["562"] == "root;Bacteria;Escherichia;Escherichia coli"
    g = taxa.parse_report(tmp_path / "a.kreport", "G", report_format="kraken2")
    assert g.counts == {"561": 60, "590": 6} and g.ledger["unresolved_above_rank"] == 14   # root 1 + Bacteria 4 + phage 9


def test_dropped_lineage_takes_its_whole_subtree_by_name_or_taxid(tmp_path):
    f = write(tmp_path / "a.kreport", K2)
    by_name = taxa.parse_report(f, "S", drop=("viruses",), report_format="kraken2")
    by_id = taxa.parse_report(f, "S", drop=("10239",), report_format="kraken2")
    for r in (by_name, by_id):
        assert "70001" not in r.counts and r.ledger["dropped_lineages"] == 9 and r.ledger["kept_at_rank"] == 54
        assert sum(r.ledger[c] for c in taxa.LEDGER_COLS[1:]) == 100
    below = taxa.parse_report(f, "G", drop=("Escherichia coli",), report_format="kraken2")    # a lineage below the chosen rank
    assert below.counts == {"561": 20, "590": 6} and below.ledger["dropped_lineages"] == 40
    assert sum(below.ledger[c] for c in taxa.LEDGER_COLS[1:]) == 100


def test_report_with_minimizer_columns_and_centrifuger_format_3(tmp_path):
    k8 = "".join("\t".join(f[:3] + ["500", "400"] + f[3:]) + "\n"
                 for f in (line.split("\t") for line in K2.splitlines()))
    r8 = taxa.parse_report(write(tmp_path / "m.kreport", k8), "S", report_format="kraken2")
    assert r8.counts == {"562": 40, "564": 14, "70001": 9} and r8.ledger["rounding_residual"] == 0
    # centrifuger-quant --output-format 3: first column = estimated abundance (not a read share), no
    # 'unclassified' row, non-canonical ranks as <letter><distance>, every row rounded
    cq = ("100.00\t101\t0\tR\t1\troot\n"
          "100.00\t101\t0\tR1\t131567\t  cellular organisms\n"
          "97.10\t91\t2\tD\t2\t    Bacteria\n"
          "80.00\t45\t45\tS\t562\t      Escherichia coli\n"
          "17.10\t45\t45\tS\t564\t      Escherichia fergusonii\n"      # 44.5 + 44.5 + 2 = 91: rows rounded up
          "2.90\t10\t10\tS\t9606\t    Homo sapiens\n")
    rc = taxa.parse_report(write(tmp_path / "q.kreport", cq), "S", drop=("Homo sapiens",), report_format="centrifuger-quant")
    assert rc.counts == {"562": 45, "564": 45}
    assert rc.ledger["unclassified"] == 0 and rc.ledger["dropped_lineages"] == 10
    assert rc.ledger["rounding_residual"] == -1 and sum(rc.ledger[c] for c in taxa.LEDGER_COLS[1:]) == 101


def test_centrifuger_kreport_layout_and_default_quant_format(tmp_path):
    # centrifuger-kreport: non-canonical ranks are '-', the unclassified row comes first
    ck = ("  0.12\t10\t10\tU\t0\tunclassified\n"
          " 99.88\t90\t5\t-\t1\troot\n"
          " 90.00\t85\t15\t-\t131567\t  cellular_organisms\n"
          " 70.00\t70\t10\tG\t561\t    Escherichia\n"
          " 60.00\t60\t60\tS\t562\t      Escherichia_coli\n")
    r = taxa.parse_report(write(tmp_path / "c.kreport", ck), "S", report_format="centrifuger-kreport")
    assert r.counts == {"562": 60} and r.ledger["unclassified"] == 10 and r.ledger["unresolved_above_rank"] == 30
    assert sum(r.ledger[c] for c in taxa.LEDGER_COLS[1:]) == 100
    # the default centrifuger-quant table (format 0) is not a Kraken-style report
    f0 = ("name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n"
          "Escherichia coli\t562\tspecies\t4600000\t50\t48\t0.39\n")
    with pytest.raises(ValueError, match="not a Kraken-style report row"):
        taxa.parse_report(write(tmp_path / "f0.tsv", f0), "S", report_format="kraken2")


def test_broken_reports_are_refused(tmp_path):
    with pytest.raises(ValueError, match="per-read output"):
        taxa.parse_report(write(tmp_path / "p.tsv", "C\tr1\t562\t100\t562:66\n"), "S", report_format="kraken2")
    with pytest.raises(ValueError, match="clade reads|does not add up"):
        taxa.parse_report(write(tmp_path / "b.kreport", " 100\t100\t0\tR\t1\troot\n 10\t10\t10\tS\t5\t  x\n"), "S", report_format="kraken2")
    with pytest.raises(ValueError, match="whole numbers"):
        taxa.parse_report(write(tmp_path / "f.kreport", " 100\t10.5\t10.5\tS\t5\tx\n"), "S", report_format="kraken2")
    with pytest.raises(ValueError, match="rank must be one of"):
        taxa.parse_report(write(tmp_path / "ok.kreport", K2), "X", report_format="kraken2")


def test_suffix_that_matches_two_reports_of_one_sample_is_an_error(tmp_path):
    for name in ("s1.kreport", "s1.quant.kreport", "s2.kreport", "s2.quant.kreport"):
        write(tmp_path / name, K2)
    with pytest.raises(ValueError, match="more than one report for sample 's1'"):
        taxa.find_reports(tmp_path, ".kreport")
    assert sorted(taxa.find_reports(tmp_path, ".quant.kreport")) == ["s1", "s2"]
    with pytest.raises(ValueError, match="give --suffix"):
        taxa.find_reports(tmp_path, "")
    with pytest.raises(ValueError, match="no file ends with"):
        taxa.find_reports(tmp_path, ".report")


def test_report_sheet_refuses_duplicate_samples_and_reports(tmp_path):
    write(tmp_path / "a.kreport", K2)
    write(tmp_path / "b.kreport", K2)
    with pytest.raises(ValueError, match="duplicate sample"):
        taxa.read_report_sheet(write(tmp_path / "d.tsv", "sample\treport\ns\ta.kreport\ns\tb.kreport\n"))
    with pytest.raises(ValueError, match="same report"):
        taxa.read_report_sheet(write(tmp_path / "e.tsv", "sample\treport\ns\ta.kreport\nt\ta.kreport\n"))
    ok = taxa.read_report_sheet(write(tmp_path / "ok.tsv", "sample\treport\ns\ta.kreport\nt\tb.kreport\n"))
    assert list(ok) == ["s", "t"]


def two_reports(tmp_path):
    write(tmp_path / "s1.kreport", K2)
    write(tmp_path / "s2.kreport", K2.replace("\t      Escherichia fergusonii", "\t      Escherichia fergusonii")
          .replace(" 14.00\t14\t14\tS\t564", " 14.00\t14\t14\tS\t565"))      # s2 has another species instead
    return taxa.find_reports(tmp_path, ".kreport")


def test_table_has_namespaced_ids_integer_counts_and_a_ledger_that_matches(tmp_path):
    out = tmp_path / "taxa"
    taxa.build_table(two_reports(tmp_path), out, "S", ("Viruses",), "cfr_test", report_format="kraken2")
    counts = pd.read_csv(out / "counts.tsv", sep="\t").set_index("feature_id")
    feats = pd.read_csv(out / "features.tsv", sep="\t")
    led = pd.read_csv(out / "ledger.tsv", sep="\t").set_index("sample")
    assert list(counts.index) == ["cfr_test|562|S", "cfr_test|564|S", "cfr_test|565|S"]
    assert counts.loc["cfr_test|564|S"].tolist() == [14, 0] and (counts.dtypes == "int64").all()
    assert (counts.sum(axis=0) == led.loc[counts.columns, "kept_at_rank"]).all()
    assert (led[list(taxa.LEDGER_COLS[1:])].sum(axis=1) == led["input_fragments"]).all()
    assert feats["identity_class"].eq("taxon").all() and feats.loc[0, "label"] == "Escherichia coli"
    with pytest.raises(ValueError, match="not found in any report: archaea"):
        taxa.build_table(two_reports(tmp_path), tmp_path / "t2", "S", ("Archaea",), report_format="kraken2")


def panel_dir(tmp_path, scope, samples=("s1", "s2"), in_frag=(30, 30), summary=None):
    q = tmp_path / f"quant_{scope}"
    q.mkdir()
    pd.DataFrame({"feature_id": ["ns|V1|like", "ns|V2|same"], samples[0]: [10, 5], samples[1]: [8, 2]}) \
        .to_csv(q / "counts.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ["ns|V1|like", "ns|V2|same"], "label": ["V1-like", "V2"],
                  "identity_class": ["like", "same"]}).to_csv(q / "features.tsv", sep="\t", index=False)
    pd.DataFrame({"sample": samples, "input_fragments": in_frag}).to_csv(q / "ledger.tsv", sep="\t", index=False)
    panel.write_quant_provenance(q, {"reads_scope": scope, "input_summary": summary})
    return q / "counts.tsv"


def test_merge_with_a_panel_counted_on_all_reads_needs_the_panel_lineage(tmp_path):
    pc = panel_dir(tmp_path, "all")
    with pytest.raises(ValueError, match="--panel-lineage"):
        taxa.check_merge_mode(pc, None)
    assert taxa.check_merge_mode(pc, "Viruses") == ["Viruses"]
    assert taxa.check_merge_mode(pc, "none") == []
    out = tmp_path / "taxa"
    taxa.build_table(two_reports(tmp_path), out, "S", ("Viruses",), report_format="kraken2")
    merged = taxa.merge_with_panel(out, pc, None, out / "merged")
    counts = pd.read_csv(merged, sep="\t").set_index("feature_id")
    assert list(counts.index[:2]) == ["ns|V1|like", "ns|V2|same"] and len(counts) == 5
    feats = pd.read_csv(out / "merged" / "features.tsv", sep="\t")
    assert feats["source"].tolist() == ["panel"] * 2 + ["background"] * 3
    rec = pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t")
    assert rec["status"].str.startswith("not_applicable").all()


def test_merge_refuses_other_samples(tmp_path):
    pc = panel_dir(tmp_path, "all", samples=("s1", "s3"))
    out = tmp_path / "taxa"
    taxa.build_table(two_reports(tmp_path), out, "S", ("Viruses",), report_format="kraken2")
    with pytest.raises(ValueError, match="only in the panel: s3; only in the reports: s2"):
        taxa.merge_with_panel(out, pc, None, out / "merged")


def test_disjoint_mode_needs_no_lineage_and_reconciles_to_the_input_reads(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [80, 80],
                  "unclassified_reads": [20, 20]}).to_csv(summ, sep="\t", index=False)
    pc = panel_dir(tmp_path, "unclassified", in_frag=(20, 20), summary=str(summ))
    assert taxa.check_merge_mode(pc, None) == []                      # overlapping lineages are allowed
    out = tmp_path / "taxa"
    taxa.build_table(two_reports(tmp_path), out, "S", report_format="kraken2")
    taxa.merge_with_panel(out, pc, None, out / "merged")
    rec = pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t").set_index("sample")
    assert (rec["status"] == "ok").all() and (rec["difference"] == 0).all()
    r = rec.loc["s1"]
    assert r.classified_kept_at_rank + r.classified_dropped_lineages + r.classified_unresolved_above_rank \
        + r.classified_rounding_residual + r.panel_features + r.panel_rest == r.total_input_reads == 100
    counts = pd.read_csv(out / "merged" / "counts.tsv", sep="\t").set_index("feature_id")
    assert "taxid|70001|S" in counts.index                            # the classifier's own virus stays


def test_disjoint_mode_fails_when_the_two_sides_do_not_add_up(tmp_path):
    summ = tmp_path / "input_summary.tsv"
    pd.DataFrame({"sample": ["s1", "s2"], "total_reads": [100, 100], "classified_reads": [70, 70],
                  "unclassified_reads": [30, 30]}).to_csv(summ, sep="\t", index=False)     # extract used a threshold
    pc = panel_dir(tmp_path, "unclassified", in_frag=(30, 30), summary=str(summ))           # the report did not
    out = tmp_path / "taxa"
    taxa.build_table(two_reports(tmp_path), out, "S", report_format="kraken2")
    with pytest.raises(ValueError, match="do not add up"):
        taxa.merge_with_panel(out, pc, None, out / "merged")
    assert (pd.read_csv(out / "merged" / "reconciliation.tsv", sep="\t")["status"] == "MISMATCH").all()
    # without the extract summary the report's own unclassified row is the reference
    pc2 = panel_dir(tmp_path / "x" if (tmp_path / "x").mkdir() is None else tmp_path, "unclassified", in_frag=(20, 20))
    taxa.merge_with_panel(out, pc2, None, out / "merged2")
    assert (pd.read_csv(out / "merged2" / "reconciliation.tsv", sep="\t")["status"] == "ok").all()


def test_reads_scope_is_recognised_from_the_extract_folder(tmp_path):
    un = tmp_path / "res" / "unclassified_reads"
    (tmp_path / "res" / "tables").mkdir(parents=True)
    un.mkdir()
    fastq(un / "s_R1.fastq.gz", "ab")
    s = [Sample("s", un / "s_R1.fastq.gz")]
    assert panel.detect_reads_scope(s) == ("all", None)               # no input summary: not an extract output
    (tmp_path / "res" / "tables" / "input_summary.tsv").write_text("sample\ttotal_reads\n")
    scope, summ = panel.detect_reads_scope(s)
    assert scope == "unclassified" and summ.name == "input_summary.tsv"
    fastq(tmp_path / "t.fastq.gz", "ab")
    assert panel.detect_reads_scope([Sample("t", tmp_path / "t.fastq.gz")]) == ("all", None)


# ------------------------------------------------------------------ reversed-read control

def test_control_reads_are_reversed_not_complemented(tmp_path):
    with gzip.open(tmp_path / "r.fastq.gz", "wt") as fh:
        fh.write("@r1/1 extra\nAACGT\n+\nABCDE\n@r2/1\nGGGTC\n+\nIIIII\n")
    n = control.reverse_fastq(tmp_path / "r.fastq.gz", tmp_path / "c.fastq.gz", None, 11)
    lines = gzip.open(tmp_path / "c.fastq.gz", "rt").read().splitlines()
    assert n == 2 and lines[:4] == ["@CTRL_r1/1 extra", "TGCAA", "+", "EDCBA"]     # TGCAA, not the complement ACGTT
    assert lines[5] == "CTGGG"


def test_make_controls_writes_paired_files_and_a_sheet(tmp_path):
    fastq(tmp_path / "s_R1.fastq.gz", ["a/1", "b/1"])
    fastq(tmp_path / "s_R2.fastq.gz", ["a/2", "b/2"])
    df = control.make_controls([Sample("s", tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz")],
                               tmp_path / "out", None)
    assert df.loc[0, "control_reads"] == 2 and df.loc[0, "sample"] == "s_ctrl"
    assert (tmp_path / "out" / "control_reads" / "s_ctrl_R2.fastq.gz").exists()
    assert (tmp_path / "out" / "control_samples.tsv").exists()


def control_output(tmp_path):
    f = tmp_path / "ctrl.cfr.tsv"
    lines = [CFR_HEADER]
    for i in range(100):                                   # control reads: 10 weak chance hits, 1 strong one
        if i < 10:
            lines.append(cfr(f"CTRL_c{i}", 562, 25, 20))
        elif i == 10:
            lines.append(cfr(f"CTRL_c{i}", 562, 16, 19) + cfr(f"CTRL_c{i}", 1280, 121, 26))    # -k 2
        else:
            lines.append(cfr(f"CTRL_c{i}", 0, 0, 0))
    for i in range(50):                                    # real reads: 30 strong, 10 weak, 10 unclassified
        lines.append(cfr(f"r{i}", 562, 400, 60) if i < 30 else cfr(f"r{i}", 562, 25, 20) if i < 40
                     else cfr(f"r{i}", 0, 0, 0))
    f.write_text("".join(lines))
    return f


def test_control_summary_reports_chance_hits_per_threshold_and_a_recommendation(tmp_path):
    df = control.summarise([control_output(tmp_path)], tmp_path / "out", scores=(0, 25, 36, 121, 200),
                           hitlens=(0, 22), target_pct=1.0)
    cell = df.set_index(["min_hitlen", "min_score"])
    assert cell.loc[(0, 0), "control_classified_pct"] == 11.0 and cell.loc[(0, 0), "real_classified_pct"] == 80.0
    assert cell.loc[(0, 25), "control_classified"] == 11              # thresholds are inclusive
    assert cell.loc[(0, 36), "control_classified"] == 1 and cell.loc[(0, 121), "control_classified"] == 1
    assert cell.loc[(0, 200), "control_classified"] == 0
    assert cell.loc[(22, 0), "control_classified"] == 1               # only the 26-long hit of the -k 2 read is left
    assert cell.loc[(0, 36), "real_classified_pct"] == 60.0
    assert cell.loc[(0, 0), "real_beyond_chance_pct"] == pytest.approx(100 * (80 - 11) / (100 - 11), abs=1e-3)
    best = df[df["recommended"]].iloc[0]
    assert (best.min_score, best.min_hitlen) == (0, 22) and best.control_classified_pct == 1.0
    for f in ("control_summary.tsv", "control_summary.png", "control_summary.pdf"):
        assert (tmp_path / "out" / f).exists()


def test_control_summary_needs_control_reads_and_handles_kraken2(tmp_path):
    f = tmp_path / "real.tsv"
    f.write_text(CFR_HEADER + cfr("r1", 562, 100, 30))
    with pytest.raises(ValueError, match="no control reads"):
        control.summarise([f], tmp_path / "o1")
    assert control.summarise([f], tmp_path / "o2", all_control=True)["control_classified_pct"].iloc[0] == 100.0
    k = tmp_path / "k.out"
    k.write_text("C\tCTRL_a\t2\t100\t2:66\nU\tCTRL_b\t0\t100\t0:66\nU\tCTRL_c\t0\t100\t0:66\nU\tCTRL_d\t0\t100\t0:66\n")
    dk = control.summarise([k], tmp_path / "o3", classifier="kraken2")
    assert len(dk) == 1 and dk["control_classified_pct"].iloc[0] == 25.0 and not dk["recommended"].any()


# ------------------------------------------------------------------ taxa-table from per-read output + taxonomy

NODES = """1\t|\t1\t|\tno rank\t|\t\t|
131567\t|\t1\t|\tno rank\t|\t\t|
2\t|\t131567\t|\tsuperkingdom\t|\t\t|
543\t|\t2\t|\tfamily\t|\t\t|
561\t|\t543\t|\tgenus\t|\t\t|
562\t|\t561\t|\tspecies\t|\t\t|
83333\t|\t562\t|\tstrain\t|\t\t|
564\t|\t561\t|\tspecies\t|\t\t|
590\t|\t543\t|\tgenus\t|\t\t|
28901\t|\t590\t|\tspecies\t|\t\t|
10239\t|\t1\t|\tsuperkingdom\t|\t\t|
70001\t|\t10239\t|\tspecies\t|\t\t|
"""
NAMES = "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in (
    (1, "root"), (131567, "cellular organisms"), (2, "Bacteria"), (543, "Enterobacteriaceae"), (561, "Escherichia"),
    (562, "Escherichia coli"), (83333, "Escherichia coli K-12"), (564, "Escherichia fergusonii"),
    (590, "Salmonella"), (28901, "Salmonella enterica"), (10239, "Viruses"), (70001, "Test phage"))) \
    + "562\t|\tE. coli\t|\t\t|\tsynonym\t|\n"


def taxonomy_files(tmp_path):
    return write(tmp_path / "nodes.dmp", NODES), write(tmp_path / "names.dmp", NAMES)


def per_read_file(tmp_path):
    f = tmp_path / "s.cfr.tsv"
    f.write_text(CFR_HEADER
                 + cfr("u1", 0, 0, 0) + cfr("u2", 0, 0, 0)                        # unclassified
                 + cfr("low1", 562, 48, 30) + cfr("low2", 562, 400, 21)           # just below one threshold each
                 + cfr("edge", 562, 49, 22)                                       # exactly at both thresholds: kept
                 + cfr("strain", 83333, 400, 60)                                  # below species: projected to 562
                 + cfr("sp", 564, 400, 60)
                 + cfr("genus", 561, 400, 60)                                     # above species: unresolved
                 + cfr("k_same", 562, 400, 60) + cfr("k_same", 83333, 400, 60)    # k=2 within one species -> 562
                 + cfr("k_genus", 562, 400, 60) + cfr("k_genus", 564, 400, 60)    # k=2, two species -> LCA genus
                 + cfr("k_fam", 562, 400, 60) + cfr("k_fam", 28901, 400, 60) + cfr("k_fam", 564, 400, 60)
                 + cfr("k_weak", 564, 400, 60) + cfr("k_weak", 28901, 10, 16)     # weak second line is ignored
                 + cfr("phage", 70001, 400, 60)
                 + cfr("alien", 999999, 400, 60))                                 # taxID missing from the taxonomy
    return f


def test_taxonomy_lca_and_both_dump_formats(tmp_path):
    nodes, _ = taxonomy_files(tmp_path)
    tax = taxa.Taxonomy(nodes)
    assert tax.path("83333") == ("83333", "562", "561", "543", "2", "131567", "1")
    assert tax.lca(["562", "83333"]) == "562" and tax.lca(["562", "564"]) == "561"
    assert tax.lca(["562", "28901"]) == "543" and tax.lca(["562", "70001"]) == "1" and tax.lca(["562", "5"]) is None
    inspect = write(tmp_path / "inspect.dmp", "1\t|\t1\t|\tno rank\t|\n561\t|\t1\t|\tgenus\t|\n562\t|\t561\t|\tspecies\t|\n")
    t2 = taxa.Taxonomy(inspect)                       # `centrifuger-inspect --taxonomy-tree` layout
    assert t2.path("562") == ("562", "561", "1") and t2.rank["562"] == "species"
    ids, by_name = taxa.scan_names(write(tmp_path / "n.dmp", "562\t|\tEscherichia coli\t|\tscientific name\t|\n"),
                                   {"562"}, {"escherichia coli"})
    assert ids == {"562": "Escherichia coli"} and by_name == {"escherichia coli": {"562"}}


def test_per_read_table_projects_to_the_rank_and_closes_the_ledger(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    f = per_read_file(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses",), "ncbi_test", 49, 22)
    counts = pd.read_csv(out / "counts.tsv", sep="\t").set_index("feature_id")["s"].to_dict()
    assert counts == {"ncbi_test|562|S": 3, "ncbi_test|564|S": 2}      # 562: edge, strain, k_same; 564: sp, k_weak
    led = pd.read_csv(out / "ledger.tsv", sep="\t").iloc[0]
    assert (led.input_fragments, led.unclassified, led.below_threshold, led.dropped_lineages, led.kept_at_rank,
            led.unresolved_above_rank) == (14, 2, 2, 1, 5, 4)         # unresolved: genus, k_genus, k_fam, alien
    assert led.not_in_taxonomy == 1 and led.reads_assigned_by_lca == 3
    assert sum(led[c] for c in taxa.LEDGER_COLS[1:]) == led.input_fragments
    feats = pd.read_csv(out / "features.tsv", sep="\t").set_index("feature_id")
    assert feats.loc["ncbi_test|562|S", "label"] == "Escherichia coli"           # scientific name, not the synonym
    assert feats.loc["ncbi_test|562|S", "lineage"].endswith("Enterobacteriaceae;Escherichia;Escherichia coli")
    g = tmp_path / "taxa_g"
    taxa.build_table_from_per_read({"s": f}, g, nodes, names, "G", (), "ncbi_test", 49, 22)
    cg = pd.read_csv(g / "counts.tsv", sep="\t").set_index("feature_id")["s"].to_dict()
    assert cg == {"ncbi_test|561|G": 7}                                          # + genus, k_genus; phage has no genus
    lg = pd.read_csv(g / "ledger.tsv", sep="\t").iloc[0]
    assert lg.unresolved_above_rank == 3 and lg.dropped_lineages == 0            # k_fam, phage, alien
    with pytest.raises(ValueError, match="not found in the taxonomy: archaea"):
        taxa.build_table_from_per_read({"s": f}, tmp_path / "x", nodes, names, "S", ("Archaea",), "ncbi_test")


def test_per_read_table_and_extract_split_the_reads_identically(tmp_path):
    """below_threshold + unclassified of taxa-table are exactly the reads `extract` writes, for any thresholds."""
    nodes, names = taxonomy_files(tmp_path)
    f = per_read_file(tmp_path)
    ids = [line.split("\t")[0] for line in f.read_text().splitlines()[1:]]
    fastq(tmp_path / "s.fastq.gz", list(dict.fromkeys(ids)))
    s = Sample("s", tmp_path / "s.fastq.gz", None, f, "centrifuger")
    for ms, hl in ((0, 0), (49, 22), (50, 22), (49, 23), (400, 60)):
        ex = extract_all([s], tmp_path / f"ex_{ms}_{hl}", ms, hl).iloc[0]
        out = tmp_path / f"t_{ms}_{hl}"
        taxa.build_table_from_per_read({"s": f}, out, nodes, names, "G", (), "ncbi_test", ms, hl)
        led = pd.read_csv(out / "ledger.tsv", sep="\t").iloc[0]
        assert led.below_threshold + led.unclassified == ex.unclassified_reads, (ms, hl)
        assert led.below_threshold == ex.low_score_reads and led.unclassified == ex.no_hit_reads
        assert led.kept_at_rank + led.dropped_lineages + led.unresolved_above_rank == ex.classified_reads
        assert led.input_fragments == ex.total_reads
    plain = pd.read_csv(tmp_path / "t_0_0" / "ledger.tsv", sep="\t").iloc[0]
    edge = pd.read_csv(tmp_path / "t_50_22" / "ledger.tsv", sep="\t").iloc[0]
    assert plain.below_threshold == 0 and edge.below_threshold == 3              # 'edge' (score 49) now falls below


def test_taxa_table_cli_per_read_needs_the_taxonomy_and_its_name(tmp_path):
    import subprocess
    import sys
    nodes, names = taxonomy_files(tmp_path)
    d = tmp_path / "cfr"
    d.mkdir()
    per_read_file(tmp_path).rename(d / "s1.cfr.tsv")
    base = [sys.executable, "-m", "unclassified_meta.cli", "taxa-table", "--input-dir", str(d), "--suffix", ".cfr.tsv",
            "-o", str(tmp_path / "out")]

    def run(*extra):
        p = subprocess.run(base + list(extra), capture_output=True, text=True)
        return p.returncode, p.stdout + p.stderr

    rc, out = run("--per-read")
    assert rc != 0 and "--nodes" in out
    rc, out = run("--per-read", "--nodes", str(nodes), "--names", str(names))
    assert rc != 0 and "--taxonomy-name" in out
    rc, out = run("--min-score", "49")                                     # thresholds without --per-read
    assert rc != 0 and "belong to --per-read" in out
    rc, out = run("--per-read", "--nodes", str(nodes), "--names", str(names), "--taxonomy-name", "ncbi_test",
                  "--min-score", "49", "--min-hitlen", "22", "--drop-lineage", "Viruses")
    assert rc == 0, out
    counts = pd.read_csv(tmp_path / "out" / "counts.tsv", sep="\t").set_index("feature_id")["s1"].to_dict()
    assert counts == {"ncbi_test|562|S": 3, "ncbi_test|564|S": 2}
