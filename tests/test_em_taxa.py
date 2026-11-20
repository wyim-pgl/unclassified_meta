"""taxa-table --per-read --assign em / --imputations."""
import json

import pandas as pd
import pytest

from unclassified_meta import manifest, taxa

from test_ported_units import CFR_HEADER, cfr, taxonomy_files


def em_file(tmp_path):
    f = tmp_path / "s.cfr.tsv"
    body = CFR_HEADER
    body += "".join(cfr(f"a{i}", 562, 400, 60) for i in range(9))        # 9 reads only E. coli
    body += cfr("b0", 564, 400, 60)                                       # 1 read only E. fergusonii
    body += "".join(cfr(f"m{i}", 562, 400, 60) + cfr(f"m{i}", 564, 400, 60) for i in range(10))   # 10 shared
    body += "".join(cfr(f"p{i}", 70001, 400, 60) for i in range(4))      # 4 phage reads (dropped lineage)
    body += "".join(cfr(f"x{i}", 70001, 400, 60) + cfr(f"x{i}", 28901, 400, 60) for i in range(2))  # phage or Salmonella
    body += cfr("u", 0, 0, 0)
    f.write_text(body)
    return f


def table(path, col="s"):
    return pd.read_csv(path, sep="\t").set_index("feature_id")[col].to_dict()


def test_em_shares_ambiguous_reads_by_abundance_and_closes_the_ledger(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em")
    c = table(out / "counts.tsv")
    assert c["t|562|S"] == pytest.approx(18, abs=1e-4) and c["t|564|S"] == pytest.approx(2, abs=1e-4)
    # the two phage-or-Salmonella reads go to the phage: Salmonella has no read of its own (MLE)
    assert c.get("t|28901|S", 0) == pytest.approx(0, abs=1e-4)
    led = pd.read_csv(out / "ledger.tsv", sep="\t").iloc[0]
    assert led.dropped_lineages == pytest.approx(6, abs=1e-4)
    assert led.reads_shared_by_em == 12 and led.assignment == "em" and led.reads_assigned_by_lca == 0
    assert sum(led[k] for k in taxa.LEDGER_COLS[1:]) == pytest.approx(led.input_fragments, abs=1e-6)


def test_lca_is_still_the_default(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    taxa.build_table_from_per_read({"s": em_file(tmp_path)}, out, nodes, names, "S", ("Viruses",), "t")
    c = table(out / "counts.tsv")
    assert c == {"t|562|S": 9, "t|564|S": 1}          # the 10 shared reads go to the genus, above the rank
    led = pd.read_csv(out / "ledger.tsv", sep="\t").iloc[0]
    assert led.unresolved_above_rank == 12 and led.assignment == "lca"


def test_imputations_are_sealed_integer_tables_that_each_close_their_ledger(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    out = tmp_path / "taxa"
    f = em_file(tmp_path)
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=5, seed=3)
    rec = manifest.verify(out)
    assert rec["imputations"] == 5 and rec["assignment"] == "em" and rec["seed"] == 3
    rows = list(pd.read_csv(out / "counts.tsv", sep="\t")["feature_id"])
    for m in range(1, 6):
        assert f"imputations/counts_{m:03d}.tsv" in rec["artifacts"]
        ic = pd.read_csv(out / f"imputations/counts_{m:03d}.tsv", sep="\t")
        il = pd.read_csv(out / f"imputations/ledger_{m:03d}.tsv", sep="\t").iloc[0]
        assert list(ic["feature_id"]) == rows and ic["s"].dtype.kind == "i"
        assert ic["s"].sum() == il.kept_at_rank
        assert sum(il[k] for k in taxa.LEDGER_COLS[1:]) == il.input_fragments
        assert 9 <= ic.set_index("feature_id").loc["t|562|S", "s"] <= 19
    again = tmp_path / "again"
    taxa.build_table_from_per_read({"s": f}, again, nodes, names, "S", ("Viruses",), "t", 0, 0,
                                   assign="em", imputations=5, seed=3)
    assert (again / "imputations/counts_004.tsv").read_bytes() == (out / "imputations/counts_004.tsv").read_bytes()
    # a rerun without imputations removes the old ones, and the new seal lists none
    taxa.build_table_from_per_read({"s": f}, out, nodes, names, "S", ("Viruses",), "t", 0, 0, assign="em")
    assert not list((out / "imputations").glob("*.tsv"))
    assert manifest.verify(out)["imputations"] == 0


def test_imputations_need_em(tmp_path):
    nodes, names = taxonomy_files(tmp_path)
    with pytest.raises(ValueError, match="--imputations needs --assign em"):
        taxa.build_table_from_per_read({"s": em_file(tmp_path)}, tmp_path / "o", nodes, names, "S", (), "t",
                                       imputations=3)
