"""Unit tests for the fragment policy, feature binning and step caching (no external tools)."""
from collections import Counter

import pandas as pd
import pytest

from unclassified_meta import fragments as fr
from unclassified_meta import panel
from unclassified_meta.runner import Steps, params_key

S2F = {"contig|a1": "ns|A|like", "contig|a2": "ns|A|like", "ref|B": "ns|B|same",
       "decoy|host:h|chr1": "decoy_host:h", "decoy|rrna|r1": "decoy_rrna"}


def sam(name, flag, ref, score=0):
    tags = f"\tAS:i:{score}" if ref != "*" else ""
    return f"{name}\t{flag}\t{ref}\t1\t42\t100M\t*\t0\t0\tACGT\tIIII{tags}\n"


def outcome(*records):
    feats, ledger, n = fr.count_features([sam("f1", *r) for r in records], S2F)
    assert n == 1 and sum(feats.values()) + sum(ledger.values()) == 1
    return next(iter(feats or ledger))


def test_concordant_pair_counts_once_for_its_feature():
    assert outcome((99, "contig|a1"), (147, "contig|a1")) == "ns|A|like"


def test_only_mate_one_aligned_counts():
    assert outcome((73, "contig|a1"), (133, "*")) == "ns|A|like"


def test_only_mate_two_aligned_counts():
    assert outcome((69, "*"), (137, "ref|B")) == "ns|B|same"


def test_mates_in_different_features_are_ambiguous():
    assert outcome((65, "contig|a1"), (129, "ref|B")) == fr.AMBIGUOUS


def test_multimapping_within_one_feature_counts():
    assert outcome((99, "contig|a1"), (355, "contig|a2"), (147, "contig|a1"), (403, "contig|a2")) == "ns|A|like"


def test_equal_scores_in_two_features_are_ambiguous_but_a_worse_hit_is_ignored():
    assert outcome((0, "contig|a1", -5), (256, "ref|B", -5)) == fr.AMBIGUOUS
    assert outcome((0, "contig|a1", 0), (256, "ref|B", -40)) == "ns|A|like"


def test_target_and_decoy_tie_is_ambiguous_decoy_and_decoy_only_goes_to_its_class():
    assert outcome((0, "contig|a1", -6), (256, "decoy|host:h|chr1", -6)) == fr.AMBIGUOUS_DECOY
    assert outcome((0, "decoy|host:h|chr1")) == "decoy_host:h"
    assert outcome((0, "decoy|host:h|chr1", -3), (256, "decoy|rrna|r1", -3)) == fr.DECOY_MIXED
    assert outcome((0, "decoy|host:h|chr1", 0), (256, "contig|a1", -30)) == "decoy_host:h"


def test_unaligned_fragment_is_unmapped_and_ledger_sums_to_input():
    lines = [sam("p1", 77, "*"), sam("p1", 141, "*"),
             sam("p2", 99, "contig|a1"), sam("p2", 147, "contig|a1"),
             sam("p3", 65, "contig|a1"), sam("p3", 129, "ref|B"),
             sam("p4", 0, "decoy|rrna|r1")]
    feats, ledger, n = fr.count_features(lines, S2F)
    assert n == 4 and feats == Counter({"ns|A|like": 1})
    assert ledger == Counter({fr.UNMAPPED: 1, fr.AMBIGUOUS: 1, "decoy_rrna": 1})
    assert sum(feats.values()) + sum(ledger.values()) == n


def test_unknown_reference_sequence_is_an_error():
    with pytest.raises(KeyError):
        fr.count_features([sam("x", 0, "not_in_map")], S2F)


def test_primary_policy_counts_r2_only_fragments():
    lines = [sam("p1", 99, "c1"), sam("p1", 147, "c1"),          # concordant
             sam("p2", 73, "c1"), sam("p2", 133, "*"),           # mate 1 only
             sam("p3", 69, "*"), sam("p3", 137, "c2"),           # mate 2 only (lost by the old -f 64 filter)
             sam("p4", 65, "c1"), sam("p4", 129, "c2"),          # discordant: mate 1 decides
             sam("p5", 77, "*"), sam("p5", 141, "*"),            # unaligned
             sam("s1", 0, "c2"), sam("s1", 256, "c1")]           # single read with a secondary alignment
    counts, n = fr.count_primary(lines)
    assert n == 6 and counts == Counter({"c1": 3, "c2": 2})


def contig_table(rows):
    cols = ["contig", "length", "category", "layer", "method", "best_hit", "best_hit_title", "identity"]
    return pd.DataFrame(rows, columns=cols)


NS = {"virus": "refseq_viral_r1", "virus_protein": "refseq_viral_protein_r1"}


def test_two_proteins_of_one_virus_form_one_protein_only_feature():
    t = contig_table([
        ("c1", 900, "virus", "virus", "diamond_blastx", "YP_1.1", "YP_1.1 replicase [Virus Q]", 55.0),
        ("c2", 700, "virus", "virus", "diamond_blastx", "YP_2.1", "YP_2.1 coat protein [Virus Q]", 60.0),
        ("c3", 800, "virus", "virus", "diamond_blastx", "YP_9.1", "YP_9.1 coat protein [Virus R]", 70.0)])
    feats, c2f = panel.define_features(t, None, {"virus"}, NS)
    assert c2f["c1"] == c2f["c2"] == "refseq_viral_protein_r1|org:Virus_Q|protein_like"
    assert c2f["c3"] != c2f["c1"] and len(feats) == 2
    q = feats.set_index("feature_id").loc[c2f["c1"]]
    assert q["label"] == "Virus Q-like (protein only)" and q["n_contigs"] == 2 and "YP_1.1" in q["evidence"]


def test_identity_class_uses_read_weighted_identity_and_keeps_like_apart_from_same():
    t = contig_table([
        ("c1", 1000, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X, complete genome", 99.0),
        ("c2", 1000, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X, complete genome", 91.0),
        ("c3", 1000, "virus", "virus", "blastn", "NC_2.1", "NC_2.1 Virus Y, complete genome", 81.0),
        ("c4", 1000, "host:h", "host", "minimap2", "chr1", "chr1", 99.9)])
    feats, c2f = panel.define_features(t, {"c1": 900, "c2": 100, "c3": 50}, {"virus"}, NS)
    f = feats.set_index("feature_id")
    assert c2f["c1"] == c2f["c2"] == "refseq_viral_r1|NC_1.1|same"          # (99*900 + 91*100) / 1000 = 98.2
    assert f.loc[c2f["c1"], "label"] == "Virus X" and f.loc[c2f["c1"], "weighted_identity"] == 98.2
    assert bool(f.loc[c2f["c1"], "mixed_identity"]) is True                 # 10 % of the reads are below 95 %
    assert c2f["c3"] == "refseq_viral_r1|NC_2.1|like" and f.loc[c2f["c3"], "label"] == "Virus Y-like"
    assert "c4" not in c2f                                                   # host contigs are not selected
    # the same contigs weighted the other way round fall below the threshold
    _, c2f2 = panel.define_features(t, {"c1": 100, "c2": 900, "c3": 50}, {"virus"}, NS)
    assert c2f2["c1"] == "refseq_viral_r1|NC_1.1|like"


def test_unmatched_contigs_become_one_feature_per_sequence_cluster_never_one_pool():
    t = contig_table([("u1", 500, "unclassified", "-", "-", None, None, None),
                      ("u2", 400, "unclassified", "-", "-", None, None, None),
                      ("u3", 300, "unclassified", "-", "-", None, None, None)])
    clusters = panel.single_linkage(["u1", "u2", "u3"], [("u1", "u2")])
    assert clusters == {"u1": 1, "u2": 1, "u3": 2}
    feats, c2f = panel.define_features(t, None, {"unclassified"}, NS, clusters=clusters)
    assert c2f["u1"] == c2f["u2"] == "denovo|cluster_0001|unmatched" and c2f["u3"] == "denovo|cluster_0002|unmatched"
    assert set(feats["identity_class"]) == {"unmatched"} and len(feats) == 2


def test_step_marker_is_invalidated_when_its_parameters_change(tmp_path):
    steps = Steps(tmp_path)
    k_sub = params_key("readback", "asm1", 1_000_000, "--very-fast")
    assert not steps.done("readback_S1", k_sub)
    out = tmp_path / "S1.counts.tsv"
    out.write_text("contig\tcount\n")
    steps.mark("readback_S1", k_sub, [out])                 # a marker records the step's outputs
    assert steps.done("readback_S1", k_sub)
    assert not steps.done("readback_S1", params_key("readback", "asm1", None, "--very-fast"))   # --subsample all
    (tmp_path / ".done" / "old_step").write_text("")       # marker of a version without keys: not trusted
    assert not steps.done("old_step", "anything")


def test_step_keys_carry_the_cache_version(monkeypatch):
    from unclassified_meta import runner
    k = params_key("readback", "asm1")
    monkeypatch.setattr(runner, "CACHE_VERSION", "counting-policy-next")
    assert params_key("readback", "asm1") != k


def test_decoy_fasta_with_repeated_identifiers_gets_unique_panel_names(tmp_path):
    src = tmp_path / "rrna.fa"
    src.write_text(">r1 first\nACGT\n>r1 second\nGGCC\n>r2\nTTAA\n")
    seq2feature = {}
    with open(tmp_path / "panel.fa", "w") as fo:
        n, total = panel._append_fasta(src, fo, "decoy|rrna|", seq2feature, lambda _s: "decoy_rrna")
    assert (n, total) == (3, 12)
    assert list(seq2feature) == ["decoy|rrna|r1", "decoy|rrna|r1~2", "decoy|rrna|r2"]
    assert (tmp_path / "panel.fa").read_text().count(">") == 3


def test_reference_titles_are_cleaned_for_labels():
    assert panel.clean_title("NC_077021.1 UNVERIFIED: Cactus carlavirus 1 genomic sequence", "NC_077021.1") == "Cactus carlavirus 1"
    assert panel.clean_title("NC_040842.1 Potexvirus sp., complete genome", "NC_040842.1") == "Potexvirus sp."
    assert panel.clean_title("NC_077221.1 Prymnesium kappa virus genome assembly, chromosome: 1", "NC_077221.1") == "Prymnesium kappa virus"
    assert panel.clean_title("", "NC_1.1") == "NC_1.1"
