"""strain-evidence: pure parts (cs parsing, overlap comparison, clustering, linkage, blocks, error flags,
classification) on hand-made inputs, the panel-build split, and one small end-to-end run on synthetic data."""
import argparse
import gzip
import json
import random
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unclassified_meta import manifest, panel
from unclassified_meta import strains as st

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from make_synthetic import rand_seq, revcomp  # noqa: E402


def track(name, intervals, diffs=(), flags=(), ref="R", length=None):
    """A contig track from (pos, kind, ref, alt) differences."""
    t = st.Track(name, ref, length or sum(e - s for s, e in intervals), intervals=list(intervals),
                 flags=list(flags))
    for pos, kind, r, a in diffs:
        d = st.Diff(pos, kind, r, a)
        t.alleles[d.site] = d.allele
        t.ref_alleles[d.site] = d.ref if kind == "snv" else "."
        if kind == "del":
            t.deleted.append((pos, pos + len(r)))
    return t


# ---------------------------------------------------------------- cs / PAF

def test_cs_string_gives_reference_positions_of_snvs_and_indels():
    diffs, end = st.parse_cs(":10*ag:5+tt:3-ccc:2", 100)
    assert end == 100 + 10 + 1 + 5 + 3 + 3 + 2
    assert diffs == [st.Diff(110, "snv", "A", "G"), st.Diff(116, "ins", "", "TT"), st.Diff(119, "del", "CCC", "")]
    assert [d.allele for d in diffs] == ["G", "+TT", "-3"]
    assert diffs[1].site == (116, "indel")


def test_unreadable_cs_string_is_refused():
    with pytest.raises(ValueError, match="cs string"):
        st.parse_cs(":10?x", 0)


def test_paf_secondary_alignments_are_skipped_and_cs_end_is_checked():
    rows = ["c1\t500\t0\t100\t+\tR\t1000\t200\t300\t99\t100\t60\ttp:A:P\tcs:Z::50*ac:49",
            "c1\t500\t0\t100\t+\tR2\t1000\t200\t300\t99\t100\t0\ttp:A:S\tcs:Z::100"]
    segs = st.parse_paf(rows)
    assert len(segs) == 1 and segs[0].diffs == [st.Diff(250, "snv", "A", "C")]
    with pytest.raises(ValueError, match="ends at"):
        st.parse_paf(["c1\t500\t0\t100\t+\tR\t1000\t200\t301\t99\t100\t60\ttp:A:P\tcs:Z::100"])


# ---------------------------------------------------------------- flags

def seg(contig, qs, qe, rs, re_, strand="+", ref="R", diffs=(), qlen=1000):
    return st.Segment(contig, qlen, qs, qe, strand, ref, 5000, rs, re_, 60, list(diffs))


def test_two_collinear_pieces_are_a_split_alignment_and_reversed_pieces_are_non_collinear():
    t = st.build_track("c", 1000, "R", [seg("c", 0, 500, 100, 600), seg("c", 500, 1000, 900, 1400)])
    assert "split_alignment" in t.flags and t.intervals == [(100, 600), (900, 1400)]
    t = st.build_track("c", 1000, "R", [seg("c", 0, 500, 900, 1400), seg("c", 500, 1000, 100, 600)])
    assert "non_collinear" in t.flags


def test_a_contig_half_on_another_reference_is_multi_reference_and_low_aligned():
    t = st.build_track("c", 1000, "R", [seg("c", 0, 500, 0, 500), seg("c", 500, 1000, 0, 500, ref="Q")])
    assert {"multi_reference", "low_aligned_fraction"} <= set(t.flags) and t.other_refs == ["Q"]
    assert not t.linkable


def test_identity_jump_along_the_contig_is_flagged():
    diffs = [st.Diff(p, "snv", "A", "C") for p in range(0, 200, 10)]            # 90 % in the first window
    t = st.build_track("c", 1000, "R", [seg("c", 0, 1000, 0, 1000, diffs=diffs)])
    assert t.flags == ["identity_jump"]
    t = st.build_track("c", 1000, "R", [seg("c", 0, 1000, 0, 1000, diffs=diffs[:2])])
    assert t.flags == []


def test_unaligned_contig_is_flagged_not_aligned():
    assert st.build_track("c", 1000, "R", []).flags == ["not_aligned"]


# ---------------------------------------------------------------- overlap comparison

def test_differences_outside_the_overlap_do_not_count():
    a = track("a", [(0, 1000)], [(100, "snv", "A", "G"), (900, "snv", "C", "T")])
    b = track("b", [(500, 1500)], [(1200, "snv", "A", "G")])
    r = st.compare_tracks(a, b)
    assert r["overlap_bp"] == 500 and r["overlap"] == [(500, 1000)]
    assert r["n_snv"] == 1 and r["diffs"] == [((900, "snv"), "T", "C")]
    assert r["identity"] == pytest.approx(99.8)


def test_shared_differences_from_the_reference_are_not_differences_between_contigs():
    a = track("a", [(0, 1000)], [(300, "snv", "A", "G"), (400, "del", "AC", "")])
    b = track("b", [(0, 1000)], [(300, "snv", "A", "G"), (400, "del", "AC", "")])
    r = st.compare_tracks(a, b)
    assert r["n_snv"] == r["n_indel"] == 0 and r["identity"] == 100.0


def test_an_indel_counts_its_length_and_an_snv_inside_a_deletion_is_not_counted_again():
    a = track("a", [(0, 1000)], [(400, "del", "ACG", "")])
    b = track("b", [(0, 1000)], [(401, "snv", "C", "T"), (600, "ins", "", "GG")])
    r = st.compare_tracks(a, b)
    assert (r["n_snv"], r["n_indel"], r["indel_bp"]) == (0, 2, 5)
    assert r["identity"] == pytest.approx(99.5)


def test_disjoint_contigs_have_no_identity_and_relation_disjoint():
    r = st.compare_tracks(track("a", [(0, 500)]), track("b", [(600, 900)]))
    assert r["overlap_bp"] == 0 and r["identity"] is None
    assert st.relation(0, 0, 200) == "disjoint"
    assert st.relation(150, 3, 200) == "disjoint"
    assert st.relation(300, 0, 200) == "overlapping_identical"
    assert st.relation(300, 2, 200) == "overlapping_divergent"


# ---------------------------------------------------------------- clusters and classification

def pair(a, b, overlap, identity, pid="P1", diffs=()):
    return {"pair_id": pid, "contig_a": a, "contig_b": b, "overlap_bp": overlap, "identity": identity,
            "diffs": list(diffs)}


def test_clusters_link_only_long_enough_overlaps_above_the_threshold_and_never_flagged_contigs():
    names = ["a", "b", "c", "d", "x"]
    pairs = [pair("a", "b", 800, 99.5), pair("b", "c", 150, 100.0), pair("a", "d", 900, 97.5),
             pair("a", "x", 600, 100.0)]
    linkable = {"a", "b", "c", "d"}
    assert st.cluster_members(names, pairs, 97, 200, linkable) == [["a", "b", "d"], ["c"], ["x"]]
    assert st.cluster_members(names, pairs, 98, 200, linkable) == [["a", "b"], ["c"], ["d"], ["x"]]


def test_distinguishing_sites_come_from_divergent_overlaps_with_other_clusters_only():
    t = {n: track(n, [(0, 1000)]) for n in "abc"}
    t["c"].flags.append("multi_reference")
    pairs = [pair("a", "b", 900, 99.0, "P1", [((100, "snv"), "A", "G"), ((150, "indel"), ".", "-2")]),
             pair("a", "c", 900, 97.0, "P2", [((300, "snv"), "C", "T")])]
    dist, used = st.distinguishing_sites({"a"}, t, pairs, 200)
    assert dist == {(100, "snv"): "A", (150, "indel"): "."} and used == ["P1"]
    dist_b, _ = st.distinguishing_sites({"b"}, t, pairs, 200)
    assert dist_b == {(100, "snv"): "G", (150, "indel"): "-2"}


def test_consensus_identity_uses_the_majority_of_covering_members_and_cover_the_union():
    m = [track("a", [(0, 600)], [(100, "snv", "A", "G")]), track("b", [(400, 1000)], [(500, "snv", "A", "G")]),
         track("c", [(0, 1000)], [(500, "snv", "A", "G")])]
    ident, cover = st.consensus_stats(m, 2000)
    assert cover == 0.5
    assert ident == pytest.approx(100 * (1 - 1 / 1000))        # site 100: 1 G vs 1 ref -> tie keeps ref


def test_haplotype_carries_a_cluster_only_with_two_matching_sites_and_no_conflict():
    dist = {(10, "snv"): "A", (40, "snv"): "G", (90, "indel"): "-2"}
    sites = [(10, "snv"), (40, "snv"), (60, "snv")]
    assert st.haplotype_carries(dist, sites, ("A", "G", "T"))
    assert not st.haplotype_carries(dist, sites, ("A", "C", "T"))
    assert not st.haplotype_carries(dist, [(10, "snv"), (60, "snv")], ("A", "T"))


def test_classification_needs_the_evidence_of_each_level_and_flags_cap_it():
    assert st.classify(False, False, False, False) == "operational_sequence_cluster"
    assert st.classify(False, True, False, False) == "candidate_strain_lineage"
    assert st.classify(False, True, True, False) == "phased_viral_haplotype"
    assert st.classify(False, False, True, False) == "operational_sequence_cluster"
    assert st.classify(False, True, True, True) == "reference_supported_viral_strain"
    assert st.classify(True, True, True, True) == "operational_sequence_cluster"


# ---------------------------------------------------------------- statistics, linkage, blocks, haplotypes

def test_binomial_tail_matches_direct_summation():
    from math import comb
    n, p = 40, 0.05
    for k in (0, 1, 2, 5, 12, 40):
        want = sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))
        assert st.binom_sf(k, n, p) == pytest.approx(want, rel=1e-9, abs=1e-15)


def test_perfect_linkage_has_r2_and_d_prime_one_and_independence_zero():
    assert st.ld_stats(30, 0, 0, 70) == (pytest.approx(1.0), pytest.approx(1.0))
    r2, dp = st.ld_stats(25, 25, 25, 25)
    assert r2 == pytest.approx(0.0) and dp == pytest.approx(0.0)
    assert st.ld_stats(10, 5, 0, 0) == (None, None)               # site 1 monomorphic


def test_a_rare_pattern_one_step_from_a_frequent_one_is_a_possible_error():
    pats = {("A", "C"): 200, ("A", "T"): 1, ("G", "T"): 60, ("G", "C"): 30}
    err = st.possible_errors(pats, 0.005)
    assert err[("A", "T")] == ("A", "C")
    assert err[("G", "T")] is None and err[("A", "C")] is None
    assert err[("G", "C")] is None                                 # 30 reads are far beyond error of 200 or 60


def test_blocks_grow_while_enough_fragments_carry_every_site_and_single_sites_are_dropped():
    frags = [{0: "A", 1: "C"}] * 6 + [{1: "C", 2: "G"}] * 2 + [{3: "T", 4: "A"}] * 5
    assert st.find_blocks(5, frags, 5) == [[0, 1], [3, 4]]
    frags += [{0: "A", 1: "C", 2: "G"}] * 5
    assert st.find_blocks(5, frags, 5) == [[0, 1, 2], [3, 4]]


def test_block_haplotypes_count_complete_patterns_with_a_seeded_interval():
    frags = [{0: "A", 1: "C"}] * 70 + [{0: "G", 1: "T"}] * 30 + [{0: "A"}] * 10
    haps = st.block_haplotypes([0, 1], frags, 0.005, 3, np.random.default_rng(1))
    assert [(h["pattern"], h["count"]) for h in haps] == [(("A", "C"), 70), (("G", "T"), 30)]
    assert haps[0]["n"] == 100 and haps[0]["ci_low"] < 0.7 < haps[0]["ci_high"]
    again = st.block_haplotypes([0, 1], frags, 0.005, 3, np.random.default_rng(1))
    assert again[0]["ci_low"] == haps[0]["ci_low"]
    assert all(h["supported"] for h in haps)


def test_linkage_rows_tabulate_major_and_other_alleles():
    sites = [{"pos": 9, "kind": "snv", "major": "A"}, {"pos": 49, "kind": "snv", "major": "C"}]
    frags = [{0: "A", 1: "C"}] * 7 + [{0: "G", 1: "T"}] * 3
    (row,) = st.linkage_rows(sites, frags, 5)
    assert (row["site_a"], row["n_major_major"], row["n_other_other"], row["r2"]) == ("10:snv", 7, 3, 1.0)
    assert st.linkage_rows(sites, frags[:4], 5) == []


# ---------------------------------------------------------------- reads

def test_read_alleles_follow_the_cigar_and_respect_base_quality():
    seq = "ACGTAC" + "GG" + "TACG"            # 6M 2I 4M from 100, then 2D, then 4M
    seq += "TTTT"
    qual = "I" * 5 + "#" + "II" + "I" * 8
    cigar = "6M2I4M2D4M"
    sites = [100, 105, 106, 110, 112, 114]
    kinds = ["snv", "snv", "indel", "indel", "indel", "snv"]
    alleles = [{"A", "C"}, {"C", "T"}, {"+GG", "."}, {"-2", "."}, {"-2", "."}, {"T", "A"}]
    got = st.read_site_alleles(100, cigar, seq, qual, sites, kinds, alleles, 20)
    # 105 has quality '#' (2); 106 carries the insertion; 110 is the deletion; 112 has no aligned base before it
    assert got == {0: "A", 2: "+GG", 3: "-2", 5: "T"}


def test_read_without_an_indel_gets_the_reference_allele_at_an_indel_site():
    got = st.read_site_alleles(0, "20M", "A" * 20, "I" * 20, [10], ["indel"], [{".", "-3"}], 20)
    assert got == {0: "."}


def test_variant_calls_need_both_strands_and_pass_the_bonferroni_binomial_test():
    L = 50
    counts = np.zeros((L, 4, 2), dtype=np.int32)
    counts[:, 0, :] = 100                                      # reference A everywhere, 200 reads
    counts[10, 2, :] = [20, 20]                                # real minor G on both strands
    counts[20, 3, 0] = 40                                      # T on the forward strand only
    counts[30, 0, :] = 0
    counts[30, 1, :] = [100, 100]                              # fixed difference C
    sites, m = st.call_sites(counts, {}, "A" * L, 10, 5, 0.02, 0.005)
    by = {s["pos"]: s for s in sites}
    assert m == 3 * L
    assert by[10]["status"] == "polymorphic" and by[10]["alleles"] == ["A", "G"]
    assert 20 not in by
    assert by[30]["status"] == "fixed" and by[30]["major"] == "C"


def test_indel_sites_are_called_against_reads_without_the_indel():
    L = 30
    counts = np.zeros((L, 4, 2), dtype=np.int32)
    counts[:, 0, :] = 50
    sites, _ = st.call_sites(counts, {(12, "-2"): [15, 15]}, "A" * L, 10, 5, 0.02, 0.005)
    (s,) = [x for x in sites if x["kind"] == "indel"]
    assert s["pos"] == 12 and s["alleles"] == [".", "-2"] and s["depth"] == 130


# ---------------------------------------------------------------- panel-build split

NS = {"virus": "refseq_viral_r1"}


def contig_table(rows):
    cols = ["contig", "length", "category", "layer", "method", "best_hit", "best_hit_title", "identity"]
    return pd.DataFrame(rows, columns=cols)


def test_strain_clusters_split_a_reference_bin_and_keep_the_parent_feature():
    t = contig_table([
        ("c1", 1000, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X, complete genome", 99.0),
        ("c2", 900, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X, complete genome", 99.0),
        ("c3", 800, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X, complete genome", 90.0),
        ("c4", 700, "virus", "virus", "blastn", "NC_2.1", "NC_2.1 Virus Y, complete genome", 99.0)])
    feats, c2f = panel.define_features(t, None, {"virus"}, NS, strain_clusters={"c1": "c98_001", "c2": "c98_001",
                                                                                 "c3": "c98_002"})
    f = feats.set_index("feature_id")
    assert c2f["c1"] == c2f["c2"] == "refseq_viral_r1|NC_1.1|same|c98_001"
    assert c2f["c3"] == "refseq_viral_r1|NC_1.1|like|c98_002"
    assert f.loc[c2f["c3"], "label"] == "Virus X-like [c98_002]"
    # the parent is what the whole bin would be (read-weighted identity over c1..c3, here by length: 96.4 %)
    assert set(f.loc[[c2f["c1"], c2f["c3"]], "parent_feature"]) == {"refseq_viral_r1|NC_1.1|same"}
    assert c2f["c4"] == "refseq_viral_r1|NC_2.1|same" and f.loc[c2f["c4"], "parent_feature"] == ""
    assert list(feats.columns) == panel.FEATURE_COLUMNS


def test_without_strain_clusters_features_are_unchanged():
    t = contig_table([("c1", 1000, "virus", "virus", "blastn", "NC_1.1", "NC_1.1 Virus X", 99.0)])
    feats, c2f = panel.define_features(t, None, {"virus"}, NS)
    assert c2f == {"c1": "refseq_viral_r1|NC_1.1|same"} and feats.loc[0, "parent_feature"] == ""


def test_strain_cluster_table_must_be_sealed_and_have_the_resolution(tmp_path):
    f = tmp_path / "viral_strain_clusters.tsv"
    f.write_text("# rule\nresolution\tcluster_id\treference\tmembers\n98\tc98_001\tR\ta,b\n99\tc99_001\tR\ta\n")
    with pytest.raises(manifest.ManifestError):
        st.read_strain_clusters(f, 98)
    manifest.seal(tmp_path, ["viral_strain_clusters.tsv"], {"step": "test"})
    m, _ = st.read_strain_clusters(f, 98)
    assert m == {"a": "c98_001", "b": "c98_001"}
    with pytest.raises(ValueError, match="resolutions: 98, 99"):
        st.read_strain_clusters(f, 97)


# ---------------------------------------------------------------- end to end

TOOLS = ("minimap2", "bowtie2", "bowtie2-build")


def _mutate_at(seq, positions, rng):
    s = list(seq)
    for p in positions:
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def _write_fastq(path, reads):
    with gzip.open(path, "wt") as fo:
        for name, s in reads:
            fo.write(f"@{name}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.fixture(scope="module")
def e2e(tmp_path_factory):
    if any(shutil.which(t) is None for t in TOOLS):
        pytest.skip("minimap2 / bowtie2 not on PATH")
    rng = random.Random(11)
    d = tmp_path_factory.mktemp("strains")
    ref1, ref2 = rand_seq(rng, 3200, 0.5), rand_seq(rng, 2000, 0.5)
    strain_a = _mutate_at(ref1, rng.sample(range(3200), 16), rng)             # 99.5 % to the reference
    b_sites = list(range(620, 1900, 40))                                        # 32 sites, 40 bp apart
    strain_b = _mutate_at(strain_a, b_sites, rng)
    contigs = {"A1": strain_a[0:1500], "A2": strain_a[1200:2600], "B1": strain_b[500:2000],
               "D1": strain_a[2700:3200], "CH": strain_a[0:600] + ref2[1200:1800], "E1": ref2[0:1000]}
    (d / "run" / "contigs").mkdir(parents=True)
    (d / "run" / "tables").mkdir()
    with open(d / "run" / "contigs" / "contigs.fa", "w") as fo:
        for n, s in contigs.items():
            fo.write(f">{n}\n{s}\n")
    hit = {n: "REF1" for n in contigs}
    hit["E1"] = "REF2"
    rows = [(n, len(s), "virus", "virus", "blastn", hit[n], f"{hit[n]} test virus", 99.0, 1.0)
            for n, s in contigs.items()]
    rows.append(("P1", 900, "virus", "virus", "diamond_blastx", "YP_1.1", "YP_1.1 rep [Virus Q]", 50.0, 1.0))
    pd.DataFrame(rows, columns=["contig", "length", "category", "layer", "method", "best_hit", "best_hit_title",
                                "identity", "cover"]).to_csv(d / "run" / "tables" / "contig_category.tsv", sep="\t",
                                                             index=False)
    with open(d / "refs.fa", "w") as fo:
        fo.write(f">REF1 test virus 1\n{ref1}\n>REF2 test virus 2\n{ref2}\n>REF3 unused\n{rand_seq(rng, 500)}\n")
    (d / "run" / "references.used.json").write_text(json.dumps({"references": {}}))
    (d / "run" / "config.used.yaml").write_text("references: {}\n")
    manifest.seal(d / "run", ["contigs/contigs.fa", "tables/contig_category.tsv", "references.used.json",
                              "config.used.yaml"], {"step": "test fixture"})
    # read pairs from two haplotypes, 70 % strain A and 30 % strain B
    r1, r2 = [], []
    for k in range(2400):
        g = strain_a if k % 10 < 7 else strain_b
        L = rng.randint(250, 350)
        s0 = rng.randint(0, len(g) - L)
        f = g[s0:s0 + L]
        if rng.random() < 0.5:
            f = revcomp(f)
        r1.append((f"f{k}/1", f[:100]))
        r2.append((f"f{k}/2", revcomp(f)[:100]))
    _write_fastq(d / "s1_R1.fastq.gz", r1)
    _write_fastq(d / "s1_R2.fastq.gz", r2)
    (d / "samples.tsv").write_text(f"sample\tr1\tr2\ns1\t{d / 's1_R1.fastq.gz'}\t{d / 's1_R2.fastq.gz'}\n")
    p = argparse.ArgumentParser()
    st.add_arguments(p)
    args = p.parse_args(["--run", str(d / "run"), "-o", str(d / "out"), "--reference", str(d / "refs.fa"),
                         "--samples", str(d / "samples.tsv"), "-t", "2"])
    assert st.main(args) == 0
    return d, b_sites


def test_end_to_end_clusters_flags_and_classes(e2e):
    d, _ = e2e
    out = d / "out"
    manifest.verify(out, required=st.OUTPUT_TABLES)
    stat = st.read_commented_tsv(out / "contig_status.tsv").set_index("contig")
    assert stat.loc["P1", "status"].startswith("not_analysable")
    assert {"multi_reference", "low_aligned_fraction"} <= set(stat.loc["CH", "flags"].split(","))
    assert stat.loc["A1", "flags"] == ""
    ov = st.read_commented_tsv(out / "viral_cluster_overlap.tsv")
    ab = ov[(ov.contig_a == "A1") & (ov.contig_b == "B1")].iloc[0]
    assert ab.relation == "overlapping_divergent" and int(ab.n_snv) == 22          # sites 620..1460 inside 500..1500
    assert ov[(ov.contig_a == "A1") & (ov.contig_b == "A2")].iloc[0].relation == "overlapping_identical"
    cl = st.read_commented_tsv(out / "viral_strain_clusters.tsv")
    c98 = cl[cl.resolution == "98"].set_index("members")
    assert c98.loc["A1,A2", "class"] == "reference_supported_viral_strain"
    assert c98.loc["B1", "class"] == "phased_viral_haplotype" and c98.loc["B1", "haplotype_blocks"]
    assert "unresolved_disjoint" in c98.loc["D1", "cluster_flags"]
    assert c98.loc["CH", "class"] == "operational_sequence_cluster" and "CH:" in c98.loc["CH", "chimera_flags"]
    c97 = cl[cl.resolution == "97"]
    assert "A1,A2,B1" in set(c97.members) or "A1,B1,A2" in set(c97.members)
    assert list(out.glob("figures/overlap_tracks_REF1.png")) and list(out.glob("figures/linkage_REF1_*.png"))


def test_end_to_end_local_haplotype_frequency_interval_covers_the_true_share(e2e):
    d, b_sites = e2e
    hap = st.read_commented_tsv(d / "out" / "viral_haplotype_evidence.tsv")
    sites = st.read_commented_tsv(d / "out" / "variant_sites.tsv")
    poly = sites[sites.status == "polymorphic"]
    assert {int(p) - 1 for p in poly.position} <= set(b_sites)
    assert len(poly) >= 25
    top = hap[hap.supported == "True"].copy()
    top["n"] = top.n_fragments.astype(int)
    blk = top[top.n == top.n.max()]
    major = blk[blk["count"].astype(int) == blk["count"].astype(int).max()].iloc[0]
    assert int(major.n_sites) >= 2
    assert float(major.ci_low) <= 0.7 <= float(major.ci_high)
    lim = st.read_commented_tsv(d / "out" / "phasing_limits.tsv").set_index("reference")
    assert int(lim.loc["REF1", "n_blocks"]) >= 1 and float(lim.loc["REF1", "insert_median"]) > 250


def test_panel_build_splits_the_reference_bin_by_the_strain_clusters(e2e, tmp_path):
    if shutil.which("seqkit") is None:
        pytest.skip("seqkit not on PATH")
    d, _ = e2e
    out = panel.build_panel(d / "run", tmp_path / "panel", ["virus"], [], None, str(d / "refs.fa"), None,
                            other_contig_decoy=False, threads=2, use_run_hosts=False, use_run_rrna=False,
                            strain_clusters=str(d / "out" / "viral_strain_clusters.tsv"), strain_resolution=98)
    f = pd.read_csv(out / "features.tsv", sep="\t", keep_default_na=False).set_index("feature_id")
    kids = f[(f.parent_feature != "") & (f.accession == "REF1") & (f.identity_class != "reference")]
    assert len(kids) == 4 and all(i.split("|")[-1].startswith("c98_") for i in kids.index)
    ref_rows = f[(f.accession == "REF1") & (f.identity_class == "reference")]
    assert len(ref_rows) == 1 and ref_rows.iloc[0].parent_feature.endswith("|REF1|same")
    assert ref_rows.iloc[0].evidence == "reference genome of a bin split into operational clusters"
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["inputs"]["strain_clusters"]["sha256"] == manifest.sha256(d / "out" / "viral_strain_clusters.tsv")
