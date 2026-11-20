"""Panel workflow for differential abundance of viruses / selected taxa.

panel-build   from a finished `run`: selected contigs + the reference genomes they hit (+ extra targets) become
              *features*; host genome(s), rRNA and the other contigs of the run become *decoys*.
panel-quant   all reads of every sample (original FASTQ, full depth) are aligned to the panel with Bowtie 2 and
              counted once per fragment under the 'feature' policy of fragments.py.

Feature = reference-anchored bin. Contigs are binned by their closest reference; the identity class is part of
the feature id, so a divergent "X-like" bin is never merged with X itself:
  same          BLASTN, read-weighted nucleotide identity >= --same-identity (default 95 %)   label "X"
  like          BLASTN/TBLASTX below that                                                      label "X-like"
  protein_like  DIAMOND blastx only (binned by the organism of the protein)                   label "X-like (protein only)"
  kraken2       Kraken2 taxid of the contig
  unmatched     no reference match: sequence clusters (single linkage on BLASTN), one feature per cluster
  reference     a reference genome kept as its own feature next to a "like" bin, or an extra target sequence
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import yaml

from . import __version__
from . import manifest
from .config import Sample, check_sample_names
from .fragments import (AMBIGUOUS, AMBIGUOUS_DECOY, DECOY_MIXED, DECOY_PREFIX, UNMAPPED, count_classes, count_features,
                        outcome_of_class,
                        run_bowtie2)
from .io_utils import make_staging_dir, open_text, read_fasta, refuse_overlap, remove_staging_dir, staged
from .runner import StepError, file_id, params_key, require, sh

log = logging.getLogger("unclassified_meta")

FEATURE_COLUMNS = ["feature_id", "label", "namespace", "accession", "identity_class", "n_sequences",
                   "total_length", "evidence", "n_contigs", "weighted_identity", "identity_min", "identity_max",
                   "mixed_identity", "reference_included", "best_hit_title", "parent_feature"]


# ---------------------------------------------------------------- feature definition (pure, unit-tested)

def organism_of(title: str) -> str | None:
    """Organism in the trailing [brackets] of a RefSeq protein title, if any."""
    m = re.search(r"\[([^\[\]]+)\]\s*$", title or "")
    return m.group(1).strip() if m else None


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", text).strip("_")


def clean_title(title: str, accession: str) -> str:
    """Reference title without its accession and trailing ', complete genome' style suffixes."""
    t = (title or "").strip()
    if t.startswith(accession):
        t = t[len(accession):].strip()
    t = re.sub(r"^(MAG|UNVERIFIED|TPA|TPA_asm):\s*", "", t)
    t = re.sub(r",?\s*((complete|partial) (genome|sequence|cds)|genome assembly|genomic sequence).*$", "", t)
    return t.strip(" ,") or accession


def weighted_mean(values: list[float], weights: list[float]) -> float | None:
    w = sum(weights)
    return sum(v * x for v, x in zip(values, weights)) / w if w > 0 else (sum(values) / len(values) if values else None)


def _bin_row(ns: str, acc: str, kind: str, members: list, weights: dict[str, float] | None,
             same_identity: float) -> dict:
    """One feature row from the contigs of a bin (members: (record, method, hit, title))."""
    recs = [m[0] for m in members]
    # identity of a nucleotide bin = nucleotide identity of its BLASTN hits only; TBLASTX identities are
    # amino-acid identities and are not mixed in (a translated hit can never make a bin 'same')
    idrec = [m[0] for m in members if pd.notna(m[0].identity) and str(m[0].identity) != ""
             and (kind != "nt" or "tblastx" not in m[1])]
    ids = [float(x.identity) for x in idrec]
    wts = [float((weights or {}).get(x.contig, 0) or 0) for x in idrec]
    if not weights or sum(wts) == 0:
        wts = [float(x.length) for x in idrec]
    wid = weighted_mean(ids, wts)
    methods = sorted({m[1] for m in members})
    title = Counter(m[3] for m in members).most_common(1)[0][0]
    mixed = False
    if kind == "nt":
        cls = "same" if (wid is not None and wid >= same_identity) else "like"
        total = sum(wts) or 1.0
        hi = sum(w for i, w in zip(ids, wts) if i >= same_identity) / total
        mixed = 0.1 <= hi <= 0.9          # contigs on both sides of the threshold: possibly two entities
        name = clean_title(title, acc)
        label = name if cls == "same" else f"{name}-like"
    elif kind == "protein_like":
        cls = "protein_like"
        label = f"{organism_of(title) or clean_title(title, acc)}-like (protein only)"
    elif kind == "kraken2":
        cls, label = "kraken2", f"{title or acc} (Kraken2 contig assignment)"
    else:
        cls, label = "unmatched", f"unmatched sequence {acc.replace('_', ' ')}"
    hits = sorted({m[2] for m in members if m[2]})
    return {"feature_id": f"{ns}|{acc}|{cls}", "label": label, "namespace": ns, "accession": acc,
            "identity_class": cls, "n_sequences": len(recs), "total_length": int(sum(int(x.length) for x in recs)),
            "evidence": ";".join(methods) + (f" ({','.join(hits[:5])}{',...' if len(hits) > 5 else ''})"
                                             if kind == "protein_like" else ""),
            "n_contigs": len(recs),
            # nucleotide identity only; amino-acid identities of protein-only bins are not comparable
            "weighted_identity": round(wid, 2) if (wid is not None and kind == "nt") else "",
            "identity_min": round(min(ids), 2) if (ids and kind == "nt") else "",
            "identity_max": round(max(ids), 2) if (ids and kind == "nt") else "",
            "mixed_identity": mixed, "reference_included": False, "best_hit_title": title, "parent_feature": ""}


def define_features(contigs: pd.DataFrame, weights: dict[str, float] | None, categories: set[str],
                    namespaces: dict[str, str], same_identity: float = 95.0,
                    clusters: dict[str, int] | None = None,
                    only_contigs: set[str] | None = None,
                    strain_clusters: dict[str, str] | None = None) -> tuple[pd.DataFrame, dict[str, str]]:
    """Bin the selected contigs into features.

    contigs          contig_category.tsv (contig, length, category, layer, method, best_hit, best_hit_title, identity)
    weights          pooled read-back fragments per contig (identity is weighted by them; contig length when absent)
    namespaces       {'virus': ..., 'virus_protein': ..., 'extra_virus': ...} versioned namespace strings
    clusters         contig -> cluster number for contigs without a reference match
    strain_clusters  contig -> operational cluster id (strain-evidence, e.g. c98_003): the nucleotide bin of a
                     reference is split, one feature per cluster (feature_id ns|acc|class|c98_003, label suffixed
                     [c98_003], parent_feature = the feature id the whole bin would have had). Contigs of the bin
                     without a cluster stay in an ordinary bin of the remaining contigs.
    Returns (features table, contig -> feature_id).
    """
    if only_contigs is not None:
        sel = contigs[contigs["contig"].isin(only_contigs)]
    else:
        sel = contigs[contigs["category"].isin(categories)]
    bins: dict[tuple, list] = defaultdict(list)
    for r in sel.itertuples(index=False):
        method = r.method if isinstance(r.method, str) else "-"
        hit = r.best_hit if isinstance(r.best_hit, str) and r.best_hit not in ("", "nan") else ""
        title = r.best_hit_title if isinstance(r.best_hit_title, str) else ""
        if not hit:
            cl = (clusters or {}).get(r.contig)
            key = ("denovo", f"cluster_{cl:04d}" if cl is not None else f"contig_{r.contig}", "unmatched")
        elif method == "diamond_blastx":
            org = organism_of(title)
            ns = namespaces.get("virus_protein" if r.category == "virus" else "protein", "protein_ref")
            key = (ns, f"org:{slug(org)}" if org else hit, "protein_like")
        elif method == "kraken2":
            key = (namespaces.get("kraken2", "kraken2"), f"taxid:{hit}", "kraken2")
        else:
            layer = r.layer if isinstance(r.layer, str) else "ref"
            if method in ("blastn_extra", "tblastx_extra"):
                ns = namespaces.get("extra_virus", "extra_virus")
            elif layer == "virus":
                ns = namespaces.get("virus", "virus_ref")
            else:
                ns = f"{layer}_ref"
            key = (ns, hit, "nt")
        bins[key].append((r, method, hit, title))

    rows, c2f = [], {}
    for (ns, acc, kind), members in bins.items():
        groups: dict[str, list] = defaultdict(list)
        for m in members:
            groups[(strain_clusters or {}).get(m[0].contig, "") if kind == "nt" else ""].append(m)
        split = len(groups) > 1 or "" not in groups
        parent = _bin_row(ns, acc, kind, members, weights, same_identity)["feature_id"] if split else ""
        for cid in sorted(groups):
            row = _bin_row(ns, acc, kind, groups[cid], weights, same_identity)
            if cid:
                row["feature_id"] = f"{row['feature_id']}|{cid}"
                row["label"] = f"{row['label']} [{cid}]"
                row["parent_feature"] = parent
            rows.append(row)
            for m in groups[cid]:
                c2f[m[0].contig] = row["feature_id"]
    feats = pd.DataFrame(rows, columns=FEATURE_COLUMNS)
    if len(feats):
        feats = feats.sort_values(["total_length", "feature_id"], ascending=[False, True]).reset_index(drop=True)
    return feats, c2f


def single_linkage(names: list[str], pairs: list[tuple[str, str]]) -> dict[str, int]:
    """Union-find clusters; cluster numbers start at 1 in the order of `names`."""
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    number, out = {}, {}
    for n in names:
        root = find(n)
        number.setdefault(root, len(number) + 1)
        out[n] = number[root]
    return out


# ---------------------------------------------------------------- panel-build

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True).stdout.strip().splitlines()
        return out[0] if out else "unknown"
    except OSError:
        return "not found"


def _version_rows(text: str) -> dict[str, tuple[str, str]]:
    """reference -> (version, date) from the REFS_VERSION.tsv text."""
    lines = [l.split("\t") for l in (text or "").splitlines() if l.strip()]
    if not lines or "version" not in lines[0] or not ({"reference", "name"} & set(lines[0])):
        return {}
    h = lines[0]
    ri, vi = h.index("reference" if "reference" in h else "name"), h.index("version")
    di = h.index("date") if "date" in h else None
    return {r[ri]: (r[vi] if len(r) > vi else "", r[di] if di is not None and len(r) > di else "")
            for r in lines[1:] if len(r) > max(ri, vi)}


def refs_namespaces(refs_dir: str | None, extra: str | None, text: str | None = None) -> tuple[dict[str, str], str]:
    """Versioned namespaces, each from the version of its own reference in REFS_VERSION.tsv (read from refs_dir
    when `text` is not given): RefSeq viral genomes and proteins by their release, Swiss-Prot by its fetch date.
    A reference whose version is not recorded, or recorded as unknown, gets 'unversioned'."""
    if text is None:
        vt = Path(refs_dir) / "REFS_VERSION.tsv" if refs_dir else None
        text = vt.read_text() if vt and vt.exists() else ""
    rows = _version_rows(text)

    def release(key: str) -> str:
        v = rows.get(key, ("", ""))[0] if rows else (text or "")       # no table: a free-text version note
        m = re.search(r"RefSeq release (\d+)", v)
        return f"r{m.group(1)}" if m else "unversioned"

    sp_ver, sp_date = rows.get("sprot_viridiplantae", ("", ""))
    sp = re.search(r"fetched (\d{4}-\d{2}-\d{2})", sp_ver)
    sp_tag = f"d{sp.group(1)}" if sp else (f"d{sp_date}" if re.fullmatch(r"\d{4}-\d{2}-\d{2}", sp_date or "") else "unversioned")
    ns = {"virus": f"refseq_viral_{release('refseq_viral_genomic')}",
          "virus_protein": f"refseq_viral_protein_{release('refseq_viral_protein')}",
          "protein": f"protein_ref_{sp_tag}", "kraken2": "kraken2_taxid"}
    if extra:
        ns["extra_virus"] = f"extra_virus:{slug(Path(extra).name)}"
    return ns, text


def _has_sequence(fasta: str | Path) -> bool:
    """At least one record with a non-empty sequence (reads only until the first one)."""
    return any(seq for _, seq in read_fasta(fasta))


def cluster_unmatched(contig_fa: Path, names: list[str], work: Path, threads: int,
                      min_identity: float = 90.0, min_cover: float = 0.5) -> dict[str, int]:
    """Single-linkage clusters of unmatched contigs: BLASTN >= min_identity over >= min_cover of the shorter one."""
    if not names:
        return {}
    require("blastn", "makeblastdb")
    seqs = {n: s for n, s in read_fasta(contig_fa) if n in set(names)}
    fa = work / "unmatched.fa"
    with open(fa, "w") as fo:
        for n in names:
            fo.write(f">{n}\n{seqs[n]}\n")
    sh(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(work / "unmatched_db")])
    out = work / "unmatched_self.tsv"
    sh(["blastn", "-task", "megablast", "-query", str(fa), "-db", str(work / "unmatched_db"),
        "-perc_identity", str(min_identity), "-evalue", "1e-10", "-num_threads", str(threads),
        "-outfmt", "6 qseqid sseqid pident length", "-out", str(out)])
    pairs = []
    for line in open(out):
        q, s, _pid, alen = line.rstrip("\n").split("\t")
        if q != s and int(alen) >= min_cover * min(len(seqs[q]), len(seqs[s])):
            pairs.append((q, s))
    order = sorted(names, key=lambda n: (-len(seqs[n]), n))
    return single_linkage(order, pairs)


def _append_fasta(src: Path, fo, prefix: str, seq2feature: dict, value_of) -> tuple[int, int]:
    """Copy `src` into the open panel FASTA with renamed headers; record sequence -> feature/decoy. A repeated
    identifier gets a ~2, ~3 ... suffix so that panel sequence names stay unique."""
    n = total = 0
    name = None
    with open_text(src) as fh:
        for line in fh:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                name, k = f"{prefix}{sid}", 1
                while name in seq2feature:        # reference sets can repeat an identifier (e.g. SILVA SSU/LSU)
                    k += 1
                    name = f"{prefix}{sid}~{k}"
                seq2feature[name] = value_of(sid)
                fo.write(f">{name}\n")
                n += 1
            elif line.strip():
                fo.write(line if line.endswith("\n") else line + "\n")
                total += len(line.strip())
    return n, total


def check_references_against_run(recorded: dict, virus, rrna, extra_virus, hosts) -> None:
    """Every reference the panel consumes (the run's own and any override) must have the content the run
    recorded for it in references.used.json; a reference the run did not use has no record to match."""
    from .classify import source_digest
    rec_hosts = {r[0]: r[1] for r in (recorded.get("host") or []) if isinstance(r, list) and len(r) == 2}
    checks = [("virus", virus, recorded.get("virus")), ("rrna", rrna, recorded.get("rrna")),
              ("extra_virus_fasta", extra_virus, recorded.get("extra_virus_fasta"))]
    checks += [(f"host {h['name']}", h["fasta"], rec_hosts.get(h["name"])) for h in hosts]
    bad = []
    for label, path, rec in checks:
        if not path or not Path(path).is_file() or not isinstance(rec, list) or len(rec) != 2 \
                or not isinstance(rec[1], str):
            continue                                     # a missing file is reported as missing further on
        if source_digest(path) != rec[1]:
            bad.append(f"{label} ({path})")
    if bad:
        raise ValueError(f"reference content differs from what the run recorded in references.used.json: "
                         f"{'; '.join(bad)}. Build the panel with the references the run classified with, or run "
                         "`run` again with the new references")


def build_panel(run_dir: Path, outdir: Path, categories: list[str], hosts: list[dict], rrna: str | None,
                virus: str | None, refs_dir: str | None, extra_targets: str | None = None,
                contig_list: str | None = None, same_identity: float = 95.0, add_reference_genomes: bool = True,
                other_contig_decoy: bool = True, threads: int = 4, use_run_hosts: bool = True,
                use_run_rrna: bool = True, strain_clusters: str | None = None,
                strain_resolution: float = 98.0) -> Path:
    require("bowtie2-build", "seqkit")
    run_dir, outdir = Path(run_dir), Path(outdir)
    ctab, cfa = run_dir / "tables" / "contig_category.tsv", run_dir / "contigs" / "contigs.fa"
    for f in (ctab, cfa):
        if not f.exists():
            raise FileNotFoundError(f"{f} not found; --run must be the output folder of a finished `unclassified-meta run`")
    # the run must be a complete generation: contigs, their classification and counts sealed together
    manifest.verify(run_dir, required=("contigs/contigs.fa", "tables/contig_category.tsv", "references.used.json"),
                    what="run outputs")
    used_cfg = {}
    if (run_dir / "config.used.yaml").exists():
        used_cfg = yaml.safe_load((run_dir / "config.used.yaml").read_text()) or {}
    uref = used_cfg.get("references", {})
    hosts = hosts or (uref.get("host") if use_run_hosts else None) or []
    rrna = rrna or (uref.get("rrna") if use_run_rrna else None)
    virus = virus or uref.get("virus")
    # namespaces come from the references the run classified with (references.used.json, written by `run`),
    # not from whatever the refs folder holds today
    used_refs = run_dir / "references.used.json"
    recorded = json.loads(used_refs.read_text()) if used_refs.exists() else None
    if recorded and refs_dir and recorded.get("refs_dir") and \
            Path(refs_dir).resolve() != Path(recorded["refs_dir"]).resolve():
        raise ValueError(f"--refs {refs_dir} differs from the references the run used ({recorded['refs_dir']}); "
                         "build the panel with the run's references (leave --refs out) or rerun `run`")
    refs_dir = refs_dir or uref.get("refs_dir")
    extra_virus = uref.get("extra_virus_fasta")
    # a configured reference that is missing is an error, never a silent omission; leave one out explicitly
    # with --no-host-decoy, --no-rrna-decoy or --no-reference-genomes
    for h in hosts:
        if not Path(h["fasta"]).exists():
            raise FileNotFoundError(f"host genome {h['name']}: {h['fasta']} not found (pass --host NAME=FASTA, or "
                                    "--no-host-decoy to build the panel without host decoys)")
    if rrna and not Path(rrna).exists():
        raise FileNotFoundError(f"rRNA reference {rrna} not found (pass --rrna FASTA, or --no-rrna-decoy)")
    if recorded is not None:
        check_references_against_run(recorded.get("references") or {}, virus, rrna, extra_virus, hosts)
    # a decoy that loads nothing must not silently remove the run's host/rRNA contigs from the decoys
    for label, f in [(f"host genome {h['name']}", h["fasta"]) for h in hosts] + ([("rRNA reference", rrna)] if rrna else []):
        if not _has_sequence(f):
            raise ValueError(f"{label} {f} has no sequence; give a non-empty FASTA or leave this decoy out "
                             "(--no-host-decoy / --no-rrna-decoy)")
    refuse_overlap([outdir / n for n in PANEL_ARTIFACTS + ("provenance.json",)],
                   [ctab, cfa, run_dir / "tables" / "contig_counts.tsv", virus, rrna, extra_targets, uref.get("extra_virus_fasta"),
                    strain_clusters] + [h["fasta"] for h in hosts], "panel-build")
    strain_map = None
    if strain_clusters:
        # operational clusters of strain-evidence (sealed; made from this very run's contig table)
        from .strains import read_strain_clusters
        strain_map, sc_rec = read_strain_clusters(strain_clusters, strain_resolution)
        made_from = ((sc_rec.get("inputs") or {}).get("contig_category.tsv") or {}).get("sha256")
        if made_from != sha256(ctab):
            raise ValueError(f"--strain-clusters {strain_clusters} was made from another contig table than "
                             f"{ctab}; run strain-evidence on this run")
    outdir.mkdir(parents=True, exist_ok=True)
    # a rebuild first invalidates the old panel: without provenance.json panel-quant refuses the folder, so an
    # interrupted rebuild can never be mistaken for the previous panel
    manifest.begin(outdir, "panel-build", name="provenance.json")
    work = outdir / "work"
    work.mkdir(exist_ok=True)
    stage = make_staging_dir(outdir, "panel-build")

    contigs = pd.read_csv(ctab, sep="\t", dtype={"best_hit": str, "best_hit_title": str, "method": str, "layer": str})
    weights = None
    cc = run_dir / "tables" / "contig_counts.tsv"
    if cc.exists():
        m = pd.read_csv(cc, sep="\t").set_index("contig")
        weights = m.sum(axis=1).to_dict()
    else:
        log.warning("%s not found (run made by an older version): identity is weighted by contig length", cc)
    only = None
    if contig_list:
        only = {l.strip().split()[0] for l in open(contig_list) if l.strip() and not l.startswith("#")}
        missing = only - set(contigs["contig"])
        if missing:
            raise ValueError(f"--contigs: {len(missing)} names are not in {ctab}, e.g. {sorted(missing)[:3]}")
    cats = set(categories)
    picked = contigs[contigs["contig"].isin(only)] if only is not None else contigs[contigs["category"].isin(cats)]
    if picked.empty and not extra_targets:
        raise ValueError(f"no contigs in categories {sorted(cats)}; categories present: "
                         f"{sorted(contigs['category'].unique())}")
    unmatched = [r.contig for r in picked.itertuples() if not (isinstance(r.best_hit, str) and r.best_hit not in ("", "nan"))]
    clusters = cluster_unmatched(cfa, unmatched, work, threads)
    if recorded is not None:
        namespaces, refs_version = refs_namespaces(None, extra_virus, recorded.get("refs_version") or "")
    else:
        log.warning("%s not found (run made by an older version): namespaces are taken from the current %s",
                    used_refs, refs_dir or "references")
        namespaces, refs_version = refs_namespaces(refs_dir, extra_virus)
    if strain_map is not None:
        unknown = sorted(set(strain_map) - set(contigs["contig"]))
        if unknown:
            raise ValueError(f"--strain-clusters names {len(unknown)} contigs that are not in {ctab}: "
                             f"{', '.join(unknown[:3])}")
    feats, c2f = define_features(contigs, weights, cats, namespaces, same_identity, clusters, only, strain_map)
    # a reference split into operational clusters: its genome stays one feature of its own next to the clusters
    split_parent = {(r.namespace, r.accession): r.parent_feature for r in feats.itertuples()
                    if isinstance(r.parent_feature, str) and r.parent_feature}

    # reference genomes hit by nucleotide bins: RefSeq viral genomes and the run's extra virus FASTA, each from
    # its own file and namespace, so every "X-like" bin competes with the genome it was named after
    ref_rows = []
    ref_feature: dict[tuple, str] = {}          # (source, accession) -> feature id
    ref_files: dict[str, Path] = {}
    if add_reference_genomes:
        for src_key, fasta in (("virus", virus), ("extra_virus", extra_virus)):
            ns = namespaces.get(src_key)
            want = {r.accession: r for r in feats.itertuples() if ns and r.namespace == ns
                    and r.identity_class in ("same", "like")}
            if not want:
                continue
            if not fasta or not Path(fasta).exists():
                opt = "--virus FASTA" if src_key == "virus" else "the run's extra_virus_fasta"
                raise FileNotFoundError(f"{src_key} reference {fasta or '(not set)'} not found, but {len(want)} "
                                        f"features are named after its genomes (provide {opt}, or "
                                        "--no-reference-genomes)")
            ids = work / f"reference_ids_{src_key}.txt"
            ids.write_text("\n".join(sorted(want)) + "\n")
            out = work / f"reference_genomes_{src_key}.fa"
            sh(f"seqkit grep -f {shlex.quote(str(ids))} {shlex.quote(str(fasta))} > {shlex.quote(str(out))}")
            ref_files[src_key] = out
            found = {acc: len(seq) for acc, seq in read_fasta(out)}
            bad = sorted(a for a in want if not found.get(a))
            if bad:
                raise ValueError(f"{len(bad)} reference genome(s) named by the features were not found, or have no "
                                 f"sequence, in {fasta}: {', '.join(bad[:5])}. Use the reference FASTA the run "
                                 "classified with, or --no-reference-genomes.")
            for acc, seq in read_fasta(out):
                r = want[acc]
                parent = split_parent.get((r.namespace, acc), "")
                if r.identity_class == "same" and not parent:
                    ref_feature[(src_key, acc)] = r.feature_id
                else:
                    fid = f"{r.namespace}|{acc}|reference"
                    ref_feature[(src_key, acc)] = fid
                    ref_rows.append({"feature_id": fid, "label": f"{clean_title(r.best_hit_title, acc)} (reference genome)",
                                     "namespace": r.namespace, "accession": acc, "identity_class": "reference",
                                     "n_sequences": 0, "total_length": 0,
                                     "evidence": ("reference genome of a bin split into operational clusters"
                                                  if parent else "reference genome of a 'like' bin"),
                                     "n_contigs": 0, "weighted_identity": "", "identity_min": "", "identity_max": "",
                                     "mixed_identity": False, "reference_included": True,
                                     "best_hit_title": r.best_hit_title, "parent_feature": parent})
    feats = pd.concat([feats, pd.DataFrame(ref_rows, columns=FEATURE_COLUMNS)], ignore_index=True) if ref_rows else feats

    # write the panel FASTA
    seq2feature: dict[str, str] = {}
    source: dict[str, str] = {}
    extra_len: dict[str, list] = defaultdict(lambda: [0, 0])
    panel_fa = stage / "panel.fa"
    target_contigs = set(c2f)
    # contigs already represented by a decoy reference (host genome, rRNA) are not added a second time
    covered = {f"host:{h['name']}" for h in hosts}           # only categories whose reference is loaded
    if rrna:
        covered.add("rRNA")
    skip_contigs = set(contigs.loc[contigs["category"].isin(covered), "contig"])
    with open(panel_fa, "w") as fo:
        for name, seq in read_fasta(cfa):
            if name in target_contigs:
                key, val = f"contig|{name}", c2f[name]
            elif other_contig_decoy and name not in skip_contigs:
                key, val = f"decoy|other_contigs|{name}", f"{DECOY_PREFIX}other_contigs"
            else:
                continue
            seq2feature[key] = val
            source[key] = "contig" if name in target_contigs else "decoy"
            fo.write(f">{key}\n{seq}\n")
        for src_key, path in ref_files.items():
            for acc, seq in read_fasta(path):
                fid = ref_feature[(src_key, acc)]
                key = f"ref|{acc}" if src_key == "virus" else f"ref|{src_key}|{acc}"
                seq2feature[key], source[key] = fid, "reference"
                extra_len[fid][0] += 1
                extra_len[fid][1] += len(seq)
                fo.write(f">{key}\n{seq}\n")
        extra_rows = []
        if extra_targets:
            ns = f"extra:{slug(Path(extra_targets).name)}"
            titles = {}
            with open_text(extra_targets) as fh:
                for line in fh:
                    if line.startswith(">"):
                        titles[line[1:].split()[0]] = line[1:].strip()
            for sid, seq in read_fasta(extra_targets):
                fid = f"{ns}|{sid}|reference"
                key = f"extra|{sid}"
                seq2feature[key], source[key] = fid, "extra"
                fo.write(f">{key}\n{seq}\n")
                extra_rows.append({"feature_id": fid, "label": clean_title(titles[sid], sid), "namespace": ns,
                                   "accession": sid, "identity_class": "reference", "n_sequences": 1,
                                   "total_length": len(seq), "evidence": "extra target (user FASTA)", "n_contigs": 0,
                                   "weighted_identity": "", "identity_min": "", "identity_max": "",
                                   "mixed_identity": False, "reference_included": True, "best_hit_title": titles[sid],
                                   "parent_feature": ""})
        if extra_rows:
            feats = pd.concat([feats, pd.DataFrame(extra_rows, columns=FEATURE_COLUMNS)], ignore_index=True)
        decoys = []
        for h in hosts:
            cls = f"{DECOY_PREFIX}host:{h['name']}"
            n, total = _append_fasta(Path(h["fasta"]), fo, f"decoy|host:{h['name']}|", seq2feature, lambda _s, c=cls: c)
            decoys.append({"decoy": cls, "source": str(h["fasta"]), "n_sequences": n, "total_length": total})
        if rrna:
            n, total = _append_fasta(Path(rrna), fo, "decoy|rrna|", seq2feature, lambda _s: f"{DECOY_PREFIX}rrna")
            decoys.append({"decoy": f"{DECOY_PREFIX}rrna", "source": str(rrna), "n_sequences": n, "total_length": total})
    if other_contig_decoy:
        n = sum(1 for v in seq2feature.values() if v == f"{DECOY_PREFIX}other_contigs")
        decoys.append({"decoy": f"{DECOY_PREFIX}other_contigs", "source": str(cfa), "n_sequences": n, "total_length": ""})
    for k in seq2feature:
        source.setdefault(k, "decoy")

    for fid, (n, total) in extra_len.items():
        i = feats.index[feats["feature_id"] == fid][0]
        feats.loc[i, "n_sequences"] = int(feats.loc[i, "n_sequences"]) + n
        feats.loc[i, "total_length"] = int(feats.loc[i, "total_length"]) + total
        feats.loc[i, "reference_included"] = True
    if feats["feature_id"].duplicated().any():
        raise StepError("duplicate feature ids: " + ", ".join(feats.loc[feats["feature_id"].duplicated(), "feature_id"]))
    feats.to_csv(stage / "features.tsv", sep="\t", index=False)
    with open(stage / "seq2feature.tsv", "w") as fo:
        fo.write("sequence\tfeature_or_decoy\tkind\tsource\n")
        for k, v in seq2feature.items():
            fo.write(f"{k}\t{v}\t{'decoy' if v.startswith(DECOY_PREFIX) else 'target'}\t{source[k]}\n")
    pd.DataFrame(decoys).to_csv(stage / "decoys.tsv", sep="\t", index=False)

    idx = stage / "index"
    idx.mkdir()
    q = shlex.quote
    sh(f"bowtie2-build --threads {threads} {q(str(panel_fa))} {q(str(idx / 'panel'))} > {q(str(idx / 'build.log'))} 2>&1")
    if not list(idx.glob("panel.*bt2*")):
        raise StepError(f"bowtie2-build wrote no index files in {idx}")
    # the panel is complete: move it into place, then write the provenance (with artifact identities) last
    for name in PANEL_ARTIFACTS:
        dst = outdir / name
        if dst.is_symlink():
            dst.unlink()
        elif dst.is_dir():
            shutil.rmtree(dst)                    # the previous index of this panel folder (owned output)
        os.replace(stage / name, dst)
    remove_staging_dir(stage)
    panel_fa = outdir / "panel.fa"

    inputs = {"contigs.fa": cfa, "contig_category.tsv": ctab}
    if cc.exists():
        inputs["contig_counts.tsv"] = cc
    if "virus" in ref_files:
        inputs["virus_reference"] = Path(virus)
    if "extra_virus" in ref_files:
        inputs["extra_virus_reference"] = Path(extra_virus)
    if rrna:
        inputs["rrna_reference"] = Path(rrna)
    for h in hosts:
        inputs[f"host:{h['name']}"] = Path(h["fasta"])
    if extra_targets:
        inputs["extra_targets"] = Path(extra_targets)
    if strain_clusters:
        inputs["strain_clusters"] = Path(strain_clusters)
    artifacts = artifact_identities(outdir)
    prov = {
        "unclassified_meta_version": __version__, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "panel_id": hashlib.sha256(json.dumps(artifacts, sort_keys=True).encode()).hexdigest()[:32],
        "artifacts": artifacts,
        "parameters": {"categories": sorted(cats), "contig_list": contig_list, "same_identity": same_identity,
                       "add_reference_genomes": add_reference_genomes, "other_contig_decoy": other_contig_decoy,
                       "host_decoys": [h["name"] for h in hosts], "rrna_decoy": bool(rrna),
                       "unmatched_cluster_min_identity": 90.0, "unmatched_cluster_min_cover": 0.5,
                       "strain_clusters": str(Path(strain_clusters).resolve()) if strain_clusters else None,
                       "strain_resolution": strain_resolution if strain_clusters else None},
        "namespaces": namespaces, "refs_version": refs_version,
        "inputs": {k: {"path": str(Path(v).resolve()), "sha256": sha256(Path(v))} for k, v in inputs.items()},
        "tools": {"bowtie2": tool_version(["bowtie2", "--version"]).split("version")[-1].strip(),
                  "blastn": tool_version(["blastn", "-version"]), "seqkit": tool_version(["seqkit", "version"])},
        "n_features": int(len(feats)), "n_target_sequences": sum(1 for v in seq2feature.values() if not v.startswith(DECOY_PREFIX)),
        "n_decoy_sequences": sum(1 for v in seq2feature.values() if v.startswith(DECOY_PREFIX)),
    }
    with staged(outdir / "provenance.json") as tmp:           # the seal, written last
        tmp.write_text(json.dumps(prov, indent=2) + "\n")
    log.info("panel: %d features, %d target sequences, %d decoy sequences -> %s", prov["n_features"],
             prov["n_target_sequences"], prov["n_decoy_sequences"], outdir)
    if len(feats) < 50:
        log.warning("the panel has %d features; ANCOM-BC2 estimates its sample bias from the features and warns "
                    "below 50 (see README, panel workflow, limits)", len(feats))
    return outdir


# ---------------------------------------------------------------- panel identity

PANEL_ARTIFACTS = ("panel.fa", "features.tsv", "seq2feature.tsv", "decoys.tsv", "index")


def artifact_identities(panel: Path) -> dict[str, str]:
    """sha256 of the panel FASTA, its tables and every Bowtie 2 index file."""
    panel = Path(panel)
    out = {}
    for name in PANEL_ARTIFACTS[:-1]:
        out[name] = sha256(panel / name)
    for f in sorted((panel / "index").glob("panel.*bt2*")):
        out[f"index/{f.name}"] = sha256(f)
    return out


def verify_panel(panel: Path) -> dict:
    """The provenance of a complete panel whose files are unchanged since panel-build wrote it."""
    panel = Path(panel)
    pf = panel / "provenance.json"
    if not pf.exists():
        raise StepError(f"{pf} not found: the panel is incomplete (an interrupted panel-build?) or not a panel "
                        "folder; run panel-build again")
    prov = json.loads(pf.read_text())
    expected = prov.get("artifacts")
    if not expected:
        raise StepError(f"{pf} records no artifact checksums (panel made by an older version); run panel-build again")
    missing = [n for n in expected if not (panel / n).exists()]
    if missing:
        raise StepError(f"panel files missing: {', '.join(missing)}; run panel-build again")
    found = artifact_identities(panel)
    changed = sorted(n for n in set(expected) | set(found) if expected.get(n) != found.get(n))
    if changed:
        raise StepError(f"panel files changed since panel-build ({', '.join(changed[:5])}); run panel-build again")
    return prov


# ---------------------------------------------------------------- quant tables: sealed provenance

QUANT_ARTIFACTS = ("counts.tsv", "ledger.tsv", "features.tsv")
# per-sample fragment classes (the candidate features/decoys of each fragment's best alignments) and the fragments
# counted for one feature only; written on every run, read by `em-uncertainty` (bootstrap, unique support)
CLASS_ARTIFACTS = ("equivalence_classes.tsv", "unique_counts.tsv")


def write_class_tables(folder: Path, samples: list, classes: dict, feature_ids: list) -> None:
    """equivalence_classes.tsv: sample, candidates (';'-joined, sorted; empty = unmapped), fragments - every
    fragment of every sample once. unique_counts.tsv: feature x sample fragments whose best alignments fall in
    that one feature only (the directly observed counts, never replaced by EM shares)."""
    rows = [(s.name, ";".join(c), int(n)) for s in samples for c, n in sorted(classes[s.name].items()) if n]
    pd.DataFrame(rows, columns=["sample", "candidates", "fragments"]).to_csv(
        folder / "equivalence_classes.tsv", sep="\t", index=False)
    uniq = pd.DataFrame({"feature_id": feature_ids})
    for s in samples:
        own = {c[0]: int(n) for c, n in classes[s.name].items() if len(c) == 1}
        uniq[s.name] = [own.get(f, 0) for f in feature_ids]
    uniq.to_csv(folder / "unique_counts.tsv", sep="\t", index=False)


def read_class_tables(qdir: Path) -> tuple[dict, "pd.DataFrame"]:
    """({sample: {candidate tuple: fragments}}, unique counts) of a sealed panel-quant folder; tuple () holds the
    unmapped fragments."""
    rec = verify_quant(Path(qdir))
    missing = [n for n in CLASS_ARTIFACTS if n not in rec["artifacts"]]
    if missing:
        raise ValueError(f"{qdir} has no {', '.join(missing)}; run panel-quant again with this version")
    ec = pd.read_csv(Path(qdir) / "equivalence_classes.tsv", sep="\t", dtype={"sample": str, "candidates": str},
                     keep_default_na=False)
    out: dict = {}
    for s, c, n in ec.itertuples(index=False):
        key = tuple(c.split(";")) if c else ()
        d = out.setdefault(s, {})
        d[key] = d.get(key, 0) + int(n)
    uniq = pd.read_csv(Path(qdir) / "unique_counts.tsv", sep="\t", dtype={"feature_id": str})
    return out, uniq


def write_quant_provenance(qdir: Path, record: dict, extra: list | None = None) -> None:
    """Seal the published quant tables: provenance.json (the generation manifest) is written last. The extract
    summary the counts are reconciled with is recorded with its checksum (verified by the merge). `extra`: further
    files (relative paths), e.g. the imputation tables."""
    qdir = Path(qdir)
    if record.get("input_summary"):
        record = {**record, "input_summary_sha256": manifest.sha256(record["input_summary"])}
    manifest.seal(qdir, [n for n in QUANT_ARTIFACTS if (qdir / n).exists()] + list(extra or []), record,
                  name="provenance.json")


def verify_quant(qdir: Path) -> dict:
    """Provenance of a quant folder whose tables are exactly the ones panel-quant published (ValueError otherwise)."""
    return manifest.verify(qdir, "provenance.json", required=("counts.tsv",), what="quant tables")


# ---------------------------------------------------------------- panel-quant

def load_seq2feature(panel: Path) -> dict[str, str]:
    out = {}
    with open(panel / "seq2feature.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = f[1]
    return out


def detect_reads_scope(samples: list[Sample]) -> tuple[str, Path | None]:
    """'unclassified' when every sample's reads lie in the unclassified_reads/ folder written by `extract` or `run`
    (recognised by tables/input_summary.tsv next to it), otherwise 'all'. Returns (scope, input summary or None)."""
    parents = {Path(s.r1).absolute().parent for s in samples}
    if len(parents) == 1:
        d = parents.pop()
        summary = d.parent / "tables" / "input_summary.tsv"
        if d.name == "unclassified_reads" and summary.exists():
            return "unclassified", summary
    return "all", None


AMBIGUOUS_MODES = ("ledger", "em")


def _ledger_policy(classes: dict) -> tuple[dict, dict]:
    """(fragments per feature, ledger outcomes) under the 'feature' policy, from the fragment classes."""
    feats, ledger = Counter(), Counter()
    for cls, n in classes.items():
        o = outcome_of_class(cls)
        if o in (UNMAPPED, AMBIGUOUS, AMBIGUOUS_DECOY, DECOY_MIXED) or o.startswith(DECOY_PREFIX):
            ledger[o] += n
        else:
            feats[o] += n
    return dict(feats), dict(ledger)


def quantify(panel: Path, samples: list[Sample], outdir: Path, threads: int = 8, max_alignments: int = 10,
             preset: str = "--sensitive", reads_scope: str = "auto", input_summary: str | Path | None = None,
             ambiguous: str = "ledger", imputations: int = 0, seed: int = 1) -> Path:
    """Count every fragment of every sample once. Writes counts.tsv, ledger.tsv, features.tsv, provenance.json.

    reads_scope records what the input reads are: 'all' reads of the samples, or the 'unclassified' reads written
    by `extract` (the reads a classifier did not assign, plus low-score reads). `taxa-table --merge-panel` needs it
    to know whether the panel and the classifier's taxa were counted on the same or on disjoint reads.

    ambiguous 'ledger' (default): a fragment whose best alignments fall in several features (or in features and
    decoys) is booked in the ledger and counted for none. 'em': it is shared among those candidates, decoys
    included, in proportion to their maximum-likelihood abundance in the sample (em.py); the counts are then
    expected fragments with 6 decimals, and `imputations` > 0 also writes that many integer tables drawn from the
    posterior assignment (imputations/counts_NNN.tsv, ledger_NNN.tsv) for `da --imputations`.
    """
    from . import em
    import numpy as np
    if ambiguous not in AMBIGUOUS_MODES:
        raise ValueError(f"--ambiguous must be one of {', '.join(AMBIGUOUS_MODES)}")
    imputations = int(imputations or 0)
    if int(seed) < 0:
        raise ValueError("--seed must be 0 or a positive whole number")
    if imputations < 0 or (imputations and ambiguous != "em"):
        raise ValueError("--imputations needs --ambiguous em and a positive number")
    panel, outdir = Path(panel), Path(outdir)
    if reads_scope not in ("auto", "all", "unclassified"):
        raise ValueError("--reads-scope must be auto, all or unclassified")
    auto_scope, auto_summary = detect_reads_scope(samples)
    if reads_scope == "auto":
        reads_scope = auto_scope
    if reads_scope == "unclassified" and input_summary is None:
        input_summary = auto_summary
    if input_summary is not None and not Path(input_summary).exists():
        raise FileNotFoundError(f"input summary not found: {input_summary}")
    if reads_scope == "all":
        input_summary = None
    log.info("reads scope: %s%s", reads_scope, f" (extract summary: {input_summary})" if input_summary else "")
    check_sample_names([s.name for s in samples])
    refuse_overlap([outdir / n for n in QUANT_ARTIFACTS + ("provenance.json",)],
                   [f for s in samples for f in (s.r1, s.r2)], "panel-quant")
    manifest.verify_consumed([f for s in samples for f in (s.r1, s.r2)], "reads")   # extracted reads: sealed set
    prov = verify_panel(panel)
    outdir.mkdir(parents=True, exist_ok=True)
    # unpublish the previous result first: until the new provenance is written, nothing may use these tables
    if (outdir / "imputations").is_symlink():
        raise ValueError(f"{outdir / 'imputations'} is a link; the imputation tables must be written inside the "
                         "output folder itself")
    previous = manifest.previous_artifacts(outdir, "provenance.json")
    manifest.begin(outdir, "panel-quant", name="provenance.json")
    feats = pd.read_csv(panel / "features.tsv", sep="\t")
    seq2feature = load_seq2feature(panel)
    work = outdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    options = f"{preset} -k {int(max_alignments)}"
    counts, ledger_rows = {}, []
    decoy_classes = sorted({v for v in seq2feature.values() if v.startswith(DECOY_PREFIX)})
    fits, inputs, em_info, all_classes = {}, {}, {}, {}
    for s in samples:
        key = params_key("panel-quant", "classes-v1", prov["panel_id"], file_id(s.r1), file_id(s.r2), options)
        cache = work / f"{s.name}.json"
        rec = json.loads(cache.read_text()) if cache.exists() else None
        if rec is not None and rec.get("key") != key:
            log.info("%s: cached counts were made with another panel, input or parameters; recounting", s.name)
            rec = None
        if rec is None:
            res = run_bowtie2(panel / "index" / "panel", s.r1, s.r2, work / f"{s.name}.bowtie2.log", threads,
                              options, lambda fh: count_classes(fh, seq2feature))
            classes, seen = res["result"]
            if res["input_reads"] is not None and res["input_reads"] != seen:
                raise StepError(f"{s.name}: {seen} fragments in the SAM stream but Bowtie 2 read {res['input_reads']}")
            rec = {"key": key, "sample": s.name, "input_fragments": seen,
                   "classes": [[list(c), n] for c, n in sorted(classes.items())],
                   "wall_s": res["wall_s"], "max_rss_gb": res["max_rss_gb"]}
            cache.write_text(json.dumps(rec))
        classes = {tuple(c): n for c, n in rec["classes"]}
        all_classes[s.name] = classes
        inputs[s.name] = (rec["input_fragments"], classes.get((), 0))
        rec["features"], rec["ledger"] = _ledger_policy(classes)
        shared = 0
        if ambiguous == "em":
            mapped = {c: n for c, n in classes.items() if c}
            fit = em.expected_counts(mapped)
            if not fit.converged:
                raise StepError(f"{s.name}: the EM assignment did not converge in {fit.iterations} iterations")
            fits[s.name] = (mapped, fit)
            shared = sum(n for c, n in mapped.items() if len(c) > 1)
            rounded = em.round_preserving(fit.counts, 6)     # the written shares add up to the mapped fragments
            rec["features"] = {t: v for t, v in rounded.items() if not t.startswith(DECOY_PREFIX)}
            rec["ledger"] = {**{t: v for t, v in rounded.items() if t.startswith(DECOY_PREFIX)},
                             UNMAPPED: classes.get((), 0)}
        assigned = round(sum(rec["features"].values()), 6)
        row = {"sample": s.name, "paired": s.paired, "input_fragments": rec["input_fragments"],
               "assigned_to_features": assigned, AMBIGUOUS: rec["ledger"].get(AMBIGUOUS, 0),
               AMBIGUOUS_DECOY: rec["ledger"].get(AMBIGUOUS_DECOY, 0)}
        for d in decoy_classes:
            row[d] = rec["ledger"].get(d, 0)
        row[DECOY_MIXED] = rec["ledger"].get(DECOY_MIXED, 0)
        row[UNMAPPED] = rec["ledger"].get(UNMAPPED, 0)
        outcomes = sum(v for k, v in row.items() if k not in ("sample", "paired", "input_fragments"))
        if abs(outcomes - rec["input_fragments"]) > (1e-6 * (len(feats) + len(decoy_classes) + 10)
                                                     if ambiguous == "em" else 0):
            raise StepError(f"{s.name}: outcomes sum to {outcomes}, input fragments {rec['input_fragments']}")
        em_info[s.name] = {"shared_fragments": shared, "iterations": fits[s.name][1].iterations,
                           "converged": fits[s.name][1].converged} if ambiguous == "em" else None
        row["wall_s"], row["max_rss_gb"] = rec["wall_s"], rec["max_rss_gb"]
        ledger_rows.append(row)
        counts[s.name] = rec["features"]
    table = pd.DataFrame({"feature_id": feats["feature_id"]})
    for s in samples:
        unknown = set(counts[s.name]) - set(feats["feature_id"])
        if unknown:
            raise StepError(f"{s.name}: counts for features missing from features.tsv: {sorted(unknown)[:3]}")
        value = (lambda v: round(float(v), 6)) if ambiguous == "em" else int
        table[s.name] = [value(counts[s.name].get(f, 0)) for f in feats["feature_id"]]
    stage = make_staging_dir(outdir, "panel-quant")
    table.to_csv(stage / "counts.tsv", sep="\t", index=False, float_format="%.6f" if ambiguous == "em" else None)
    pd.DataFrame(ledger_rows).to_csv(stage / "ledger.tsv", sep="\t", index=False)
    feats.to_csv(stage / "features.tsv", sep="\t", index=False)
    write_class_tables(stage, samples, all_classes, list(feats["feature_id"]))
    imp_names = list(CLASS_ARTIFACTS)
    if imputations:
        (stage / "imputations").mkdir()
    draws = {}
    for si, s in enumerate(samples):
        if not imputations:
            break
        mapped, fit = fits[s.name]
        # one independent stream per sample, and a different one from taxa-table's (side 2 = panel)
        draws[s.name] = em.imputations(mapped, fit, imputations, np.random.default_rng([int(seed), 2, si]))
    for m in range(imputations):
        it = pd.DataFrame({"feature_id": feats["feature_id"]})
        irows = []
        for s in samples:
            draw = draws[s.name][m]
            it[s.name] = [int(draw.get(f, 0)) for f in feats["feature_id"]]
            total, unmapped = inputs[s.name]
            irow = {"sample": s.name, "imputation": m + 1, "seed": int(seed), "input_fragments": total,
                    "assigned_to_features": int(it[s.name].sum()),
                    **{d: int(draw.get(d, 0)) for d in decoy_classes}, UNMAPPED: unmapped}
            if irow["assigned_to_features"] + sum(irow[d] for d in decoy_classes) + unmapped != total:
                raise StepError(f"{s.name}: imputation {m + 1} does not keep every fragment")
            irows.append(irow)
        for kind, df in (("counts", it), ("ledger", pd.DataFrame(irows))):
            name = f"imputations/{kind}_{m + 1:03d}.tsv"
            df.to_csv(stage / name, sep="\t", index=False)
            imp_names.append(name)
    manifest.remove_stale(outdir, previous, imp_names)
    if imputations:
        (outdir / "imputations").mkdir(exist_ok=True)
    for name in list(QUANT_ARTIFACTS) + imp_names:
        os.replace(stage / name, outdir / name)
    remove_staging_dir(stage)
    write_quant_provenance(outdir, {
        "unclassified_meta_version": __version__, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
         "panel": prov, "bowtie2_options": options, "counting_unit": "fragment (read pair or single read)",
         "reads_scope": reads_scope,
         "input_summary": str(Path(input_summary).resolve()) if input_summary else None,
         "policy": ("a fragment counts once for a feature when every reported alignment of either mate falls in "
                    "that feature; otherwise it goes to one ledger outcome") if ambiguous == "ledger" else
                   ("a fragment whose best alignments fall in several features/decoys is shared among them by EM "
                    "(expected counts)"),
         "ambiguous": ambiguous, "imputations": imputations,
         **({"em": em_info} if ambiguous == "em" else {}), **({"seed": int(seed)} if imputations else {}),
         "samples": {s.name: {"r1": file_id(s.r1), "r2": file_id(s.r2)} for s in samples}}, imp_names)
    log.info("counts: %s (%d features x %d samples); ledger: %s", outdir / "counts.tsv", len(feats), len(samples),
             outdir / "ledger.tsv")
    return outdir / "counts.tsv"
