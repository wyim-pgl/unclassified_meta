"""panel-taxonomy: one taxID per panel feature, from the same classifier and taxonomy as the background taxa.

A panel feature is built from reads the classifier left unclassified (or below its thresholds), yet the same
organism often also has classified reads, which form a taxon row of the background table. Without a link between
the two, one organism is tested twice in `da` (e.g. 'Cactus virus X-like' contigs at 70 % identity next to the
classified taxon 'Opuntia potexvirus A'). Classifying the panel's sequences themselves (contigs are long, so most
are assigned even when their reads were not) gives each feature a taxID in the background's taxonomy;
`taxa-table --merge-panel --panel-taxonomy` then adds the feature to its taxon's row. The reads of the two sides
are disjoint, so the sum is the organism's total.

Steps:
  1. panel-taxonomy --panel PANEL --write-fasta panel_targets.fa      (the target sequences, no decoys)
  2. centrifuger -x INDEX -k 5 -u panel_targets.fa > panel_targets.cfr.tsv  (same index as the reads)
  3. panel-taxonomy --panel PANEL --classified panel_targets.cfr.tsv --nodes ... --names ... -o panel_taxonomy

Consensus: each sequence gets the lowest common ancestor of its accepted taxIDs; a feature gets the deepest
taxon that holds at least --min-support (default 0.8) of the classified sequence length of the feature, provided
that at least --min-classified (default 0.5) of its sequence length was classified at all. Otherwise the feature
keeps no taxID and stays a feature of its own.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from . import manifest
from .extract import centrifuge_reads
from .fragments import DECOY_PREFIX
from .io_utils import normalise_read_id, read_fasta, refuse_overlap, staged
from .taxa import Taxonomy, scan_names

log = logging.getLogger("unclassified_meta")


def target_sequences(panel: Path) -> dict[str, str]:
    """sequence id -> feature id for every target (non-decoy) sequence of a verified panel."""
    from .panel import load_seq2feature, verify_panel
    verify_panel(panel)
    return {s: f for s, f in load_seq2feature(panel).items() if not f.startswith(DECOY_PREFIX)}


def write_fasta(panel: str | Path, out: str | Path) -> Path:
    panel, out = Path(panel), Path(out)
    from .panel import verify_panel
    tag = verify_panel(panel)["panel_id"]          # the generation first, then everything read from it
    targets = target_sequences(panel)
    if panel.resolve() in out.resolve().parents or out.resolve() == panel.resolve():
        raise ValueError(f"--write-fasta {out} lies inside the panel folder {panel}; write it elsewhere")
    # also every file the panel uses through a link (index shards, or a whole linked folder, on other storage)
    import os
    used = [Path(root) / f for root, _, files in os.walk(panel, followlinks=True) for f in files]
    refuse_overlap([out], used, "panel-taxonomy --write-fasta")
    n = 0
    with staged(out) as tmp:
        with open(tmp, "w") as fo:
            for sid, seq in read_fasta(panel / "panel.fa"):
                if sid in targets:
                    fo.write(f">{tag}:{sid}\n")     # binds the classifier output to this panel generation
                    for i in range(0, len(seq), 80):
                        fo.write(seq[i:i + 80] + "\n")
                    n += 1
        # checked before the new file replaces anything: on a failure the previous file stays as it was
        if n != len(targets):
            raise ValueError(f"{panel / 'panel.fa'} holds {n} of the {len(targets)} target sequences")
        if verify_panel(panel)["panel_id"] != tag:        # the panel was rebuilt while it was being read
            raise ValueError(f"{panel} changed while its sequences were written; run --write-fasta again")
    log.info("panel-taxonomy: %d target sequences of %d features -> %s", n, len(set(targets.values())), out)
    return out


def consensus(lengths: dict[str, int], seq_taxid: dict[str, str | None], tax: Taxonomy,
              min_support: float, min_classified: float) -> tuple[str | None, float, float]:
    """(taxID, support, classified fraction) of one feature from its sequences' lengths and taxIDs."""
    total = sum(lengths.values())
    classified = {s: t for s, t in seq_taxid.items() if t is not None and tax.path(t)}
    cl_len = sum(lengths[s] for s in classified)
    frac = cl_len / total if total else 0.0
    if not classified or frac < min_classified:
        return None, 0.0, frac
    support: dict[str, int] = {}
    for s, t in classified.items():
        for node in tax.path(t):
            support[node] = support.get(node, 0) + lengths[s]
    ok = [n for n, v in support.items() if v / cl_len >= min_support]
    best = max(ok, key=lambda n: (len(tax.path(n)), support[n]))     # deepest supported node
    return best, support[best] / cl_len, frac


def build(panel: str | Path, classified: str | Path, nodes: str | Path, names: str | Path, outdir: str | Path,
          min_score: int = 0, min_hitlen: int = 0, min_support: float = 0.8, min_classified: float = 0.5) -> Path:
    panel, outdir = Path(panel), Path(outdir)
    if not (0 <= min_classified <= 1):              # also refuses nan
        raise ValueError("--min-classified must lie between 0 and 1")
    if not (0.5 < min_support <= 1):
        raise ValueError("--min-support must be above 0.5 and at most 1 (a majority of the classified length)")
    import os
    from .panel import verify_panel
    panel_id = verify_panel(panel).get("panel_id")
    taxonomy_sha = {k: manifest.sha256(v) for k, v in (("nodes", nodes), ("names", names))}   # before reading them
    targets = target_sequences(panel)
    lengths = {sid: len(seq) for sid, seq in read_fasta(panel / "panel.fa") if sid in targets}
    # the classifier output carries read IDs as `centrifuge_reads` normalises them (a trailing /1 or /2 is removed)
    by_norm: dict[str, str] = {}
    for sid in targets:
        key = normalise_read_id(sid)
        if key in by_norm:
            raise ValueError(f"target sequences {by_norm[key]!r} and {sid!r} cannot be told apart in classifier "
                             "output (they differ only by a /1 or /2 ending); rename one of them")
        by_norm[key] = sid
    classified_sha = manifest.sha256(classified)                # the calls used are exactly the sealed file
    tax = Taxonomy(nodes)
    seq_taxid: dict[str, str | None] = {}
    for rid, status, taxids in centrifuge_reads(classified, min_score, min_hitlen):
        tag, sep, rest = rid.partition(":")
        if not sep or tag != panel_id:
            raise ValueError(f"{classified}: sequence {rid!r} was not written by --write-fasta for this panel "
                             f"generation (expected the prefix {panel_id}:); write the FASTA again and classify it")
        sid = by_norm.get(rest)
        if sid is None:
            raise ValueError(f"{classified}: sequence {rid!r} is not a target sequence of {panel} (classify the "
                             "FASTA written by --write-fasta)")
        seq_taxid[sid] = (taxids[0] if len(taxids) == 1 else tax.lca(taxids)) if status == 2 else None
    missing = sorted(set(targets) - set(seq_taxid))
    if missing:
        raise ValueError(f"{classified} has no line for {len(missing)} target sequence(s), e.g. {missing[0]}: was "
                         "it made from the FASTA written by --write-fasta?")
    rows = []
    by_feature: dict[str, list[str]] = {}
    for sid, f in targets.items():
        by_feature.setdefault(f, []).append(sid)
    for f, sids in sorted(by_feature.items()):
        t, sup, frac = consensus({s: lengths[s] for s in sids}, {s: seq_taxid[s] for s in sids}, tax,
                                 min_support, min_classified)
        rows.append({"feature_id": f, "taxid": t or "", "rank": tax.rank.get(t, "") if t else "",
                     "support": round(sup, 4), "classified_fraction": round(frac, 4), "n_sequences": len(sids),
                     "n_classified": sum(seq_taxid[s] is not None and bool(tax.path(seq_taxid[s])) for s in sids),
                     "status": "assigned" if t else ("unclassified" if frac == 0 else "insufficient")})
    df = pd.DataFrame(rows)
    ids = {t for t in df["taxid"] if t}
    nm, _ = scan_names(names, want_ids=ids)
    df.insert(2, "name", [nm.get(t, "") if t else "" for t in df["taxid"]])
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "panel_taxonomy.tsv"
    panel_files = [Path(root) / f for root, _, files in os.walk(panel, followlinks=True) for f in files]
    refuse_overlap([out, outdir / manifest.NAME], [classified, nodes, names, *panel_files], "panel-taxonomy")
    if {k: manifest.sha256(v) for k, v in (("nodes", nodes), ("names", names))} != taxonomy_sha:
        raise ValueError("nodes/names changed while panel-taxonomy was reading them; run it again")
    if manifest.sha256(classified) != classified_sha or verify_panel(panel).get("panel_id") != panel_id:
        raise ValueError(f"{classified} or the panel changed while panel-taxonomy was reading them; run it again")
    refuse_overlap([outdir / manifest.MARKER, outdir / "panel-taxonomy.log"], [classified, nodes, names, *panel_files],
                   "panel-taxonomy")
    manifest.begin(outdir, "panel-taxonomy")
    with staged(out) as tmp:
        df.to_csv(tmp, sep="\t", index=False)
    manifest.seal(outdir, ["panel_taxonomy.tsv"], {
        "step": "panel-taxonomy", "panel": str(panel.resolve()), "panel_id": panel_id,
        "nodes_sha256": taxonomy_sha["nodes"], "names_sha256": taxonomy_sha["names"],
        "classified": str(Path(classified).resolve()),
        "classified_sha256": classified_sha, "min_score": min_score, "min_hitlen": min_hitlen,
        "min_support": min_support, "min_classified": min_classified})
    log.info("panel-taxonomy: %d of %d features assigned -> %s", int((df["status"] == "assigned").sum()), len(df), out)
    return out
