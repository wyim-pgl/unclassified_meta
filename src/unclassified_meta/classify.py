"""Assign each contig one category by an ordered list of evidence layers (first match wins).

Layers (each skipped when its reference is not configured):
  rrna       BLASTN (dc-megablast) to an rRNA FASTA; >= min_cover of the contig at >= min_identity
  organelle  the same against plastid / mitochondrial genomes
  host       minimap2 to each host genome, in list order; >= min_cover of the contig aligned
  virus      BLASTN (>= blastn_min_len bp) to RefSeq viral genomes [+ extra virus FASTA], then DIAMOND blastx
             of the rest to RefSeq viral proteins; optional TBLASTX against the small extra virus FASTA
  kraken2    Kraken2 on the contigs; category from the taxonomy lineage (virus, bacteria, ...)
  protein    DIAMOND blastx to Swiss-Prot Viridiplantae (+ optional host proteome)
Contigs left after all layers are 'unclassified'.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

from .io_utils import read_fasta, write_fasta
from .runner import require, sh

log = logging.getLogger("unclassified_meta")

BLAST_FMT = "6 qseqid sseqid pident length qlen qstart qend evalue bitscore stitle"

# Kraken2 lineage -> category; checked from the leaf upwards, first match wins.
KRAKEN_TAXA = [("9606", "human"), ("10239", "virus"), ("2", "bacteria"), ("2157", "archaea"),
               ("4751", "fungi"), ("33090", "plant"), ("33208", "animal"), ("2759", "other_eukaryote")]


def merged_cover(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by 1-based inclusive intervals (any orientation), overlaps merged."""
    iv = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_s is None or s > cur_e + 1:
            if cur_s is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_s is not None:
        total += cur_e - cur_s + 1
    return total


def parse_blast(path: Path, min_identity: float = 0.0, max_evalue: float = 10.0,
                min_len: int = 0) -> dict[str, dict]:
    """Per query: best hit (by bitscore), its identity and title, and merged query cover of passing HSPs."""
    hsps = defaultdict(list)
    best: dict[str, dict] = {}
    if not Path(path).exists():
        return {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        q, s, pid, alen, qlen, qs, qe, ev, bits = f[:9]
        title = f[9] if len(f) > 9 else s
        if float(pid) < min_identity or float(ev) > max_evalue or int(alen) < min_len:
            continue
        hsps[q].append((int(qs), int(qe)))
        if q not in best or float(bits) > best[q]["bitscore"]:
            best[q] = {"hit": s, "title": title, "identity": float(pid), "bitscore": float(bits),
                       "qlen": int(qlen)}
    for q, b in best.items():
        b["cover"] = merged_cover(hsps[q]) / b["qlen"] if b["qlen"] else 0.0
    return best


def parse_paf_cover(path: Path) -> dict[str, dict]:
    """minimap2 PAF -> per query merged cover fraction, best target and identity (matches / block length)."""
    iv = defaultdict(list)
    best: dict[str, dict] = {}
    if not Path(path).exists():
        return {}
    for line in open(path):
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, qlen, qs, qe, t = f[0], int(f[1]), int(f[2]), int(f[3]), f[5]
        matches, block = int(f[9]), int(f[10])
        iv[q].append((qs + 1, qe))
        if q not in best or matches > best[q]["matches"]:
            best[q] = {"hit": t, "title": t, "matches": matches,
                       "identity": round(100 * matches / block, 2) if block else 0.0, "qlen": qlen}
    for q, b in best.items():
        b["cover"] = merged_cover(iv[q]) / b["qlen"] if b["qlen"] else 0.0
    return best


def load_ktaxonomy(path: Path) -> tuple[dict, dict]:
    """Kraken2 ktaxonomy.tsv ('taxid | parent | rank | depth | name') -> (parent, name)."""
    parent, name = {}, {}
    for line in open(path):
        f = [x.strip() for x in line.split("|")]
        if len(f) >= 2:
            parent[f[0]] = f[1]
            name[f[0]] = f[4] if len(f) > 4 else f[0]
    return parent, name


def kraken_category(taxid: str, parent: dict) -> str:
    seen, t = 0, taxid
    while t in parent and seen < 200:
        for tid, cat in KRAKEN_TAXA:
            if t == tid:
                return cat
        if parent[t] == t:
            break
        t, seen = parent[t], seen + 1
    return "other"


def blast_db(fasta: str | Path, fallback: Path) -> Path:
    """Use a BLAST database built next to the FASTA (download-refs does this); otherwise build one in the workdir."""
    fasta = Path(fasta)
    for ext in (".nal", ".nsq", ".00.nsq", ".ndb"):
        if Path(str(fasta) + ext).exists():
            return fasta
    if not Path(str(fallback) + ".ndone").exists():
        sh(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(fallback)])
        Path(str(fallback) + ".ndone").touch()
    return fallback


def assign_by_priority(contigs: list[str], layer_hits: list[tuple[str, dict[str, dict]]]) -> dict[str, dict]:
    """layer_hits: [(category, {contig: info}), ...] in priority order. First layer that hit wins."""
    out = {}
    for c in contigs:
        for category, hits in layer_hits:
            if c in hits:
                out[c] = {"category": category, **hits[c]}
                break
        else:
            out[c] = {"category": "unclassified"}
    return out


class Classifier:
    """Runs the layers on the contigs still unassigned, so later (slower) layers see fewer contigs."""

    def __init__(self, cfg: dict, workdir: Path, threads: int):
        self.cfg, self.work, self.threads = cfg, workdir, threads
        self.work.mkdir(parents=True, exist_ok=True)

    def _blast(self, query: Path, ref: str, tag: str, task: str = "dc-megablast",
               program: str = "blastn", evalue: float = 1e-10) -> Path:
        require(program, "makeblastdb")
        db = blast_db(ref, self.work / f"{tag}_db")
        out = self.work / f"{tag}.{program}.tsv"
        cmd = [program, "-query", str(query), "-db", str(db), "-evalue", str(evalue),
               "-max_target_seqs", "5", "-num_threads", str(self.threads), "-outfmt", BLAST_FMT,
               "-out", str(out)]
        if program == "blastn":
            cmd[1:1] = ["-task", task]
        sh(cmd)
        return out

    def _diamond(self, query: Path, db: str | list[str], tag: str, evalue: float) -> Path:
        """DIAMOND blastx; `db` is a .dmnd, a protein FASTA, or a list of them (FASTA/.dmnd mixed is built into one)."""
        require("diamond")
        dbs = [db] if isinstance(db, str) else [d for d in db if d]
        if len(dbs) == 1 and str(dbs[0]).endswith(".dmnd"):
            dmnd = Path(dbs[0])
        else:
            dmnd = self.work / f"{tag}.dmnd"
            if not dmnd.exists():
                fasta = self.work / f"{tag}.faa"
                with open(fasta, "w") as fo:
                    for d in dbs:
                        if str(d).endswith(".dmnd"):
                            sh(["diamond", "getseq", "-d", str(d)], stdout=self.work / "_getseq.faa")
                            fo.write((self.work / "_getseq.faa").read_text())
                        else:
                            fo.writelines(open(d))
                sh(["diamond", "makedb", "--in", str(fasta), "-d", str(dmnd.with_suffix("")), "-p", str(self.threads)])
        out = self.work / f"{tag}.diamond.tsv"
        sh(["diamond", "blastx", "-q", str(query), "-d", str(dmnd), "-p", str(self.threads),
            "-e", str(evalue), "-k", "1", "--sensitive", "--quiet",
            "-f", "6", "qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend", "evalue",
            "bitscore", "stitle", "-o", str(out)])
        return out

    def layer(self, name: str, query: Path) -> list[tuple[str, dict]]:
        """Return [(category, hits)] for one configured layer (host can return several)."""
        refs, th = self.cfg["references"], self.cfg["thresholds"]
        if name in ("rrna", "organelle"):
            if not refs.get(name):
                return []
            t = th[name]
            out = self._blast(query, refs[name], name, evalue=t["evalue"])
            hits = {q: {**b, "method": "blastn"} for q, b in parse_blast(out, t["min_identity"], t["evalue"]).items()
                    if b["cover"] >= t["min_cover"]}
            return [("rRNA" if name == "rrna" else "organelle", hits)]
        if name == "host":
            res = []
            for h in refs.get("host") or []:
                require("minimap2")
                paf = self.work / f"host_{h['name']}.paf"
                sh(["minimap2", "-x", th["host"]["preset"], "-t", str(self.threads), "--secondary=no", "-c",
                    str(h["fasta"]), str(query)], stdout=paf)
                hits = {q: {**b, "method": "minimap2"} for q, b in parse_paf_cover(paf).items()
                        if b["cover"] >= th["host"]["min_cover"]}
                res.append((f"host:{h['name']}", hits))
            return res
        if name == "virus":
            t = th["virus"]
            hits: dict[str, dict] = {}
            if refs.get("virus"):
                for q, h in parse_blast(self._blast(query, refs["virus"], "virus_genomes", evalue=t["evalue"]),
                                        max_evalue=t["evalue"], min_len=t["blastn_min_len"]).items():
                    hits[q] = {**h, "method": "blastn"}
            if refs.get("extra_virus_fasta"):
                for q, h in parse_blast(self._blast(query, refs["extra_virus_fasta"], "virus_extra",
                                                    evalue=t["evalue"]),
                                        max_evalue=t["evalue"], min_len=t["blastn_min_len"]).items():
                    if q not in hits or h["bitscore"] > hits[q]["bitscore"]:
                        hits[q] = {**h, "method": "blastn_extra"}
            if refs.get("virus_protein"):
                rest = self._subset(query, set(hits), "virus_blastx_query.fa")
                if rest:
                    for q, h in parse_blast(self._diamond(rest, refs["virus_protein"], "virus_protein", t["evalue"]),
                                            max_evalue=t["evalue"]).items():
                        hits.setdefault(q, {**h, "method": "diamond_blastx"})
                self._drop_cellular_homologs(hits, query)
            if refs.get("extra_virus_fasta") and t.get("tblastx_extra"):
                rest = self._subset(query, set(hits), "virus_tblastx_query.fa")
                if rest:
                    for q, h in parse_blast(self._blast(rest, refs["extra_virus_fasta"], "virus_extra",
                                                        program="tblastx", evalue=t["evalue"]),
                                            max_evalue=t["evalue"]).items():
                        hits.setdefault(q, {**h, "method": "tblastx_extra"})
            return [("virus", hits)] if (refs.get("virus") or refs.get("virus_protein")
                                         or refs.get("extra_virus_fasta")) else []
        if name == "kraken2":
            db = refs.get("kraken2_db")
            if not db:
                return []
            require("kraken2")
            out = self.work / "kraken2_contigs.tsv"
            sh(["kraken2", "--db", str(db), "--threads", str(self.threads), "--output", str(out),
                "--report", str(self.work / "kraken2_contigs.report"), str(query)])
            tax = refs.get("kraken2_taxonomy") or str(Path(db) / "ktaxonomy.tsv")
            parent, tname = load_ktaxonomy(Path(tax)) if Path(tax).exists() else ({}, {})
            if not parent:
                log.warning("no Kraken2 taxonomy file at %s; Kraken2 contig hits are reported as 'kraken2_hit'", tax)
            per_cat = defaultdict(dict)
            for line in open(out):
                f = line.split("\t")
                if f[0] == "C":
                    cat = kraken_category(f[2], parent) if parent else "kraken2_hit"
                    per_cat[cat][f[1]] = {"hit": f[2], "title": tname.get(f[2], f[2]), "method": "kraken2"}
            merged = {}
            for cat, hits in per_cat.items():
                for q, h in hits.items():
                    merged[q] = {**h, "kraken_category": cat}
            return [("kraken2", merged)]
        if name == "protein":
            dbs = [refs.get("protein"), refs.get("host_proteome")]
            if not any(dbs):
                return []
            dbs = [d for d in dbs if d]
            out = self._diamond(query, dbs[0] if len(dbs) == 1 else dbs, "protein", th["protein"]["evalue"])
            return [("protein", {q: {**h, "method": "diamond_blastx"}
                                 for q, h in parse_blast(out, max_evalue=th["protein"]["evalue"]).items()})]
        raise ValueError(name)

    def _drop_cellular_homologs(self, hits: dict[str, dict], query: Path) -> None:
        """Viral protein sets (giant viruses especially) carry homologs of cellular genes. A contig assigned to a
        virus only by DIAMOND keeps the label when its viral bitscore beats its best plant-protein bitscore."""
        cand = {q for q, h in hits.items() if h.get("method") == "diamond_blastx"}
        dbs = [d for d in (self.cfg["references"].get("protein"), self.cfg["references"].get("host_proteome")) if d]
        if not cand or not dbs:
            return
        sub = self.work / "virus_blastx_hits.fa"
        write_fasta([(n, s) for n, s in read_fasta(query) if n in cand], sub)
        out = self._diamond(sub, dbs[0] if len(dbs) == 1 else dbs, "virus_vs_plant",
                            self.cfg["thresholds"]["protein"]["evalue"])
        plant = parse_blast(out, max_evalue=self.cfg["thresholds"]["protein"]["evalue"])
        for q in cand:
            if q in plant and plant[q]["bitscore"] >= hits[q]["bitscore"]:
                log.info("%s: plant protein %s (bits %.0f) beats viral %s (bits %.0f); not called virus",
                         q, plant[q]["hit"], plant[q]["bitscore"], hits[q]["hit"], hits[q]["bitscore"])
                del hits[q]

    def _subset(self, fasta: Path, exclude: set[str], name: str) -> Path | None:
        recs = [(n, s) for n, s in read_fasta(fasta) if n not in exclude]
        if not recs:
            return None
        path = self.work / name
        write_fasta(recs, path)
        return path

    def run(self, contigs: Path) -> dict[str, dict]:
        names = [n for n, _ in read_fasta(contigs)]
        lengths = {n: len(s) for n, s in read_fasta(contigs)}
        assigned: dict[str, dict] = {}
        query = contigs
        for layer in self.cfg["priority"]:
            if query is None:
                break
            for category, hits in self.layer(layer, query):
                for q, info in hits.items():
                    if q in assigned:
                        continue
                    cat = info.pop("kraken_category", None) if category == "kraken2" else None
                    assigned[q] = {"category": cat or category, "layer": layer, **info}
                log.info("layer %-9s %-22s %6d contigs", layer, category, len(hits))
            query = self._subset(contigs, set(assigned), f"remaining_after_{layer}.fa")
        out = {}
        for n in names:
            info = assigned.get(n, {"category": "unclassified", "layer": "-"})
            out[n] = {"length": lengths[n], **info}
        return out


def write_contig_table(result: dict[str, dict], path: Path) -> None:
    cols = ["contig", "length", "category", "layer", "method", "best_hit", "best_hit_title", "identity", "cover"]
    with open(path, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for c, r in result.items():
            ident = r.get("identity")
            cov = r.get("cover")
            fo.write("\t".join([
                c, str(r["length"]), r["category"], r.get("layer", "-"),
                r.get("method", r.get("layer", "-")), str(r.get("hit", "")), str(r.get("title", "")).replace("\t", " "),
                f"{ident:.2f}" if isinstance(ident, float) else "", f"{cov:.3f}" if isinstance(cov, float) else "",
            ]) + "\n")
