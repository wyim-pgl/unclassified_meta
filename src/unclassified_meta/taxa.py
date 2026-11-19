"""taxa-table: single-rank taxon x sample counts from Kraken-style reports, readable by `da`.

Accepted reports: Kraken2 `--report` (6 columns, or 8 with `--report-minimizer-data`), `centrifuger-kreport`
and `centrifuger-quant --output-format 3`. Their rows are a tree: column 2 holds the reads of the whole clade,
column 3 the reads assigned to the node itself, the name is indented by two spaces per level. The first column
(a percentage; an estimated abundance in centrifuger-quant) is never used.

One rank only: rows at different ranks overlap (a species' reads are also in its genus), so the table takes the
clade reads of the rows at ONE rank code. Reads assigned above that rank are reported in the ledger, not added.

Origin: the single-rank table, the lineage drop list and the use of the clade-read column are ideas contributed
by Eli Holschbach. The duplicate-sample
and report-pattern checks, the per-sample ledger and the merge with a panel table were added here.
"""
from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .config import check_sample_names
from .extract import centrifuge_reads
from . import manifest
from .io_utils import make_staging_dir, open_text, refuse_overlap, remove_staging_dir, staged

log = logging.getLogger("unclassified_meta")
RANKS = ("D", "K", "P", "C", "O", "F", "G", "S")
# every read of the input has exactly one of the outcomes after input_fragments
LEDGER_COLS = ("input_fragments", "unclassified", "below_threshold", "dropped_lineages", "kept_at_rank",
               "unresolved_above_rank", "rounding_residual")
RANK_NAMES = {"D": {"superkingdom", "domain", "acellular root"}, "K": {"kingdom"}, "P": {"phylum"}, "C": {"class"},
              "O": {"order"}, "F": {"family"}, "G": {"genus"}, "S": {"species"}}


@dataclass
class Report:
    """What one report contributes: counts per taxID at the chosen rank and the read ledger."""
    counts: dict = field(default_factory=dict)       # taxid -> reads (clade reads at the rank, drops removed)
    names: dict = field(default_factory=dict)        # taxid -> name
    lineages: dict = field(default_factory=dict)     # taxid -> "root;...;name"
    ledger: dict = field(default_factory=dict)
    dropped_found: set = field(default_factory=set)  # drop-list entries seen in this report


def _rows(path: str | Path):
    """Yield (pct, clade_reads, own_reads, rank_code, taxid, depth, name) for every row of a Kraken-style report."""
    with open_text(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                raise ValueError(f"{path}: line {n} has {len(f)} tab-separated columns; a Kraken-style report has 6 "
                                 "(or 8 with minimizer data). Pass the report file, not the per-read output.")
            try:
                pct = float(f[0])
                clade, own = float(f[1]), float(f[2])
            except ValueError:
                if n == 1:
                    continue                     # header line
                raise ValueError(f"{path}: line {n} is not a Kraken-style report row: {line[:80]!r}") from None
            if clade != int(clade) or own != int(own) or clade < 0 or own < 0:
                raise ValueError(f"{path}: line {n}: read counts must be non-negative whole numbers "
                                 f"(found {f[1]!r}, {f[2]!r})")
            raw = f[-1]
            lead = len(raw) - len(raw.lstrip(" "))
            if lead % 2:
                raise ValueError(f"{path}: line {n}: the name is indented by {lead} spaces; expected two per level")
            yield pct, int(clade), int(own), f[-3].strip(), f[-2].strip(), lead // 2, raw.strip()


REPORT_FORMATS = ("kraken2", "centrifuger-kreport", "centrifuger-quant")


def parse_report(path: str | Path, rank: str = "S", drop: tuple = (), panel_lineage: tuple = (), *,
                 report_format: str) -> Report:
    """Counts at one rank plus a ledger in which every read of the report has exactly one outcome.

    report_format is given, never guessed: 'kraken2' (Kraken2 --report, 6 or 8 columns), 'centrifuger-kreport'
    or 'centrifuger-quant' (centrifuger-quant --output-format 3). In every row the clade reads must equal the
    row's own reads plus the clade reads of its children: exactly for kraken2 and centrifuger-kreport, within
    half a read per number involved for centrifuger-quant, which rounds every row on its own (its first column
    is an estimated abundance, not a read share). The rounding left over is the ledger's rounding_residual.

    `panel_lineage` entries are dropped like `drop`, but must not lie below the chosen rank: reads assigned to a
    node between the rank and the lineage (e.g. a genus when the panel covers one of its species) could belong
    to the panel's organism and would be counted in both tables.
    """
    if rank not in RANKS:
        raise ValueError(f"rank must be one of {', '.join(RANKS)} (Kraken rank codes), not {rank!r}")
    if report_format not in REPORT_FORMATS:
        raise ValueError(f"report format must be one of {', '.join(REPORT_FORMATS)}, not {report_format!r}")
    rounds = report_format == "centrifuger-quant"
    strict = {str(d).strip().lower() for d in panel_lineage if str(d).strip()}
    drop_low = {str(d).strip().lower() for d in drop if str(d).strip()} | strict
    rep = Report()
    path_names: list[str] = []
    stack_drop: list[bool] = []          # per depth: is this node inside a dropped lineage
    stack_rank: list[str | None] = []    # per depth: taxid of the enclosing row at the chosen rank
    total = unclassified = dropped = unresolved = n_rows = 0
    seen_any = False
    open_rows: list[list] = []           # per depth: [name, clade, own, sum of children's clade reads, n children]

    def close(level: int) -> None:
        """Check clade = own + children for the rows deeper than `level`, which are now complete."""
        while len(open_rows) > level:
            name_, clade_, own_, kids, n_kids = open_rows.pop()
            diff = clade_ - own_ - kids
            bound = 0.5 * (n_kids + 2) if rounds else 0
            if abs(diff) > bound:
                raise ValueError(f"{path}: row {name_!r}: clade reads {clade_} != own reads {own_} + children's clade "
                                 f"reads {kids} ({report_format} report); is the file complete, is column 2 the clade "
                                 "read count and is --report-format right?")
            if open_rows:
                open_rows[-1][3] += clade_
                open_rows[-1][4] += 1

    for pct, clade, own, code, taxid, depth, name in _rows(path):
        n_rows += 1
        seen_any = True
        if depth > len(path_names):
            raise ValueError(f"{path}: row {name!r} is indented deeper than its predecessor allows")
        close(depth)
        open_rows.append([name, clade, own, 0, 0])
        del path_names[depth:], stack_drop[depth:], stack_rank[depth:]
        is_unclassified = code == "U" or (taxid == "0" and name.lower() == "unclassified")
        if depth == 0:
            total += clade
        if is_unclassified:
            if depth != 0:
                raise ValueError(f"{path}: an 'unclassified' row inside the tree")
            unclassified += clade
            path_names.append(name), stack_drop.append(False), stack_rank.append(None)
            continue
        parent_drop = stack_drop[-1] if stack_drop else False
        parent_rank = stack_rank[-1] if stack_rank else None
        here_drop = name.lower() in drop_low or taxid.lower() in drop_low
        if here_drop:
            rep.dropped_found.update({name.lower(), taxid.lower()} & drop_low)
        if here_drop and parent_rank is not None and ({name.lower(), taxid.lower()} & strict):
            raise ValueError(f"{path}: panel lineage {name!r} lies below rank {rank}: reads assigned between the two "
                             f"could belong to the panel's organisms and be counted twice. Use a rank at or below "
                             f"the lineage, or name the enclosing rank-{rank} clade as the panel lineage.")
        if here_drop and not parent_drop:
            dropped += clade
            if parent_rank is not None:          # a dropped lineage below a kept row: take its reads out
                rep.counts[parent_rank] -= clade
        in_drop = parent_drop or here_drop
        rank_here = parent_rank
        if code == rank and not in_drop:
            if parent_rank is not None:
                raise ValueError(f"{path}: row {name!r} at rank {rank} lies inside another row at rank {rank}")
            if taxid in rep.counts:
                raise ValueError(f"{path}: taxID {taxid} occurs twice at rank {rank}")
            rep.counts[taxid] = clade
            rep.names[taxid] = name
            rep.lineages[taxid] = ";".join(path_names + [name])
            rank_here = taxid
        if not in_drop and rank_here is None:
            unresolved += own                    # assigned above the chosen rank (or in a clade without it)
        path_names.append(name), stack_drop.append(in_drop), stack_rank.append(rank_here)
    close(0)
    if not seen_any:
        raise ValueError(f"{path}: no report rows found")
    kept = sum(rep.counts.values())
    if any(v < 0 for v in rep.counts.values()):
        raise ValueError(f"{path}: a dropped lineage holds more reads than the row above it")
    residual = total - unclassified - dropped - kept - unresolved
    allowed = 0.5 * n_rows + 1 if rounds else 0
    if abs(residual) > allowed:
        raise ValueError(f"{path}: the tree does not add up ({residual} reads unaccounted for in {n_rows} rows of a "
                         f"{report_format} report); is the file complete and is column 2 the clade read count?")
    # below_threshold is 0 here: a report made with thresholds no longer says which reads they removed
    rep.ledger = {"input_fragments": total, "unclassified": unclassified, "below_threshold": 0,
                  "dropped_lineages": dropped, "kept_at_rank": kept, "unresolved_above_rank": unresolved,
                  "rounding_residual": residual, "report_rows": n_rows, "report_kind": report_format}
    assert sum(rep.ledger[c] for c in LEDGER_COLS[1:]) == rep.ledger["input_fragments"]
    return rep


def find_reports(input_dir: str | Path, suffix: str) -> dict[str, Path]:
    """Reports in a folder by file-name ending; sample = file name without the ending.

    Refuses endings that catch two reports of one sample (X.kreport and X.quant.kreport with '.kreport').
    """
    if not suffix:
        raise ValueError("give --suffix, the file-name ending of the reports to use (for example .quant.kreport)")
    files = sorted(p for p in Path(input_dir).iterdir() if p.is_file() and p.name.endswith(suffix))
    if not files:
        raise ValueError(f"{input_dir}: no file ends with {suffix!r}")
    found: dict[str, Path] = {}
    for p in files:
        name = p.name[: -len(suffix)]
        if not name:
            raise ValueError(f"{p}: the file name is only the suffix")
        if name in found:
            raise ValueError(f"duplicate sample {name!r}: {found[name].name} and {p.name}")
        found[name] = p
    names = sorted(found)
    for a in names:
        clash = [b for b in names if b != a and b.startswith(a + ".")]
        if clash:
            raise ValueError(f"--suffix {suffix!r} matches more than one report for sample {a!r}: "
                             f"{found[a].name} and {', '.join(found[b].name for b in clash)}. "
                             "Give a more specific --suffix (for example .quant.kreport) or use --reports.")
    return found


def read_report_sheet(path: str | Path) -> dict[str, Path]:
    """TSV with columns sample and report."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"sample", "report"} <= set(df.columns):
        raise ValueError(f"{path}: needs the columns 'sample' and 'report'")
    dup = df["sample"][df["sample"].duplicated()].unique()
    if len(dup):
        raise ValueError(f"{path}: duplicate sample name(s): {', '.join(dup)}")
    dupf = df["report"][df["report"].duplicated()].unique()
    if len(dupf):
        raise ValueError(f"{path}: the same report is listed for two samples: {', '.join(dupf)}")
    base = Path(path).parent
    out = {}
    for s, r in zip(df["sample"], df["report"]):
        p = Path(r) if Path(r).is_absolute() else base / r
        if not p.exists():
            raise FileNotFoundError(f"{path}: report of sample {s} not found: {p}")
        out[s] = p
    return out


def _check_namespace(namespace: str) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_.\-]+", namespace):
        raise ValueError("--namespace / --taxonomy-name may hold letters, digits, '_', '.' and '-' only")


def _write_table(parsed: dict, sources: dict, outdir: Path, rank: str, namespace: str, evidence: str) -> Path:
    samples = list(parsed)
    names, lineages = {}, {}
    for s in samples:
        for t, n in parsed[s].names.items():
            if t in names and names[t] != n:
                log.warning("taxID %s is named %r and %r in different samples; keeping the first", t, names[t], n)
            names.setdefault(t, n)
            lineages.setdefault(t, parsed[s].lineages[t])
    taxids = sorted(names, key=lambda t: -sum(parsed[s].counts.get(t, 0) for s in samples))
    if not taxids:
        raise ValueError(f"no read at rank {rank} in any sample")
    fid = {t: f"{namespace}|{t}|{rank}" for t in taxids}
    counts = pd.DataFrame({s: [int(parsed[s].counts.get(t, 0)) for t in taxids] for s in samples},
                          index=[fid[t] for t in taxids])
    counts.index.name = "feature_id"
    outdir.mkdir(parents=True, exist_ok=True)
    refuse_overlap([outdir / n for n in ("counts.tsv", "features.tsv", "ledger.tsv", manifest.NAME)],
                   list(sources.values()), "taxa-table")
    manifest.begin(outdir, "taxa-table")
    with staged(outdir / "counts.tsv") as tmp:
        counts.to_csv(tmp, sep="\t")
    feats = pd.DataFrame({
        "feature_id": [fid[t] for t in taxids], "label": [names[t] for t in taxids], "namespace": namespace,
        "accession": taxids, "identity_class": "taxon", "rank": rank, "evidence": evidence,
        "lineage": [lineages[t] for t in taxids]})
    with staged(outdir / "features.tsv") as tmp:
        feats.to_csv(tmp, sep="\t", index=False)
    ledger = pd.DataFrame([{"sample": s, "source": str(sources[s]), **parsed[s].ledger} for s in samples])
    assert (ledger[list(LEDGER_COLS[1:])].sum(axis=1) == ledger["input_fragments"]).all()
    assert (counts.sum(axis=0).to_numpy() == ledger["kept_at_rank"].to_numpy()).all()
    with staged(outdir / "ledger.tsv") as tmp:
        ledger.to_csv(tmp, sep="\t", index=False)
    manifest.seal(outdir, ["counts.tsv", "features.tsv", "ledger.tsv"], {"step": "taxa-table", "rank": rank,
                                                                         "namespace": namespace})
    log.info("taxa-table: %d features at rank %s x %d samples -> %s", len(taxids), rank, len(samples), outdir)
    return outdir / "counts.tsv"


def build_table(reports: dict[str, Path], outdir: str | Path, rank: str = "S", drop: tuple = (),
                namespace: str = "taxid", panel_lineage: tuple = (), *, report_format: str) -> Path:
    """Write counts.tsv, features.tsv and ledger.tsv from one Kraken-style report per sample.
    `panel_lineage`: lineages a panel covers (see parse_report); dropped, and refused below the rank."""
    if not reports:
        raise ValueError("no reports given")
    _check_namespace(namespace)
    check_sample_names(list(reports), "reports")
    parsed = {}
    for sample, path in reports.items():
        log.info("taxa-table: %s <- %s", sample, path)
        parsed[sample] = parse_report(path, rank, drop, panel_lineage, report_format=report_format)
    wanted = {str(d).strip().lower() for d in tuple(drop) + tuple(panel_lineage) if str(d).strip()}
    found = set().union(*(r.dropped_found for r in parsed.values()))
    missing = sorted(wanted - found)
    if missing:
        raise ValueError(f"--drop-lineage entr{'y' if len(missing) == 1 else 'ies'} not found in any report: "
                         f"{', '.join(missing)} (names must match the report exactly, case aside; a taxID also works)")
    return _write_table(parsed, reports, Path(outdir), rank, namespace, f"clade reads of a {report_format} report")


# ---------------------------------------------------------------- per-read classifier output + taxonomy

def _dmp_fields(line: str) -> list[str]:
    f = line.rstrip("\n").split("\t|\t")
    f[-1] = f[-1].rstrip("|").rstrip("\t")
    return f


class Taxonomy:
    """Parent and rank of every taxID from nodes.dmp (NCBI taxdump, or `centrifuger-inspect --taxonomy-tree`)."""

    def __init__(self, nodes: str | Path):
        self.parent: dict[str, str] = {}
        self.rank: dict[str, str] = {}
        with open_text(nodes) as fh:
            for line in fh:
                f = _dmp_fields(line)
                if len(f) < 3:
                    raise ValueError(f"{nodes}: not a nodes.dmp line (taxID | parent | rank): {line[:80]!r}")
                self.parent[f[0].strip()] = f[1].strip()
                self.rank[f[0].strip()] = f[2].strip().lower()
        if not self.parent:
            raise ValueError(f"{nodes}: empty taxonomy")
        self._path: dict[str, tuple] = {}

    def path(self, taxid: str) -> tuple:
        """taxid, its parent, ..., root. Empty when the taxID is not in the taxonomy."""
        hit = self._path.get(taxid)
        if hit is not None:
            return hit
        out, t = [], taxid
        while t in self.parent and len(out) < 300:
            out.append(t)
            p = self.parent[t]
            if p == t:
                break
            t = p
        self._path[taxid] = tuple(out)
        return self._path[taxid]

    def lca(self, taxids: list) -> str | None:
        """Lowest common ancestor; None when a taxID is unknown."""
        paths = [self.path(t) for t in taxids]
        if any(not p for p in paths):
            return None
        common = set(paths[0]).intersection(*map(set, paths[1:]))
        return next((t for t in paths[0] if t in common), None)


def scan_names(names: str | Path, want_ids: set | None = None, want_names: set | None = None) -> tuple[dict, dict]:
    """One pass over names.dmp: scientific names of `want_ids`, and the taxIDs of the (lower-case) `want_names`."""
    by_id, by_name = {}, {}
    want_ids, want_names = want_ids or set(), want_names or set()
    with open_text(names) as fh:
        for line in fh:
            f = _dmp_fields(line)
            if len(f) < 3 or f[-1].strip() != "scientific name":
                continue
            tid, name = f[0].strip(), f[1].strip()
            if tid in want_ids:
                by_id[tid] = name
            if name.lower() in want_names:
                by_name.setdefault(name.lower(), set()).add(tid)
    return by_id, by_name


def build_table_from_per_read(per_read: dict[str, Path], outdir: str | Path, nodes: str | Path, names: str | Path,
                              rank: str = "S", drop: tuple = (), namespace: str = "taxid", min_score: int = 0,
                              min_hitlen: int = 0, panel_lineage: tuple = ()) -> Path:
    """Counts at one rank straight from Centrifuge / Centrifuger per-read output (no centrifuger-quant needed).

    Every read gets one outcome (see `centrifuge_reads` for the first three):
      unclassified       no assignment
      below_threshold    assignments, none with score >= min_score and hitLength >= min_hitlen
      (classified)       one taxID, or the lowest common ancestor of the accepted taxIDs when there are several
        dropped_lineages     the taxID lies in a dropped lineage
        kept_at_rank         the taxID is at the rank or below it: counted for its ancestor at the rank
        unresolved_above_rank  the taxID is above the rank, or not in the taxonomy
    unclassified + below_threshold are exactly the reads `extract` writes with the same thresholds.

    `panel_lineage` entries are dropped too, but refused when they lie below the rank: a read whose accepted
    taxIDs straddle the lineage gets their LCA, which can sit above the lineage and inside a kept rank clade, so
    a read of the panel's organism would be counted in both tables.
    """
    if not per_read:
        raise ValueError("no per-read files given")
    if rank not in RANKS:
        raise ValueError(f"rank must be one of {', '.join(RANKS)}, not {rank!r}")
    _check_namespace(namespace)
    check_sample_names(list(per_read), "per-read files")
    min_score, min_hitlen = int(min_score or 0), int(min_hitlen or 0)
    tax = Taxonomy(nodes)
    strict = {str(d).strip().lower() for d in panel_lineage if str(d).strip()}
    wanted = {str(d).strip().lower() for d in drop if str(d).strip()} | strict
    _, name_hits = scan_names(names, want_names=wanted)
    drop_ids = set()
    target = RANK_NAMES[rank]
    for d in sorted(wanted):
        ids = name_hits.get(d, set()) | ({d} if d in tax.parent else set())
        if not ids:
            raise ValueError(f"--drop-lineage entry not found in the taxonomy: {d} (a scientific name or a taxID)")
        if d in strict:
            for t in ids:
                above = [a for a in tax.path(t)[1:] if tax.rank.get(a) in target]
                if above:
                    raise ValueError(f"panel lineage {d!r} (taxID {t}) lies below rank {rank} (inside taxID "
                                     f"{above[0]}): reads assigned between the two could belong to the panel's "
                                     "organisms and be counted twice. Use a rank at or below the lineage, or name the "
                                     f"enclosing rank-{rank} clade as the panel lineage.")
        drop_ids |= ids
    resolved: dict[str, tuple] = {}          # taxid -> (outcome, taxid at the rank)

    def resolve(taxid: str) -> tuple:
        hit = resolved.get(taxid)
        if hit is None:
            path = tax.path(taxid)
            if not path:
                hit = ("unknown", None)
            elif drop_ids.intersection(path):
                hit = ("dropped", None)
            else:
                at = next((t for t in path if tax.rank.get(t) in target), None)
                hit = ("kept", at) if at else ("unresolved", None)
            resolved[taxid] = hit
        return hit

    parsed, lca_cache = {}, {}
    for sample, path in per_read.items():
        log.info("taxa-table: %s <- %s (per-read output)", sample, path)
        rep = Report()
        n = {"total": 0, "unclassified": 0, "below": 0, "dropped": 0, "unresolved": 0, "unknown": 0, "lca": 0}
        for _, status, taxids in centrifuge_reads(path, min_score, min_hitlen):
            n["total"] += 1
            if status == 0:
                n["unclassified"] += 1
                continue
            if status == 1:
                n["below"] += 1
                continue
            if len(taxids) == 1:
                taxid = taxids[0]
            else:
                key = tuple(sorted(taxids))
                if key not in lca_cache:
                    lca_cache[key] = tax.lca(list(key))
                taxid = lca_cache[key]
                n["lca"] += 1
            outcome, at = resolve(taxid) if taxid is not None else ("unknown", None)
            if outcome == "kept":
                rep.counts[at] = rep.counts.get(at, 0) + 1
            elif outcome == "dropped":
                n["dropped"] += 1
            else:
                n["unresolved"] += 1
                n["unknown"] += outcome == "unknown"
        if n["total"] == 0:
            raise ValueError(f"{path}: no reads found")
        if n["unknown"]:
            log.warning("%s: %d reads carry a taxID that is not in %s; counted as unresolved_above_rank "
                        "(use the taxonomy the index was built with)", sample, n["unknown"], nodes)
        rep.ledger = {"input_fragments": n["total"], "unclassified": n["unclassified"], "below_threshold": n["below"],
                      "dropped_lineages": n["dropped"], "kept_at_rank": sum(rep.counts.values()),
                      "unresolved_above_rank": n["unresolved"], "rounding_residual": 0,
                      "not_in_taxonomy": n["unknown"], "reads_assigned_by_lca": n["lca"],
                      "min_score": min_score, "min_hitlen": min_hitlen}
        parsed[sample] = rep
    need = set().union(*(r.counts for r in parsed.values()))
    lineage_ids = set().union(*(tax.path(t) for t in need)) if need else set()
    id_names, _ = scan_names(names, want_ids=lineage_ids)
    for rep in parsed.values():
        for t in rep.counts:
            rep.names[t] = id_names.get(t, f"taxID {t}")
            rep.lineages[t] = ";".join(id_names.get(x, x) for x in reversed(tax.path(t)))
    return _write_table(parsed, per_read, Path(outdir), rank, namespace,
                        "reads of the per-read classifier output projected to the rank")


def panel_scope(panel_counts: str | Path) -> tuple[str, Path | None]:
    """Reads scope recorded by panel-quant next to its counts; the quant tables must be the published ones."""
    from .panel import verify_quant
    rec = verify_quant(Path(panel_counts).parent)
    summary = rec.get("input_summary")
    return rec.get("reads_scope", "all"), Path(summary) if summary else None


def check_merge_mode(panel_counts: str | Path, panel_lineage: str | None) -> list[str]:
    """Decide which lineages must leave the background so that no read is represented twice.

    scope 'all':          the panel counted every read, also those the classifier assigned to the panel's own
                          organisms -> the panel lineage must be dropped from the reports (--panel-lineage).
    scope 'unclassified': the panel counted only reads the classifier did not assign -> the two tables hold
                          disjoint reads and nothing has to be dropped.
    """
    scope, _ = panel_scope(panel_counts)
    given = [x.strip() for x in (panel_lineage or "").split(",") if x.strip()]
    if scope == "unclassified":
        return [x for x in given if x.lower() != "none"]
    if not given:
        raise ValueError("panel-quant counted ALL reads of the samples, so reads of the panel's organisms may also sit "
                         "in the classifier reports. Name the lineage the panel covers with --panel-lineage (e.g. "
                         "Viruses) so it is left out of the background, or --panel-lineage none when the classifier "
                         "database does not hold these organisms. (Counting the panel on the unclassified_reads/ of "
                         "`extract` instead makes the two tables disjoint and needs no lineage.)")
    if [x.lower() for x in given] == ["none"]:
        log.warning("--panel-lineage none: no lineage is dropped; reads of the panel's organisms that the classifier "
                    "assigned would be counted in both tables")
        return []
    if any(x.lower() == "none" for x in given):
        raise ValueError("--panel-lineage is either 'none' or a list of lineages")
    return given


def merge_with_panel(taxa_dir: str | Path, panel_counts: str | Path, panel_features: str | Path | None,
                     outdir: str | Path) -> Path:
    """Append the taxa as background features to a panel count table and reconcile the reads of both.

    Writes counts.tsv, features.tsv, merge_summary.tsv and reconciliation.tsv. The reconciliation row of a sample is
    classified kept + dropped + unresolved (+ rounding) + panel features + panel rest; for a panel counted on the
    unclassified reads this must equal the sample's total input reads, which is checked when the total is known.
    """
    taxa_dir, outdir = Path(taxa_dir), Path(outdir)
    panel_counts = Path(panel_counts)
    panel_features = Path(panel_features) if panel_features else panel_counts.parent / "features.tsv"
    if not panel_features.exists():
        raise FileNotFoundError(f"panel features table not found: {panel_features} (pass --merge-panel-features)")
    # exactly the files a verified seal covers: the panel tables and the taxa tables this merge reads
    rec_pc = manifest.verify_file(panel_counts, what="panel count tables")
    rec_pf = manifest.verify_file(panel_features, what="panel feature tables")
    rec_t = manifest.verify(taxa_dir, required=("counts.tsv", "features.tsv", "ledger.tsv"), what="taxa tables")
    scope, summary_path = panel_scope(panel_counts)
    # read private copies, each checked against its manifest, so no table can change between check and read
    wanted = {str(panel_counts): rec_pc["artifacts"][panel_counts.name],
              str(panel_features): rec_pf["artifacts"][panel_features.name],
              **{str(taxa_dir / n): rec_t["artifacts"][n] for n in ("counts.tsv", "features.tsv", "ledger.tsv")}}
    pl_path = panel_counts.parent / "ledger.tsv"
    if pl_path.exists():
        wanted[str(pl_path)] = rec_pc["artifacts"].get("ledger.tsv")
    if summary_path:
        want = rec_pc.get("input_summary_sha256")
        if not summary_path.exists() or not want:
            raise manifest.ManifestError(f"the extract summary {summary_path} recorded by panel-quant is missing or "
                                         "has no recorded checksum; run panel-quant again")
        wanted[str(summary_path)] = want              # checked like the other copies
    outdir.mkdir(parents=True, exist_ok=True)
    stage = make_staging_dir(outdir, "taxa-table --merge-panel")
    try:
        cp = manifest.private_copies(wanted, stage)
        pc = pd.read_csv(cp[str(panel_counts)], sep="\t").set_index("feature_id")
        pf = pd.read_csv(cp[str(panel_features)], sep="\t")
        tc = pd.read_csv(cp[str(taxa_dir / "counts.tsv")], sep="\t").set_index("feature_id")
        tf = pd.read_csv(cp[str(taxa_dir / "features.tsv")], sep="\t")
        tl = pd.read_csv(cp[str(taxa_dir / "ledger.tsv")], sep="\t", dtype={"sample": str},
                         keep_default_na=False).set_index("sample")
        pl = (pd.read_csv(cp[str(pl_path)], sep="\t", dtype={"sample": str}, keep_default_na=False)
              .set_index("sample") if str(pl_path) in cp else None)
        summ = (pd.read_csv(cp[str(summary_path)], sep="\t", dtype={"sample": str}, keep_default_na=False)
                .set_index("sample") if str(summary_path) in cp else None)
    finally:
        remove_staging_dir(stage)
    only_p, only_t = sorted(set(pc.columns) - set(tc.columns)), sorted(set(tc.columns) - set(pc.columns))
    if only_p or only_t:
        raise ValueError("the panel table and the reports do not hold the same samples; "
                         f"only in the panel: {', '.join(only_p) or '-'}; only in the reports: {', '.join(only_t) or '-'}")
    both = set(pc.index) & set(tc.index)
    if both:
        raise ValueError(f"feature IDs occur in both tables: {', '.join(sorted(both)[:5])}")
    samples = list(pc.columns)
    rows, problems = [], []
    for s in samples:
        t = tl.loc[s]
        classified = int(t["kept_at_rank"] + t["dropped_lineages"] + t["unresolved_above_rank"] + t["rounding_residual"])
        feat = int(pc[s].sum())
        panel_in = int(pl.loc[s, "input_fragments"]) if pl is not None and s in pl.index else None
        row = {"sample": s, "reads_scope": scope, "classified_kept_at_rank": int(t["kept_at_rank"]),
               "classified_dropped_lineages": int(t["dropped_lineages"]),
               "classified_unresolved_above_rank": int(t["unresolved_above_rank"]),
               "classified_rounding_residual": int(t["rounding_residual"]),
               "panel_features": feat, "panel_rest": None if panel_in is None else panel_in - feat,
               "panel_input_fragments": panel_in,
               "classifier_unclassified_side": int(t["unclassified"] + t["below_threshold"]),
               "total_input_reads": None, "reconstructed_total": None, "difference": None, "status": "not_checked"}
        if scope == "unclassified" and panel_in is not None:
            row["reconstructed_total"] = classified + panel_in
            if summ is not None and s in summ.index:
                total = int(summ.loc[s, "total_reads"])
                extracted = int(summ.loc[s, "unclassified_reads"])
                if panel_in != extracted:
                    problems.append(f"{s}: panel-quant counted {panel_in} fragments but `extract` wrote {extracted}")
                row["total_input_reads"] = total
            elif int(t["unclassified"] + t["below_threshold"]) > 0:
                row["total_input_reads"] = int(t["input_fragments"])
            if row["total_input_reads"] is not None:
                row["difference"] = row["reconstructed_total"] - row["total_input_reads"]
                # exact: a report's rounding is already inside `classified` as its rounding_residual
                row["status"] = "ok" if row["difference"] == 0 else "MISMATCH"
                if row["status"] == "MISMATCH":
                    problems.append(f"{s}: classified reads in the report ({classified}) + unclassified reads counted "
                                    f"on the panel ({panel_in}) = {row['reconstructed_total']}, but the sample has "
                                    f"{row['total_input_reads']} reads")
        elif scope == "all":
            row["status"] = "not_applicable (panel counted all reads)"
        rows.append(row)
    manifest.begin(outdir, "taxa-table --merge-panel")
    with staged(outdir / "reconciliation.tsv") as tmp:
        pd.DataFrame(rows).to_csv(tmp, sep="\t", index=False)
    if problems:
        raise ValueError("the classified and unclassified sides do not add up to the input reads (were the report and "
                         "`extract` made from the same classifier run and with the same --min-score/--min-hitlen?):\n  "
                         + "\n  ".join(problems) + f"\nsee {outdir / 'reconciliation.tsv'}")
    counts = pd.concat([pc, tc[samples]])
    with staged(outdir / "counts.tsv") as tmp:
        counts.to_csv(tmp, sep="\t")
    pf = pf[pf["feature_id"].isin(pc.index)].assign(source="panel")
    feats = pd.concat([pf, tf.assign(source="background")], ignore_index=True)
    with staged(outdir / "features.tsv") as tmp:
        feats.to_csv(tmp, sep="\t", index=False)
    summary = pd.DataFrame({"sample": samples, "reads_scope": scope, "panel_features": len(pc),
                            "background_features": len(tc), "panel_fragments": pc.sum(axis=0).to_numpy(),
                            "background_fragments": tc[samples].sum(axis=0).to_numpy()})
    with staged(outdir / "merge_summary.tsv") as tmp:
        summary.to_csv(tmp, sep="\t", index=False)
    manifest.seal(outdir, ["counts.tsv", "features.tsv", "merge_summary.tsv", "reconciliation.tsv"],
                  {"step": "taxa-table --merge-panel", "reads_scope": scope})
    log.info("merged %d panel + %d background features (reads scope: %s) -> %s", len(pc), len(tc), scope, outdir)
    return outdir / "counts.tsv"
