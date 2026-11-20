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

import numpy as np
import pandas as pd

from . import em
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
    em: dict = field(default_factory=dict)           # --assign em: equivalence classes and the fitted proportions
    unique: dict = field(default_factory=dict)       # --assign em: taxid -> reads whose candidates are that taxon only


# EM targets for reads that leave the table compete with the taxa for ambiguous reads. Each keeps its own taxID
# (a dropped plant and a dropped animal are different candidates), behind a prefix that no taxID can have.
T_DROPPED, T_UNRESOLVED, T_UNKNOWN = "#dropped|", "#unresolved|", "#not_in_taxonomy|"


def _outcome_of_token(tok: str) -> str:
    for prefix, name in ((T_DROPPED, "dropped"), (T_UNRESOLVED, "unresolved"), (T_UNKNOWN, "unknown")):
        if tok.startswith(prefix):
            return name
    return "kept"
ASSIGN_MODES = ("lca", "em")
DECIMALS = 6          # expected (fractional) counts are written with this many decimals


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


def _check_ledger(ledger: pd.DataFrame, counts: pd.DataFrame) -> None:
    """Every read has one outcome, and the counts are the reads kept at the rank (exactly, or to the written
    decimals for expected counts)."""
    tol = 10 ** -DECIMALS * (len(counts) + 10)
    if not np.allclose(ledger[list(LEDGER_COLS[1:])].sum(axis=1), ledger["input_fragments"], rtol=0, atol=tol):
        raise AssertionError("ledger outcomes do not add up to the input reads")
    if not np.allclose(counts.sum(axis=0).to_numpy(dtype=float), ledger["kept_at_rank"].to_numpy(dtype=float),
                       rtol=0, atol=tol):
        raise AssertionError("the counts do not add up to kept_at_rank")


def _write_table(parsed: dict, sources: dict, outdir: Path, rank: str, namespace: str, evidence: str,
                 imputed: list | None = None, record: dict | None = None, unique: bool = False) -> Path:
    """`imputed`: M lists, one per imputation, of (counts {sample: {taxid: int}}, ledgers {sample: dict}).
    `unique` (--assign em): also write unique_counts.tsv, the rows and samples of counts.tsv with each report's
    `unique` reads (whole reads whose candidate set is that taxon alone)."""
    samples = list(parsed)
    imputed = imputed or []
    fractional = any(isinstance(v, float) and not float(v).is_integer()
                     for r in parsed.values() for v in r.counts.values())
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
    value = (lambda v: round(float(v), DECIMALS)) if fractional else int
    counts = pd.DataFrame({s: [value(parsed[s].counts.get(t, 0)) for t in taxids] for s in samples},
                          index=[fid[t] for t in taxids])
    counts.index.name = "feature_id"
    outdir.mkdir(parents=True, exist_ok=True)
    imp_names = [f"imputations/{kind}_{m + 1:03d}.tsv" for m in range(len(imputed)) for kind in ("counts", "ledger")]
    if (outdir / "imputations").is_symlink():
        raise ValueError(f"{outdir / 'imputations'} is a link; the imputation tables must be written inside the "
                         "output folder itself")
    extra_names = ["unique_counts.tsv"] if unique else []
    refuse_overlap([outdir / n for n in ("counts.tsv", "features.tsv", "ledger.tsv", manifest.NAME, *extra_names,
                                         *imp_names)], list(sources.values()), "taxa-table")
    if unique:
        extra = {t for s in samples for t, v in parsed[s].unique.items() if v} - set(taxids)
        if extra:
            raise AssertionError(f"unique reads of taxa without an expected count: {sorted(extra)[:3]}")
        uc = pd.DataFrame({s: [int(parsed[s].unique.get(t, 0)) for t in taxids] for s in samples},
                          index=[fid[t] for t in taxids])
        uc.index.name = "feature_id"
        if (uc.to_numpy(dtype=float) > counts.to_numpy(dtype=float) + 10 ** -DECIMALS * 10).any():
            raise AssertionError("a unique-read count exceeds the expected count of its taxon")
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "taxa-table")
    manifest.remove_stale(outdir, previous, imp_names + extra_names)
    with staged(outdir / "counts.tsv") as tmp:
        counts.to_csv(tmp, sep="\t", float_format=f"%.{DECIMALS}f" if fractional else None)
    if unique:
        with staged(outdir / "unique_counts.tsv") as tmp:
            uc.to_csv(tmp, sep="\t")
    feats = pd.DataFrame({
        "feature_id": [fid[t] for t in taxids], "label": [names[t] for t in taxids], "namespace": namespace,
        "accession": taxids, "identity_class": "taxon", "rank": rank, "evidence": evidence,
        "lineage": [lineages[t] for t in taxids]})
    with staged(outdir / "features.tsv") as tmp:
        feats.to_csv(tmp, sep="\t", index=False)
    ledger = pd.DataFrame([{"sample": s, "source": str(sources[s]), **parsed[s].ledger} for s in samples])
    _check_ledger(ledger, counts)
    with staged(outdir / "ledger.tsv") as tmp:
        ledger.to_csv(tmp, sep="\t", index=False)
    written = ["counts.tsv", "features.tsv", "ledger.tsv", *extra_names]
    known = set(taxids)
    for m, (icounts, iledgers) in enumerate(imputed):
        extra = {t for s in samples for t in icounts[s]} - known
        if extra:
            raise AssertionError(f"imputation {m + 1}: taxa without an expected count: {sorted(extra)[:3]}")
        ic = pd.DataFrame({s: [int(icounts[s].get(t, 0)) for t in taxids] for s in samples},
                          index=[fid[t] for t in taxids])
        ic.index.name = "feature_id"
        il = pd.DataFrame([{"sample": s, "source": str(sources[s]), **iledgers[s]} for s in samples])
        _check_ledger(il, ic)
        for kind, df, idx in (("counts", ic, True), ("ledger", il, False)):
            name = f"imputations/{kind}_{m + 1:03d}.tsv"
            with staged(outdir / name) as tmp:
                df.to_csv(tmp, sep="\t", index=idx)
            written.append(name)
    manifest.seal(outdir, written, {"step": "taxa-table", "rank": rank, "namespace": namespace,
                                    **(record or {}), "imputations": len(imputed)})
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
                              min_hitlen: int = 0, panel_lineage: tuple = (), assign: str = "lca",
                              imputations: int = 0, seed: int = 1) -> Path:
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

    assign 'lca' (default): a read with several accepted taxIDs gets their lowest common ancestor.
    assign 'em': such a read is shared among its candidates (the taxa at the rank, a dropped lineage, 'above the
    rank', 'not in the taxonomy') in proportion to their maximum-likelihood abundance in the sample (see em.py);
    counts and ledger then hold expected reads with DECIMALS decimals, and `imputations` > 0 also writes that many
    integer count tables drawn from the posterior assignment probabilities (imputations/counts_NNN.tsv and
    ledger_NNN.tsv) for `da --imputations`. Reads with one candidate are counted the same way in both modes.
    With 'em', unique_counts.tsv (same rows and samples, whole reads) holds the reads whose candidate set resolves
    to that one taxon at the rank: the evidence that does not depend on how shared reads are divided. The LCA output
    is unchanged.
    """
    if assign not in ASSIGN_MODES:
        raise ValueError(f"--assign must be one of {', '.join(ASSIGN_MODES)}")
    imputations = int(imputations or 0)
    if int(seed) < 0:
        raise ValueError("--seed must be 0 or a positive whole number")
    if imputations < 0 or (imputations and assign != "em"):
        raise ValueError("--imputations needs --assign em and a positive number")
    if not per_read:
        raise ValueError("no per-read files given")
    if rank not in RANKS:
        raise ValueError(f"rank must be one of {', '.join(RANKS)}, not {rank!r}")
    _check_namespace(namespace)
    check_sample_names(list(per_read), "per-read files")
    min_score, min_hitlen = int(min_score or 0), int(min_hitlen or 0)
    nodes_sha, names_sha = manifest.sha256(nodes), manifest.sha256(names)     # the taxonomy the table is made with
    tax = Taxonomy(nodes)
    strict = {str(d).strip().lower() for d in panel_lineage if str(d).strip()}
    wanted = {str(d).strip().lower() for d in drop if str(d).strip()} | strict
    _, name_hits = scan_names(names, want_names=wanted)
    drop_ids, strict_ids = set(), set()
    target = RANK_NAMES[rank]
    for d in sorted(wanted):
        ids = name_hits.get(d, set()) | ({d} if d in tax.parent else set())
        if d in strict:
            strict_ids |= ids
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

    token_of: dict[str, str] = {}

    def token(taxid: str) -> str:
        hit = token_of.get(taxid)
        if hit is None:
            outcome, at = resolve(taxid)
            hit = token_of[taxid] = {"kept": at, "dropped": T_DROPPED + taxid, "unresolved": T_UNRESOLVED + taxid,
                                     "unknown": T_UNKNOWN + taxid}[outcome]
        return hit

    parsed, lca_cache = {}, {}
    for sample, path in per_read.items():
        log.info("taxa-table: %s <- %s (per-read output, assignment %s)", sample, path, assign)
        rep = Report()
        n = {"total": 0, "unclassified": 0, "below": 0, "dropped": 0, "unresolved": 0, "unknown": 0, "lca": 0}
        classes: dict[tuple, int] = {}
        for _, status, taxids in centrifuge_reads(path, min_score, min_hitlen):
            n["total"] += 1
            if status == 0:
                n["unclassified"] += 1
                continue
            if status == 1:
                n["below"] += 1
                continue
            if assign == "em":
                if strict_ids and any(strict_ids.intersection(tax.path(t)) for t in taxids):
                    # a read that may be the panel's organism: the panel counted all reads, so it leaves the
                    # background whole (sharing it would count part of it twice)
                    n["dropped"] += 1
                    continue
                key = tuple(sorted({token(t) for t in taxids}))
                classes[key] = classes.get(key, 0) + 1
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
        if assign == "em":
            n_dropped_whole = n["dropped"]
            # unique-read evidence: reads whose candidate set resolves to one taxon at the rank (singleton classes)
            rep.unique = {k[0]: int(c) for k, c in classes.items() if len(k) == 1 and _outcome_of_token(k[0]) == "kept"}
            fit = em.expected_counts(classes)
            if not fit.converged:
                raise ValueError(f"{sample}: the EM assignment did not converge in {fit.iterations} iterations; the "
                                 "expected counts would not be the maximum-likelihood ones (report this case)")
            dropped, unresolved, unknown = float(n["dropped"]), 0.0, 0.0
            for tok, c in fit.counts.items():
                kind = _outcome_of_token(tok)
                if kind == "dropped":
                    dropped += c
                elif kind == "unresolved":
                    unresolved += c
                elif kind == "unknown":
                    unknown += c
                elif round(c, DECIMALS) > 0:
                    rep.counts[tok] = round(c, DECIMALS)
            n["dropped"] = round(dropped, DECIMALS)
            n["unknown"] = round(unknown, DECIMALS)
            n["unresolved"] = round(unresolved + unknown, DECIMALS)
            n["lca"] = 0
            rep.em = {"classes": classes, "theta": fit.theta, "iterations": fit.iterations,
                      "dropped_whole": int(round(n_dropped_whole)),
                      "converged": fit.converged,
                      "multi": sum(k for c, k in classes.items() if len(c) > 1)}
        if n["unknown"]:
            log.warning("%s: %d reads carry a taxID that is not in %s; counted as unresolved_above_rank "
                        "(use the taxonomy the index was built with)", sample, n["unknown"], nodes)
        kept = round(sum(rep.counts.values()), DECIMALS) if assign == "em" else sum(rep.counts.values())
        residual = 0
        if assign == "em":
            residual = round(n["total"] - (n["unclassified"] + n["below"] + n["dropped"] + kept + n["unresolved"]),
                             DECIMALS)
        rep.ledger = {"input_fragments": n["total"], "unclassified": n["unclassified"], "below_threshold": n["below"],
                      "dropped_lineages": n["dropped"], "kept_at_rank": kept,
                      "unresolved_above_rank": n["unresolved"], "rounding_residual": residual,
                      "not_in_taxonomy": n["unknown"], "reads_assigned_by_lca": n["lca"],
                      "assignment": assign, "reads_shared_by_em": rep.em.get("multi", 0),
                      "em_iterations": rep.em.get("iterations", 0),
                      "min_score": min_score, "min_hitlen": min_hitlen}
        parsed[sample] = rep
    imputed = [({}, {}) for _ in range(imputations)]
    for si, (sample, rep) in enumerate(parsed.items()):
        if not imputations:
            break
        # one independent stream per sample, and a different one from panel-quant's (side 1 = taxa)
        rng = np.random.default_rng([int(seed), 1, si])
        fit = em.Fit({}, rep.em["theta"], 0, True)
        for m, draw in enumerate(em.imputations(rep.em["classes"], fit, imputations, rng)):
            icounts, iledgers = imputed[m]
            dropped, unresolved, unknown = rep.em["dropped_whole"], 0, 0
            for tok in [t for t in draw if _outcome_of_token(t) != "kept"]:
                v = draw.pop(tok)
                kind = _outcome_of_token(tok)
                dropped += v * (kind == "dropped")
                unknown += v * (kind == "unknown")
                unresolved += v * (kind in ("unresolved", "unknown"))
            icounts[sample] = draw
            led = rep.ledger
            iledgers[sample] = {"input_fragments": led["input_fragments"], "unclassified": led["unclassified"],
                                "below_threshold": led["below_threshold"], "dropped_lineages": dropped,
                                "kept_at_rank": sum(draw.values()), "unresolved_above_rank": unresolved,
                                "rounding_residual": 0, "not_in_taxonomy": unknown, "imputation": m + 1,
                                "seed": int(seed)}
    for sample, rep in parsed.items():   # a taxon drawn in an imputation but rounded to 0 expected reads keeps a row
        for icounts, _ in imputed:
            for t in icounts[sample]:
                rep.counts.setdefault(t, 0.0)
    need = set().union(*(r.counts for r in parsed.values()))
    lineage_ids = set().union(*(tax.path(t) for t in need)) if need else set()
    id_names, _ = scan_names(names, want_ids=lineage_ids)
    if (manifest.sha256(nodes), manifest.sha256(names)) != (nodes_sha, names_sha):
        raise ValueError("nodes/names changed while taxa-table was reading them; run it again")
    for rep in parsed.values():
        for t in rep.counts:
            rep.names[t] = id_names.get(t, f"taxID {t}")
            rep.lineages[t] = ";".join(id_names.get(x, x) for x in reversed(tax.path(t)))
    return _write_table(parsed, per_read, Path(outdir), rank, namespace,
                        "reads of the per-read classifier output projected to the rank"
                        + (" (ambiguous reads shared by EM)" if assign == "em" else ""), imputed,
                        {"assignment": assign, **({"seed": int(seed)} if imputations else {}),
                         "nodes_sha256": nodes_sha, "names_sha256": names_sha,
                         "dropped_lineages": sorted(wanted)}, unique=assign == "em")


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


def _num(v):
    """An integer when the value is whole, else the value rounded to DECIMALS (expected counts)."""
    v = float(v)
    return int(round(v)) if abs(v - round(v)) < 10 ** -(DECIMALS + 2) else round(v, DECIMALS)


def read_panel_taxonomy(path: str | Path) -> dict[str, str]:
    """feature_id -> taxID from a TSV with the columns feature_id and taxid (other columns are ignored; an empty
    taxid means the feature stays a feature)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = {"feature_id", "taxid"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: needs the columns feature_id and taxid (missing: {', '.join(sorted(missing))})")
    dup = df["feature_id"][df["feature_id"].duplicated()]
    if len(dup):
        raise ValueError(f"{path}: feature IDs listed twice: {', '.join(sorted(set(dup))[:5])}")
    return {f: t.strip() for f, t in zip(df["feature_id"], df["taxid"]) if t.strip()}


def collapse_plan(features: list[str], mapping: dict[str, str], nodes: str | Path, names: str | Path, rank: str,
                  namespace: str, drop: tuple = ()) -> pd.DataFrame:
    """Which panel features join which taxon row: one row per panel feature with its taxID, the taxon at the rank
    ('into', a feature ID of the background table) or the reason it stays a feature."""
    tax = Taxonomy(nodes)
    wanted = {str(d).strip().lower() for d in drop if str(d).strip()}
    _, by_name = scan_names(names, want_names=wanted)
    drop_ids = set().union(*by_name.values()) if by_name else set()
    drop_ids |= {d for d in wanted if d in tax.parent}
    target = RANK_NAMES[rank]
    rows = []
    for f in features:
        t = mapping.get(f)
        row = {"feature_id": f, "taxid": t or "", "rank_taxid": "", "into": "", "status": "no_taxid"}
        if t:
            path = tax.path(t)
            at = next((x for x in path if tax.rank.get(x) in target), None)
            if not path:
                row["status"] = "not_in_taxonomy"
            elif drop_ids.intersection(path):
                row["status"] = "dropped_lineage"
            elif at is None:
                row["status"] = "above_rank"
            else:
                row.update(rank_taxid=at, into=f"{namespace}|{at}|{rank}", status="collapsed")
        rows.append(row)
    plan = pd.DataFrame(rows, columns=["feature_id", "taxid", "rank_taxid", "into", "status"])
    if (plan["status"] == "collapsed").any():
        ids = set(plan.loc[plan["status"] == "collapsed", "rank_taxid"])
        lineage_ids = set().union(*(tax.path(t) for t in ids))
        nm, _ = scan_names(names, want_ids=lineage_ids)
        plan["label"] = [nm.get(t, f"taxID {t}") if t else "" for t in plan["rank_taxid"]]
        plan["lineage"] = [";".join(nm.get(x, x) for x in reversed(tax.path(t))) if t else ""
                           for t in plan["rank_taxid"]]
    else:
        plan["label"], plan["lineage"] = "", ""
    return plan


def apply_collapse(counts: pd.DataFrame, plan: pd.DataFrame) -> pd.DataFrame:
    """Move the counts of collapsed panel features into their taxon rows (created when absent). The reads of the
    two sides are disjoint, so a row's total is the sum of both."""
    out = counts.copy()
    for r in plan[plan["status"] == "collapsed"].itertuples():
        if r.feature_id not in out.index:
            continue
        add = out.loc[r.feature_id]
        if r.into in out.index:
            out.loc[r.into] = out.loc[r.into] + add
        else:
            out.loc[r.into] = add
        out = out.drop(index=r.feature_id)
    return out


def _whole(v):
    """The exact integer value of a count, or None when it is not a whole number (decimal parsing, no float)."""
    from decimal import Decimal, InvalidOperation
    if isinstance(v, (int, np.integer)):
        return int(v)
    try:
        d = Decimal(str(v).strip())
    except (InvalidOperation, ValueError):
        return None
    if not d.is_finite() or d != d.to_integral_value():
        return None
    return int(d)


def _format_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Counts as text: whole counts exactly as integers, expected counts with DECIMALS decimals."""
    def cell(v):
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        f = float(v)
        return str(int(f)) if f.is_integer() and abs(f) < 2 ** 53 else f"{f:.{DECIMALS}f}"
    return df.apply(lambda col: col.map(cell))


def read_draw(path) -> pd.DataFrame:
    """An imputation table: every count must be written as a whole number of reads (checked on the text, before
    any floating-point parsing could hide a fraction)."""
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).set_index("feature_id")
    bad = ~raw.apply(lambda col: col.str.fullmatch(r"[0-9]+")).to_numpy()
    if bad.any():
        i, j = map(int, np.argwhere(bad)[0])
        raise ValueError(f"{path}: {raw.columns[j]} / {raw.index[i]} holds {raw.iat[i, j]!r}, not a whole number of "
                         "reads; an imputation table holds whole reads only")
    return raw.astype("int64")


def _same_whole(a, b) -> bool:
    """Two counts that must be the same whole number of reads (no truncation, no float rounding)."""
    wa, wb = _whole(a), _whole(b)
    return wa is not None and wb is not None and wa == wb


def fractional_side(df: pd.DataFrame) -> bool:
    return any(not float(v).is_integer() for v in df.to_numpy().ravel())


def merge_with_panel(taxa_dir: str | Path, panel_counts: str | Path, panel_features: str | Path | None,
                     outdir: str | Path, panel_taxonomy: str | Path | None = None, nodes: str | Path | None = None,
                     names: str | Path | None = None, drop: tuple = ()) -> Path:
    """Append the taxa as background features to a panel count table and reconcile the reads of both.

    Writes counts.tsv, features.tsv, merge_summary.tsv and reconciliation.tsv. The reconciliation row of a sample is
    classified kept + dropped + unresolved (+ rounding) + panel features + panel rest; for a panel counted on the
    unclassified reads this must equal the sample's total input reads, which is checked when the total is known.

    Expected (fractional) counts from EM are reconciled to the written decimals. When both sides, or one, carry
    imputation tables, the merged imputations are written too (imputations/counts_NNN.tsv; a side without
    imputations contributes its one table to every imputation).

    panel_taxonomy (feature_id -> taxid, e.g. from `panel-taxonomy`): a panel feature whose taxID lies at or below
    the rank of the taxa table, outside a dropped lineage, is added to that taxon's row instead of standing as a
    second feature for the same organism (needs the nodes/names of the taxa table). collapsed.tsv lists every panel
    feature with its outcome.
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
    if panel_taxonomy and scope != "unclassified":       # before anything in the output folder is touched
        raise ValueError("--panel-taxonomy needs a panel counted on the unclassified reads (panel-quant on the "
                         "unclassified_reads/ of `run`/`extract`): with all reads counted, a feature added to a taxon "
                         "row would count reads the classifier already put there")
    # read private copies, each checked against its manifest, so no table can change between check and read
    wanted = {str(panel_counts): rec_pc["artifacts"][panel_counts.name],
              str(panel_features): rec_pf["artifacts"][panel_features.name],
              **{str(taxa_dir / n): rec_t["artifacts"][n] for n in ("counts.tsv", "features.tsv", "ledger.tsv")}}
    pl_path = panel_counts.parent / "ledger.tsv"
    if pl_path.exists():
        wanted[str(pl_path)] = rec_pc["artifacts"].get("ledger.tsv")
    m_t, m_p = int(rec_t.get("imputations") or 0), int(rec_pc.get("imputations") or 0)
    if m_t and m_p and m_t != m_p:
        raise ValueError(f"the taxa table has {m_t} imputations and the panel table {m_p}; use the same number")
    n_imp = max(m_t, m_p)
    if n_imp:
        # a side without imputations may only be one whose counts carry no assignment uncertainty
        if not m_t and rec_t.get("assignment") == "em":
            raise ValueError("the panel table has imputations but the taxa table was made with --assign em and none: "
                             "rerun taxa-table with --imputations " + str(n_imp))
        if not m_p and rec_pc.get("ambiguous") == "em":
            raise ValueError("the taxa table has imputations but panel-quant was run with --ambiguous em and none: "
                             "rerun panel-quant with --imputations " + str(n_imp))
    imp_t = [taxa_dir / f"imputations/counts_{m + 1:03d}.tsv" for m in range(m_t)]
    imp_p = [panel_counts.parent / f"imputations/counts_{m + 1:03d}.tsv" for m in range(m_p)]
    led_t = [taxa_dir / f"imputations/ledger_{m + 1:03d}.tsv" for m in range(m_t)]
    led_p = [panel_counts.parent / f"imputations/ledger_{m + 1:03d}.tsv" for m in range(m_p)]
    for f, rec, base in ([(f, rec_t, taxa_dir) for f in imp_t + led_t]
                         + [(f, rec_pc, panel_counts.parent) for f in imp_p + led_p]):
        key = str(f.relative_to(base))
        if key not in rec["artifacts"]:
            raise manifest.ManifestError(f"{f} is not covered by the manifest of {base}; run that step again")
        wanted[str(f)] = rec["artifacts"][key]
    taxonomy_rec_sha = None
    if panel_taxonomy:
        if not (nodes and names):
            raise ValueError("--panel-taxonomy needs --nodes and --names (the taxonomy of the taxa table)")
        given, real = Path(panel_taxonomy), Path(panel_taxonomy).resolve()
        sealed_at = next((x for x in (real, given) if (x.parent / manifest.NAME).exists()
                          or manifest.is_generated(x.parent)), None)        # the canonical file's own seal first
        if sealed_at is not None:
            rec_x = manifest.verify_file(sealed_at, what="panel taxonomy")       # made by panel-taxonomy
            for kind, path in (("nodes", nodes), ("names", names)):
                want = rec_x.get(f"{kind}_sha256")
                if want and manifest.sha256(path) != want:
                    raise ValueError(f"{panel_taxonomy} was made with another {kind} file than {path}; use the "
                                     "taxonomy of the taxa table for both")
            if rec_x.get("nodes_sha256") and rec_x.get("names_sha256"):
                taxonomy_rec_sha = {str(nodes): rec_x["nodes_sha256"], str(names): rec_x["names_sha256"]}
            want_panel = (rec_pc.get("panel") or {}).get("panel_id")
            if rec_x.get("panel_id") and want_panel and rec_x["panel_id"] != want_panel:
                raise ValueError(f"{panel_taxonomy} was made for another panel (panel_id {rec_x['panel_id']}) than "
                                 f"the panel-quant counts (panel_id {want_panel}); run panel-taxonomy on that panel")
            wanted[str(panel_taxonomy)] = rec_x["artifacts"][sealed_at.name]
        else:
            log.warning("%s has no manifest (not made by `panel-taxonomy`); it is used as it is", panel_taxonomy)
            wanted[str(panel_taxonomy)] = None
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
        tl = pd.read_csv(cp[str(taxa_dir / "ledger.tsv")], sep="\t", dtype=str,
                         keep_default_na=False).set_index("sample")
        pl = (pd.read_csv(cp[str(pl_path)], sep="\t", dtype={"sample": str}, keep_default_na=False)
              .set_index("sample") if str(pl_path) in cp else None)
        summ = (pd.read_csv(cp[str(summary_path)], sep="\t", dtype={"sample": str}, keep_default_na=False)
                .set_index("sample") if str(summary_path) in cp else None)
        it = [read_draw(cp[str(f)]) for f in imp_t]
        ip = [read_draw(cp[str(f)]) for f in imp_p]
        itl = [pd.read_csv(cp[str(f)], sep="\t", dtype=str, keep_default_na=False).set_index("sample") for f in led_t]
        ipl = [pd.read_csv(cp[str(f)], sep="\t", dtype=str, keep_default_na=False).set_index("sample") for f in led_p]
        mapping = read_panel_taxonomy(cp[str(panel_taxonomy)]) if panel_taxonomy else None
        taxonomy_sha = manifest.sha256(cp[str(panel_taxonomy)]) if panel_taxonomy else None
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
    # whole reads must reconcile exactly; expected counts to well below one read (the taxa ledger carries its own
    # rounding residual, so only floating-point error remains)
    tol = 1e-4 if (fractional_side(pc) or fractional_side(tc[[c for c in tc.columns]]) or
                   any(not float(v).is_integer() for v in tl[["kept_at_rank", "dropped_lineages",
                                                              "unresolved_above_rank"]].to_numpy().ravel())) else 0
    rows, problems = [], []
    for s in samples:
        t = tl.loc[s]
        parts = [t["kept_at_rank"], t["dropped_lineages"], t["unresolved_above_rank"], t["rounding_residual"]]
        # whole reads in exact integer arithmetic (per sample); expected counts as floats
        whole = [_whole(v) for v in parts]
        tol_s = 0 if all(w is not None for w in whole) and not fractional_side(pc[[s]]) else 1e-4
        from decimal import Decimal
        classified = sum(whole) if tol_s == 0 else sum(Decimal(str(v)) for v in parts)
        feat = (sum(int(v) for v in pc[s]) if not fractional_side(pc[[s]])
                else _num(sum(__import__("decimal").Decimal(str(v)) for v in pc[s])))   # exact for whole reads
        panel_in = int(pl.loc[s, "input_fragments"]) if pl is not None and s in pl.index else None
        exact = lambda v: _whole(v) if _whole(v) is not None else _num(v)      # whole reads stay exact integers
        row = {"sample": s, "reads_scope": scope, "classified_kept_at_rank": exact(t["kept_at_rank"]),
               "classified_dropped_lineages": exact(t["dropped_lineages"]),
               "classified_unresolved_above_rank": exact(t["unresolved_above_rank"]),
               "classified_rounding_residual": exact(t["rounding_residual"]),
               "panel_features": feat,
               "panel_rest": None if panel_in is None else (panel_in - feat if isinstance(feat, int)
                                                            else _num(panel_in - feat)),
               "panel_input_fragments": panel_in,
               "classifier_unclassified_side": _whole(t["unclassified"]) + _whole(t["below_threshold"]),
               "total_input_reads": None, "reconstructed_total": None, "difference": None, "status": "not_checked"}
        if scope == "unclassified" and panel_in is not None:
            row["reconstructed_total"] = classified + panel_in      # int, or an exact Decimal
            if summ is not None and s in summ.index:
                total = int(summ.loc[s, "total_reads"])
                extracted = int(summ.loc[s, "unclassified_reads"])
                if panel_in != extracted:
                    problems.append(f"{s}: panel-quant counted {panel_in} fragments but `extract` wrote {extracted}")
                row["total_input_reads"] = total
            elif _whole(t["unclassified"]) + _whole(t["below_threshold"]) > 0:
                row["total_input_reads"] = _whole(t["input_fragments"])
            if row["total_input_reads"] is not None:
                diff = row["reconstructed_total"] - row["total_input_reads"]
                row["difference"] = diff
                # exact for whole reads (a report's rounding is already inside `classified` as its
                # rounding_residual); to the written decimals for expected counts
                row["status"] = "ok" if abs(row["difference"]) <= tol_s else "MISMATCH"
                if row["status"] == "MISMATCH":
                    problems.append(f"{s}: classified reads in the report ({classified}) + unclassified reads counted "
                                    f"on the panel ({panel_in}) = {row['reconstructed_total']}, but the sample has "
                                    f"{row['total_input_reads']} reads")
        elif scope == "all":
            row["status"] = "not_applicable (panel counted all reads)"
        for key in ("reconstructed_total", "difference"):       # exact while checked; written as numbers
            if row[key] is not None and not isinstance(row[key], int):
                row[key] = _num(float(row[key]))
        rows.append(row)
    # everything is computed and checked first; the output folder is touched only when all of it holds
    fractional = any(not float(v).is_integer() for v in pd.concat([pc, tc[samples]]).to_numpy().ravel())
    ffmt = f"%.{DECIMALS}f" if fractional else None
    counts = pd.concat([pc.astype(object), tc[samples].astype(object)])     # no promotion of whole counts to float
    pf = pf[pf["feature_id"].isin(pc.index)].assign(source="panel")
    feats = pd.concat([pf, tf.assign(source="background")], ignore_index=True)
    written = ["counts.tsv", "features.tsv", "merge_summary.tsv", "reconciliation.tsv"]
    plan, collapsed_table = None, None
    if mapping is not None:
        rec_rank, rec_ns = rec_t.get("rank"), rec_t.get("namespace")
        if not (rec_rank and rec_ns):
            raise ValueError(f"{taxa_dir}/MANIFEST.json records no rank/namespace; run taxa-table again")
        if not rec_t.get("nodes_sha256"):
            raise ValueError(f"the taxa table in {taxa_dir} records no taxonomy (made from reports); --panel-taxonomy "
                             "needs a table made with --per-read and the same --nodes/--names")
        unknown = sorted(set(mapping) - set(pc.index))
        if unknown:
            log.warning("--panel-taxonomy lists %d feature(s) not in the panel table (ignored): %s", len(unknown),
                        ", ".join(unknown[:3]))
        tax_sha = {str(nodes): manifest.sha256(nodes), str(names): manifest.sha256(names)}
        if (rec_t["nodes_sha256"], rec_t.get("names_sha256")) != (tax_sha[str(nodes)], tax_sha[str(names)]):
            raise ValueError(f"the taxa table in {taxa_dir} was made with another taxonomy than {nodes} / {names}; "
                             "the panel features would join rows of a different taxonomy")
        # the lineages the taxa table left out stay out, whatever the merge call passes
        plan = collapse_plan(list(pc.index), mapping, nodes, names, rec_rank, rec_ns,
                             tuple(drop) + tuple(rec_t.get("dropped_lineages") or ()))
        if any(manifest.sha256(f) != h for f, h in tax_sha.items()):
            raise ValueError("nodes/names changed while the merge was reading them; run it again")
        if taxonomy_rec_sha and taxonomy_rec_sha != tax_sha:
            raise ValueError(f"{panel_taxonomy} was made with another taxonomy than {nodes} / {names}")
        counts = apply_collapse(counts, plan)
        gone = set(plan.loc[plan["status"] == "collapsed", "feature_id"])
        new_rows = plan[(plan["status"] == "collapsed") & ~plan["into"].isin(feats["feature_id"])].drop_duplicates("into")
        add = pd.DataFrame({"feature_id": new_rows["into"], "label": new_rows["label"], "namespace": rec_ns,
                            "accession": new_rows["rank_taxid"], "identity_class": "taxon", "rank": rec_rank,
                            "evidence": "panel features assigned to this taxon (no classified read)",
                            "lineage": new_rows["lineage"], "source": "background"})
        feats = pd.concat([feats[~feats["feature_id"].isin(gone)], add], ignore_index=True)
        joined = plan.groupby("into")["feature_id"].apply(lambda x: ";".join(sorted(x))).to_dict()
        feats["collapsed_panel_features"] = feats["feature_id"].map(joined).fillna("")
        clash = sorted(set(samples) & set(plan.columns))
        if clash:
            raise ValueError(f"sample name(s) {', '.join(clash)} are also column names of collapsed.tsv; rename the "
                             "sample(s)")
        collapsed_table = pd.concat([plan.reset_index(drop=True), pc.reindex(plan["feature_id"]).reset_index(drop=True)],
                                    axis=1)
        written.append("collapsed.tsv")
    if set(feats["feature_id"]) != set(counts.index) or feats["feature_id"].duplicated().any():
        missing = sorted(set(counts.index) - set(feats["feature_id"]))
        extra = sorted(set(feats["feature_id"]) - set(counts.index))
        raise ValueError("the feature table does not describe exactly the counted features "
                         f"(counted but not described: {', '.join(missing[:3]) or '-'}; described but not counted: "
                         f"{', '.join(extra[:3]) or '-'}); use the features.tsv of the same panel-quant run")
    order = list(feats["feature_id"])
    counts = counts.loc[order]
    draws = []
    for m in range(n_imp):
        a = ip[m] if ip else pc
        b = it[m] if it else tc
        if set(a.index) != set(pc.index) or set(b.index) != set(tc.index) or a.index.duplicated().any() \
                or b.index.duplicated().any():
            raise ValueError(f"imputation {m + 1}: its rows are not the rows of the expected tables; rerun the steps "
                             "that made the imputations")
        for side, tab in (("panel", a), ("taxa", b)):
            vals = tab[samples].to_numpy()
            if not (np.isfinite(vals).all() and (vals >= 0).all() and (vals == np.round(vals)).all()):
                raise ValueError(f"imputation {m + 1}: the {side} draw holds negative or fractional counts")
        for s in samples:          # every draw keeps the reads its own ledger says it has
            if it and not _same_whole(b[s].sum(), itl[m].loc[s, "kept_at_rank"]):
                raise ValueError(f"imputation {m + 1}, {s}: the taxa draw holds {b[s].sum()} reads but its "
                                 f"ledger {itl[m].loc[s, 'kept_at_rank']}")
            if ip and not _same_whole(a[s].sum(), ipl[m].loc[s, "assigned_to_features"]):
                raise ValueError(f"imputation {m + 1}, {s}: the panel draw holds {a[s].sum()} fragments but "
                                 f"its ledger {ipl[m].loc[s, 'assigned_to_features']}")
        mc = pd.concat([a[samples], b[samples]])
        if plan is not None:
            mc = apply_collapse(mc, plan)
        if set(mc.index) != set(order) or mc.index.duplicated().any():
            raise ValueError(f"imputation {m + 1}: its rows are not the rows of the merged table; rerun the steps "
                             "that made the imputations")
        mc = mc.loc[order]
        if fractional_side(mc):
            raise AssertionError("an imputation table holds expected (fractional) counts")
        draws.append((f"imputations/counts_{m + 1:03d}.tsv", mc.astype("int64")))
    if (outdir / "imputations").is_symlink():
        raise ValueError(f"{outdir / 'imputations'} is a link; the imputation tables must be written inside the "
                         "output folder itself")
    previous = manifest.previous_artifacts(outdir)
    manifest.begin(outdir, "taxa-table --merge-panel")
    with staged(outdir / "reconciliation.tsv") as tmp:
        pd.DataFrame(rows).to_csv(tmp, sep="\t", index=False)
    if problems:
        raise ValueError("the classified and unclassified sides do not add up to the input reads (were the report and "
                         "`extract` made from the same classifier run and with the same --min-score/--min-hitlen?):\n  "
                         + "\n  ".join(problems) + f"\nsee {outdir / 'reconciliation.tsv'}")
    if collapsed_table is not None:
        with staged(outdir / "collapsed.tsv") as tmp:
            collapsed_table.to_csv(tmp, sep="\t", index=False, float_format=ffmt)
        log.info("panel taxonomy: %d of %d panel features joined a taxon row (%s)",
                 int((plan["status"] == "collapsed").sum()), len(plan),
                 ", ".join(f"{k} {v}" for k, v in plan["status"].value_counts().items()))
    with staged(outdir / "counts.tsv") as tmp:
        _format_counts(counts).to_csv(tmp, sep="\t")
    with staged(outdir / "features.tsv") as tmp:
        feats.to_csv(tmp, sep="\t", index=False)
    for name, mc in draws:
        with staged(outdir / name) as tmp:
            mc.to_csv(tmp, sep="\t")
        written.append(name)
    manifest.remove_stale(outdir, previous, written)
    summary = pd.DataFrame({"sample": samples, "reads_scope": scope, "panel_features": len(pc),
                            "background_features": len(tc), "panel_fragments": pc.sum(axis=0).to_numpy(),
                            "background_fragments": tc[samples].sum(axis=0).to_numpy()})
    with staged(outdir / "merge_summary.tsv") as tmp:
        summary.to_csv(tmp, sep="\t", index=False)
    manifest.seal(outdir, written,
                  {"step": "taxa-table --merge-panel", "reads_scope": scope, "imputations": n_imp,
                   **({"panel_taxonomy": str(Path(panel_taxonomy).resolve()), "panel_taxonomy_sha256": taxonomy_sha}
                      if panel_taxonomy else {})})
    log.info("merged %d panel + %d background features (reads scope: %s) -> %s", len(pc), len(tc), scope, outdir)
    return outdir / "counts.tsv"
