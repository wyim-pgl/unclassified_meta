"""Configuration defaults, YAML loading, sample sheets and input-folder detection."""
from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from pathlib import Path

import yaml

from .io_utils import open_text

CLASSIFIERS = ("kraken2", "centrifuge", "centrifuger", "none")
LAYERS = ("rrna", "organelle", "host", "virus", "kraken2", "protein")

DEFAULTS: dict = {
    "mode": "rna",                 # rna | dna
    "outdir": "results",
    "threads": 8,
    "memory_gb": 32,               # memory available to the assembler
    "subsample": 1_000_000,        # read pairs (or single reads) per sample; "all" keeps every read
    "seed": 11,
    "min_contig_len": 300,
    "min_score": 0,                # Centrifuge(r) only: assignments with a lower score count as unclassified
    "min_hitlen": 0,               # Centrifuge(r) only: assignments with a shorter hit count as unclassified
    "assembler": "auto",           # auto (rna: rnaspades, dna: megahit) | rnaspades | megahit | metaspades
    "references": {
        # REQUIRED (made by `unclassified-meta download-refs --outdir refs`; point refs_dir at that folder)
        "refs_dir": None,
        "rrna": None,              # plant rRNA FASTA (SILVA 138.2 Chloroplastida/Chloroplast/Mitochondria + Rfam 5.8S)
        "virus": None,             # RefSeq viral genomes FASTA (BLASTN)
        "virus_protein": None,     # RefSeq viral proteins, DIAMOND .dmnd or protein FASTA (blastx)
        "protein": None,           # Swiss-Prot Viridiplantae, DIAMOND .dmnd or protein FASTA (blastx)
        "host": [],                # REQUIRED, list of {name: ..., fasta: ...}, checked in list order
        # OPTIONAL
        "organelle": None,         # plastid / mitochondrial genomes FASTA (BLASTN)
        "extra_virus_fasta": None, # small custom virus genome FASTA (BLASTN, and TBLASTX when enabled)
        "host_proteome": None,     # host protein FASTA, added to the protein layer
        "kraken2_db": None,        # Kraken2 database directory (classifies contigs)
        "kraken2_taxonomy": None,  # ktaxonomy.tsv; default <kraken2_db>/ktaxonomy.tsv
    },
    "priority": list(LAYERS),
    "thresholds": {
        "rrna": {"min_identity": 85.0, "min_cover": 0.5, "evalue": 1e-20},
        "organelle": {"min_identity": 85.0, "min_cover": 0.5, "evalue": 1e-10},
        "host": {"min_cover": 0.5, "preset": "asm20"},
        "virus": {"blastn_min_len": 100, "evalue": 1e-10, "tblastx_extra": False},
        "protein": {"evalue": 1e-10},
    },
    "readback_preset": "--very-fast",
}


@dataclass
class Sample:
    name: str
    r1: Path
    r2: Path | None = None
    classifier_output: Path | None = None
    classifier: str = "none"

    @property
    def paired(self) -> bool:
        return self.r2 is not None


def _merge(base: dict, extra: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (extra or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def parse_subsample(value) -> int | None:
    """'all'/None/0 -> None (keep every read); otherwise a positive integer."""
    if value is None or (isinstance(value, str) and value.strip().lower() in ("all", "none", "0", "")):
        return None
    n = int(float(value))
    if n < 0:
        raise ValueError("subsample must be a positive number of reads or 'all'")
    return n or None


def load_config(path: str | Path | None, overrides: dict | None = None) -> dict:
    cfg = copy.deepcopy(DEFAULTS)
    if path:
        with open(path) as fh:
            cfg = _merge(cfg, yaml.safe_load(fh) or {})
    cfg = _merge(cfg, {k: v for k, v in (overrides or {}).items() if v is not None})
    if cfg["mode"] not in ("rna", "dna"):
        raise ValueError(f"mode must be rna or dna, not {cfg['mode']!r}")
    cfg["subsample"] = parse_subsample(cfg["subsample"])
    bad = [x for x in cfg["priority"] if x not in LAYERS]
    if bad:
        raise ValueError(f"unknown priority layer(s) {bad}; allowed: {LAYERS}")
    apply_refs_dir(cfg)
    hosts = cfg["references"].get("host") or []
    if isinstance(hosts, (str, dict)):
        hosts = [hosts]
    cfg["references"]["host"] = [h if isinstance(h, dict) else {"name": Path(h).stem, "fasta": h}
                                 for h in hosts]
    return cfg


REFS_LAYOUT = {
    "rrna": "rrna/plant_rrna.fasta",
    "virus": "viral/refseq_viral.fna",
    "virus_protein": "viral/refseq_viral_protein.dmnd",
    "protein": "protein/sprot_viridiplantae.dmnd",
}
REQUIRED_REFS = ("rrna", "virus", "virus_protein", "protein")


def apply_refs_dir(cfg: dict) -> None:
    """Fill unset required references from the download-refs folder layout."""
    refs = cfg["references"]
    d = refs.get("refs_dir")
    if not d:
        return
    for key, rel in REFS_LAYOUT.items():
        if not refs.get(key):
            path = Path(d) / rel
            if path.suffix == ".dmnd" and not path.exists() and path.with_suffix(".faa").exists():
                path = path.with_suffix(".faa")   # protein FASTA without a prebuilt DIAMOND database
            refs[key] = str(path)


def check_references(cfg: dict) -> None:
    """Fail early with one message listing every missing required reference."""
    refs, problems = cfg["references"], []
    for key in REQUIRED_REFS:
        v = refs.get(key)
        if not v:
            problems.append(f"references.{key} is not set")
        elif not Path(v).exists():
            problems.append(f"references.{key}: file not found: {v}")
    if not refs.get("host"):
        problems.append("references.host: at least one host genome is required (--host NAME=genome.fa)")
    for h in refs.get("host") or []:
        if not Path(h["fasta"]).exists():
            problems.append(f"host genome {h['name']}: file not found: {h['fasta']}")
    for key in ("organelle", "extra_virus_fasta", "host_proteome"):
        if refs.get(key) and not Path(refs[key]).exists():
            problems.append(f"references.{key}: file not found: {refs[key]}")
    if refs.get("kraken2_db") and not Path(refs["kraken2_db"]).is_dir():
        problems.append(f"references.kraken2_db: directory not found: {refs['kraken2_db']}")
    if problems:
        raise ValueError("missing references:\n  - " + "\n  - ".join(problems) +
                         "\nRun `unclassified-meta download-refs --outdir refs` first and pass --refs refs "
                         "(or set the paths in the config file).")


def read_samplesheet(path: str | Path) -> list[Sample]:
    """TSV with header: sample, r1, [r2], [classifier_output], [classifier]."""
    path = Path(path)
    lines = [l.rstrip("\n") for l in open(path) if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    if "sample" not in header or "r1" not in header:
        raise ValueError(f"{path}: header needs at least 'sample' and 'r1' columns, got {header}")
    base = path.parent
    samples, seen = [], set()
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))

        def p(key):
            v = (row.get(key) or "").strip()
            if not v or v in (".", "NA", "-"):
                return None
            q = Path(v)
            return q if q.is_absolute() else (base / q)
        clf = (row.get("classifier") or "").strip().lower() or ("kraken2" if p("classifier_output") else "none")
        if clf not in CLASSIFIERS:
            raise ValueError(f"{path}: classifier {clf!r} not in {CLASSIFIERS}")
        s = Sample(row["sample"].strip(), p("r1"), p("r2"), p("classifier_output"), clf)
        if s.classifier != "none" and s.classifier_output is None:
            raise ValueError(f"{path}: sample {s.name} has classifier {clf} but no classifier_output")
        if s.name in seen:
            raise ValueError(f"{path}: duplicate sample name {s.name}")
        seen.add(s.name)
        samples.append(s)
    validate_samples(samples)
    return samples


# column names of the output tables; a sample with one of these names would overwrite a column
RESERVED_SAMPLE_NAMES = frozenset({
    "feature_id", "label", "namespace", "accession", "identity_class", "n_sequences", "total_length", "evidence",
    "n_contigs", "weighted_identity", "identity_min", "identity_max", "mixed_identity", "reference_included",
    "best_hit_title", "parent_feature", "rank", "lineage", "source", "contig", "sample", "category"})


def check_sample_names(names, where: str = "samples") -> None:
    """Refuse duplicate sample names and names that are column names of the output tables."""
    names = [str(n) for n in names]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError(f"{where}: duplicate sample name(s): {', '.join(dup)}")
    bad = sorted(n for n in names if n.lower() in RESERVED_SAMPLE_NAMES)
    if bad:
        raise ValueError(f"{where}: sample name(s) {', '.join(bad)} are reserved (column names of the output "
                         "tables); rename the sample(s)")


# fixed columns of tables/reference_by_sample.tsv; a sample adds the column pct_<sample>
REFERENCE_TABLE_COLUMNS = frozenset({
    "category", "best_hit", "best_hit_title", "n_contigs", "reads_pooled", "pct_of_reads_used_pooled",
    "blastn_identity", "tblastx_identity", "diamond_identity", "minimap2_identity"})


def output_read_names(s: "Sample") -> list[str]:
    """File names `extract` writes for a sample under unclassified_reads/."""
    return [f"{s.name}_R1.fastq.gz", f"{s.name}_R2.fastq.gz"] if s.paired else [f"{s.name}.fastq.gz"]


def check_generated_names(samples: list["Sample"]) -> None:
    """Names made from sample names must not collide: output read files (paired A -> A_R1.fastq.gz, single
    A_R1 -> A_R1.fastq.gz) and the per-sample columns of the reference table (pct_<sample>)."""
    owner: dict[str, str] = {}
    for s in samples:
        for f in output_read_names(s):
            if f.lower() in owner:
                raise ValueError(f"samples {owner[f.lower()]} and {s.name} would write the same output file "
                                 f"unclassified_reads/{f}; rename one of them")
            owner[f.lower()] = s.name
        col = f"pct_{s.name}"
        if col in REFERENCE_TABLE_COLUMNS:
            raise ValueError(f"sample {s.name}: its column {col} would overwrite a column of the reference table; "
                             "rename the sample")


def validate_samples(samples: list[Sample]) -> None:
    if not samples:
        raise ValueError("no samples given")
    for s in samples:
        for f in (s.r1, s.r2, s.classifier_output):
            if f is not None and not Path(f).exists():
                raise FileNotFoundError(f"sample {s.name}: file not found: {f}")
        if not re.fullmatch(r"[A-Za-z0-9._-]+", s.name):
            raise ValueError(f"sample name {s.name!r}: use letters, digits, '.', '_' or '-' only")
    check_sample_names([s.name for s in samples])
    check_generated_names(samples)


# paired read files; the sample is the part before the mate token. Tried in this order:
PAIRED_RES = (
    # STAR --outReadsUnmapped Fastx: <sample>.Unmapped.out.mate1 / .mate2 (optionally gzipped)
    re.compile(r"^(?P<sample>.+)\.Unmapped\.out\.mate(?P<mate>[12])(?:\.gz)?$"),
    # <sample>[_clean]_R1.fastq.gz, <sample>_1.fq.gz, <sample>.R2.fastq, <sample>_R1_001.fastq.gz ...
    re.compile(r"^(?P<sample>.+?)(?:[._-]clean)?[._-]R?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$"),
)
SINGLE_RE = re.compile(r"^(?P<sample>.+?)\.(?:fastq|fq)(?:\.gz)?$")
CLASSIFIER_EXT = {
    ".kraken": "kraken2", ".kraken2": "kraken2", ".kraken.out": "kraken2", ".kraken2.out": "kraken2",
    ".k2": "kraken2", ".centrifuge": "centrifuge", ".centrifuge.out": "centrifuge", ".centrifuge.tsv": "centrifuge",
    ".cfr": "centrifuger", ".cfr.tsv": "centrifuger", ".centrifuger": "centrifuger", ".centrifuger.out": "centrifuger",
    ".centrifuger.tsv": "centrifuger",
}


def _is_int(x: str) -> bool:
    return re.fullmatch(r"-?\d+", x.strip()) is not None


def _is_num(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def _line_kind(f: list[str]) -> str | None:
    """What one tab-separated line looks like: a Kraken2 per-read line, a Centrifuge(r) per-read line, or a
    Kraken-style report row."""
    if len(f) >= 5 and f[0] in ("C", "U") and (_is_int(f[2]) or "(taxid" in f[2]):
        return "kraken2"
    if len(f) >= 8 and _is_int(f[2]) and all(_is_num(x) for x in f[3:8]):
        return "centrifuge"
    if len(f) in (6, 8) and _is_num(f[0]) and _is_int(f[1]) and _is_int(f[2]) \
            and re.fullmatch(r"[A-Z]\d*|-", f[-3].strip()):
        return "report"
    return None


def sniff_classifier(path: str | Path, n_lines: int = 20) -> str | None:
    """Kind of a classifier file from its content (gzip or plain), never from its name: 'kraken2', 'centrifuge'
    (also Centrifuger), 'report', or None when the first lines do not all agree on one kind."""
    kinds, n = set(), 0
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if n == 0 and f[0] == "readID":
                continue                              # Centrifuge(r) header
            kinds.add(_line_kind(f))
            n += 1
            if n >= n_lines:
                break
    if n == 0:
        return "empty"                                # no records (a header at most): valid output of no reads
    return kinds.pop() if len(kinds) == 1 and None not in kinds else None


def _classifier_ext(name: str) -> tuple[str, str] | None:
    """(extension incl. a trailing .gz, classifier) of a classifier file name, longest extension first."""
    gz = name.endswith(".gz")
    base = name[:-3] if gz else name
    for ext in sorted(CLASSIFIER_EXT, key=len, reverse=True):
        if base.endswith(ext) and len(base) > len(ext):
            return ext + (".gz" if gz else ""), CLASSIFIER_EXT[ext]
    return None


def detect_samples(folder: str | Path) -> list[Sample]:
    """Find read files and their per-read classifier outputs in one folder.

    Reads: paired files (`<s>_R1/_R2`, `<s>_1/_2`, `<s>_clean_R1/_R2`, STAR `<s>.Unmapped.out.mate1/2`) or single
    files (`<s>.fastq.gz`). Classifier outputs (`<s>.kraken2`, `<s>.cfr`, `<s>.cfr.tsv`, ... see CLASSIFIER_EXT) are
    matched to a sample by the exact sample name. A sample with two candidate outputs, or an output whose name
    matches no sample, is an error rather than a silent choice. Kraken-style reports are ignored.
    """
    folder = Path(folder)
    files = sorted(p for p in folder.iterdir() if p.is_file())
    pairs: dict[str, dict] = {}
    singles: dict[str, Path] = {}
    outputs: list[tuple[str, Path, str]] = []
    for p in files:
        hit = _classifier_ext(p.name)
        if hit:
            kind = sniff_classifier(p)
            if kind == "report":
                continue                              # a Kraken-style report, not per-read output: not used here
            fam = "kraken2" if hit[1] == "kraken2" else "centrifuge"
            if kind is None:
                raise ValueError(f"{folder}: {p.name} is named like {hit[1]} output, but its lines are not per-read "
                                 "classifier output (Kraken2: C/U, read ID, taxID ...; Centrifuge(r): read ID, seq ID, "
                                 "taxID, score, 2nd score, hit length, query length, matches)")
            if kind not in (fam, "empty"):
                raise ValueError(f"{folder}: {p.name} is named like {hit[1]} output but holds {kind} output")
            outputs.append((p.name[: -len(hit[0])], p, hit[1]))
            continue
        for rx in PAIRED_RES:
            m = rx.match(p.name)
            if m:
                mates = pairs.setdefault(m["sample"], {})
                if m["mate"] in mates:
                    raise ValueError(f"{folder}: two mate-{m['mate']} files for sample {m['sample']}: "
                                     f"{mates[m['mate']].name} and {p.name}")
                mates[m["mate"]] = p
                break
        else:
            m = SINGLE_RE.match(p.name)
            if m:
                if m["sample"] in singles:
                    raise ValueError(f"{folder}: two read files for sample {m['sample']}: {singles[m['sample']].name} "
                                     f"and {p.name}; keep one")
                singles[m["sample"]] = p
    samples: dict[str, Sample] = {}
    for name, mates in pairs.items():
        if "1" not in mates:
            raise ValueError(f"{folder}: sample {name} has a mate-2 file but no mate-1 file ({mates['2'].name})")
        samples[name] = Sample(name, mates["1"], mates.get("2"))          # mate 1 alone: single-end
    for name, p in singles.items():
        if name in samples:
            raise ValueError(f"{folder}: sample {name} has both paired files and a single file {p.name}")
        samples[name] = Sample(name, p)
    by_sample: dict[str, list] = {}
    for name, p, clf in outputs:
        if name not in samples:
            raise ValueError(f"{folder}: classifier output {p.name} matches no read files (no reads for sample "
                             f"{name!r}); rename it to <sample> + extension, or use a sample sheet")
        by_sample.setdefault(name, []).append((p, clf))
    for name, cands in by_sample.items():
        if len(cands) > 1:
            raise ValueError(f"{folder}: sample {name} has more than one classifier output: "
                             f"{', '.join(c[0].name for c in cands)}; keep one or use a sample sheet")
        samples[name].classifier_output, samples[name].classifier = cands[0]
    out = sorted(samples.values(), key=lambda s: s.name)
    validate_samples(out)
    return out


def detection_table(samples: list[Sample]) -> str:
    """The sample-to-file matching, one line per sample (printed by `detect` and logged by `run`)."""
    rows = [("sample", "R1", "R2", "classifier output", "classifier")]
    for s in samples:
        rows.append((s.name, Path(s.r1).name, Path(s.r2).name if s.r2 else "-",
                     Path(s.classifier_output).name if s.classifier_output else "-", s.classifier))
    width = [max(len(r[i]) for r in rows) for i in range(5)]
    return "\n".join("  ".join(c.ljust(w) for c, w in zip(r, width)).rstrip() for r in rows)


def write_samplesheet(samples: list[Sample], path: str | Path) -> None:
    """Written through a unique temporary file and renamed into place: a link at `path` is replaced, never
    followed, so the sheet can never overwrite an input."""
    from .io_utils import staged
    with staged(path) as tmp, open(tmp, "w") as fo:
        fo.write("sample\tr1\tr2\tclassifier_output\tclassifier\n")
        for s in samples:
            fo.write(f"{s.name}\t{Path(s.r1).resolve()}\t{Path(s.r2).resolve() if s.r2 else '.'}\t"
                     f"{Path(s.classifier_output).resolve() if s.classifier_output else '.'}\t{s.classifier}\n")
