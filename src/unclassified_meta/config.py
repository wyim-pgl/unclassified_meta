"""Configuration defaults, YAML loading, sample sheets and input-folder detection."""
from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from pathlib import Path

import yaml

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


def validate_samples(samples: list[Sample]) -> None:
    if not samples:
        raise ValueError("no samples given")
    for s in samples:
        for f in (s.r1, s.r2, s.classifier_output):
            if f is not None and not Path(f).exists():
                raise FileNotFoundError(f"sample {s.name}: file not found: {f}")
        if not re.fullmatch(r"[A-Za-z0-9._-]+", s.name):
            raise ValueError(f"sample name {s.name!r}: use letters, digits, '.', '_' or '-' only")


FASTQ_RE = re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])(?:R)?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$")
SINGLE_RE = re.compile(r"^(?P<sample>.+?)\.(?:fastq|fq)(?:\.gz)?$")
CLASSIFIER_EXT = {
    ".kraken": "kraken2", ".kraken2": "kraken2", ".kraken.out": "kraken2", ".kraken2.out": "kraken2",
    ".k2": "kraken2", ".centrifuge": "centrifuge", ".centrifuge.out": "centrifuge",
    ".cfr": "centrifuger", ".centrifuger": "centrifuger", ".centrifuger.out": "centrifuger",
}


def _looks_like_report(path: Path) -> bool:
    """Kraken/centrifuge *reports* (percent-first summaries) are not per-read outputs."""
    name = path.name.lower()
    if "report" in name or name.endswith((".kreport", ".kreport2")):
        return True
    with open(path) as fh:
        first = fh.readline().split("\t")
    try:
        float(first[0])
        return True          # report lines start with a percentage
    except (ValueError, IndexError):
        return False


def detect_samples(folder: str | Path) -> list[Sample]:
    """Find FASTQ files (paired *_R1/_R2, *_1/_2, or single) and matching per-read classifier outputs."""
    folder = Path(folder)
    files = sorted(p for p in folder.iterdir() if p.is_file())
    pairs: dict[str, dict] = {}
    singles: dict[str, Path] = {}
    for p in files:
        m = FASTQ_RE.match(p.name)
        if m:
            pairs.setdefault(m["sample"], {})[m["mate"]] = p
            continue
        m = SINGLE_RE.match(p.name)
        if m:
            singles[m["sample"]] = p
    samples = []
    for name, mates in pairs.items():
        if "1" not in mates:
            raise ValueError(f"{folder}: sample {name} has a mate-2 file but no mate-1 file")
        samples.append(Sample(name, mates["1"], mates.get("2")))
    for name, p in singles.items():
        if name not in pairs:
            samples.append(Sample(name, p))
    for s in samples:
        for p in files:
            low = p.name.lower()
            if not low.startswith(s.name.lower()):
                continue
            rest = low[len(s.name):]
            for ext, clf in CLASSIFIER_EXT.items():
                if rest == ext and not _looks_like_report(p):
                    s.classifier_output, s.classifier = p, clf
    samples.sort(key=lambda s: s.name)
    validate_samples(samples)
    return samples


def write_samplesheet(samples: list[Sample], path: str | Path) -> None:
    with open(path, "w") as fo:
        fo.write("sample\tr1\tr2\tclassifier_output\tclassifier\n")
        for s in samples:
            fo.write(f"{s.name}\t{Path(s.r1).resolve()}\t{Path(s.r2).resolve() if s.r2 else '.'}\t"
                     f"{Path(s.classifier_output).resolve() if s.classifier_output else '.'}\t{s.classifier}\n")
