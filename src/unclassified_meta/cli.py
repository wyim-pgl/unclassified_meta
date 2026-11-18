"""Command-line interface: run, extract, download-refs, detect, report-config."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from . import __version__
from .config import detect_samples, load_config, read_samplesheet, write_samplesheet
from .runner import StepError, setup_logging


def _samples(args):
    if bool(args.samples) == bool(args.input_dir):
        raise SystemExit("give exactly one of --samples SHEET.tsv or --input-dir FOLDER")
    return read_samplesheet(args.samples) if args.samples else detect_samples(args.input_dir)


def _hosts(values):
    out = []
    for v in values or []:
        name, _, path = v.partition("=")
        if not path:
            name, path = Path(v).name.split(".")[0], v
        out.append({"name": name, "fasta": path})
    return out


def _cfg(args) -> dict:
    refs = {}
    for key, attr in (("refs_dir", "refs"), ("rrna", "rrna"), ("virus", "virus"), ("virus_protein", "virus_protein"),
                      ("protein", "protein"), ("organelle", "organelle"), ("extra_virus_fasta", "extra_virus_fasta"),
                      ("host_proteome", "host_proteome"), ("kraken2_db", "kraken2_db")):
        v = getattr(args, attr, None)
        if v:
            refs[key] = v
    if getattr(args, "host", None):
        refs["host"] = _hosts(args.host)
    over = {"mode": getattr(args, "mode", None), "outdir": args.outdir, "threads": getattr(args, "threads", None),
            "memory_gb": getattr(args, "memory", None), "subsample": getattr(args, "subsample", None),
            "seed": getattr(args, "seed", None), "min_contig_len": getattr(args, "min_contig_len", None),
            "assembler": getattr(args, "assembler", None), "keep_work": getattr(args, "keep_work", None) or None}
    if refs:
        over["references"] = refs
    if getattr(args, "tblastx", False):
        over["thresholds"] = {"virus": {"tblastx_extra": True}}
    return load_config(getattr(args, "config", None), over)


def _input_args(p):
    g = p.add_argument_group("input (choose one)")
    g.add_argument("--samples", help="sample sheet TSV: sample, r1, [r2], [classifier_output], [classifier]")
    g.add_argument("--input-dir", help="folder with FASTQ files (+ per-read classifier outputs); samples are detected")
    p.add_argument("-o", "--outdir", default=None, help="output folder (default: results)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="unclassified-meta", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="full pipeline")
    _input_args(r)
    r.add_argument("-c", "--config", help="YAML config (command-line options override it)")
    r.add_argument("--mode", choices=["rna", "dna"], help="RNA-seq (rnaSPAdes) or DNA shotgun (MEGAHIT); default rna")
    r.add_argument("--refs", help="folder made by `download-refs` (fills --rrna/--virus/--virus-protein/--protein)")
    r.add_argument("--host", action="append", metavar="NAME=GENOME.fa",
                   help="host genome; repeat for several hosts, checked in the given order (required)")
    r.add_argument("--rrna", help="rRNA FASTA (default from --refs)")
    r.add_argument("--virus", help="virus genome FASTA (default from --refs: RefSeq viral)")
    r.add_argument("--virus-protein", help="virus protein .dmnd/FASTA (default from --refs)")
    r.add_argument("--protein", help="plant protein .dmnd/FASTA (default from --refs: Swiss-Prot Viridiplantae)")
    r.add_argument("--host-proteome", help="optional host protein FASTA added to the protein layer")
    r.add_argument("--organelle", help="optional plastid/mitochondrial genome FASTA")
    r.add_argument("--extra-virus-fasta", help="optional small custom virus genome FASTA (BLASTN)")
    r.add_argument("--tblastx", action="store_true", help="also run TBLASTX against --extra-virus-fasta")
    r.add_argument("--kraken2-db", help="optional Kraken2 database to classify contigs")
    r.add_argument("--subsample", help="reads (pairs) per sample: a number (default 1000000) or 'all'")
    r.add_argument("--seed", type=int, help="subsampling seed (default 11)")
    r.add_argument("--memory", type=float, help="memory for the assembler in GB (default 32)")
    r.add_argument("-t", "--threads", type=int, help="threads (default 8)")
    r.add_argument("--assembler", choices=["auto", "rnaspades", "megahit", "metaspades"])
    r.add_argument("--min-contig-len", type=int, help="minimum contig length (default 300)")
    r.add_argument("--keep-work", action="store_true", help="keep pooled reads in work/")

    e = sub.add_parser("extract", help="only write each sample's unclassified reads (FASTQ.gz)")
    _input_args(e)

    d = sub.add_parser("detect", help="show the samples found in --input-dir and write the sample sheet")
    d.add_argument("--input-dir", required=True)
    d.add_argument("-o", "--outdir", default=".")

    dl = sub.add_parser("download-refs", help="download and index the required references")
    dl.add_argument("-o", "--outdir", default="refs")
    dl.add_argument("-t", "--threads", type=int, default=4)
    dl.add_argument("--rrna-fasta", help="use an existing plant rRNA FASTA instead of downloading SILVA/Rfam")
    dl.add_argument("--only", nargs="+", choices=["viral", "rrna", "protein"], default=["viral", "rrna", "protein"])
    dl.add_argument("--keep-sources", action="store_true", help="keep the downloaded .gz files")
    dl.add_argument("--dry-run", action="store_true", help="list the files and sizes, download nothing")

    c = sub.add_parser("config", help="print the default configuration as YAML")
    c.add_argument("-o", "--outdir", default=None)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.cmd == "config":
            print(yaml.safe_dump(load_config(None), sort_keys=False))
        elif args.cmd == "detect":
            samples = detect_samples(args.input_dir)
            out = Path(args.outdir) / "samplesheet.detected.tsv"
            out.parent.mkdir(parents=True, exist_ok=True)
            write_samplesheet(samples, out)
            for s in samples:
                print(f"{s.name}\t{'paired' if s.paired else 'single'}\t{s.classifier}\t"
                      f"{s.classifier_output.name if s.classifier_output else '-'}")
            print(f"wrote {out}")
        elif args.cmd == "extract":
            from .pipeline import extract_all
            samples = _samples(args)
            outdir = Path(args.outdir or "results")
            setup_logging(outdir / "logs" / "pipeline.log")
            write_samplesheet(samples, outdir / "samplesheet.used.tsv")
            df = extract_all(samples, outdir)
            print(df.to_string(index=False))
            print(f"unclassified reads written to {outdir / 'unclassified_reads'}")
        elif args.cmd == "download-refs":
            from . import refs
            outdir = Path(args.outdir)
            setup_logging(outdir / "download-refs.log")
            if args.dry_run:
                rows = refs.plan(outdir)
                for name, url, dest in rows:
                    n = refs.remote_size(url)
                    print(f"{name}\t{'?' if n is None else f'{n / 1e6:.1f} MB'}\t{url}\t{dest}")
                refs.size_check(outdir, rows)
            else:
                vt = refs.download_refs(outdir, args.threads, Path(args.rrna_fasta) if args.rrna_fasta else None,
                                        args.keep_sources, tuple(args.only))
                print(open(vt).read())
        elif args.cmd == "run":
            from .pipeline import run
            cfg = _cfg(args)
            rep = run(cfg, _samples(args))
            print(f"report: {rep}")
    except (StepError, ValueError, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
