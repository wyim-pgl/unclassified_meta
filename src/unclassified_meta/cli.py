"""Command-line interface: run, extract, detect, download-refs, config, panel-build, panel-quant, taxa-table,
control-reads, control-summary, da, seal-refs.

Every command that writes runs `_guard` before its first write (log file and sample sheet included): it refuses
an output folder that holds, or lies inside, an input, and takes an exclusive lock on the output folder so a
second command on the same folder stops at once."""
from __future__ import annotations

import argparse
import contextlib
import sys
from pathlib import Path

import yaml

from . import __version__
from .config import detect_samples, detection_table, load_config, read_samplesheet, write_samplesheet
from .io_utils import folder_lock, preflight
from .runner import StepError, setup_logging

# arguments that never name an input file
_NOT_INPUTS = {"cmd", "outdir", "group", "formula", "reference_level", "rand_formula", "p_adjust", "categories",
               "suffix", "rank", "namespace", "drop_lineage", "panel_lineage", "report_format", "preset",
               "reads_scope", "prefix", "classifier", "n_reads", "mode", "subsample", "only"}


def _samples(args):
    if bool(args.samples) == bool(args.input_dir):
        raise SystemExit("give exactly one of --samples SHEET.tsv or --input-dir FOLDER")
    return read_samplesheet(args.samples) if args.samples else detect_samples(args.input_dir)


def _paths_in(value):
    """Existing paths named by an argument value (lists and NAME=PATH host entries included)."""
    if isinstance(value, (list, tuple)):
        return [x for v in value for x in _paths_in(v)]
    if isinstance(value, dict):
        return [x for v in value.values() for x in _paths_in(v)]
    if not isinstance(value, (str, Path)) or not str(value):
        return []
    out = [p for p in {str(value), str(value).partition("=")[2]} if p and Path(p).exists()]
    return out


def _inputs(args, samples=(), extra=()) -> list:
    """Every input of a command: the files and folders its arguments name, and every sample's files."""
    ins = [x for k, v in vars(args).items() if k not in _NOT_INPUTS for x in _paths_in(v)]
    ins += [x for s in samples for x in (s.r1, s.r2, s.classifier_output) if x]
    return ins + [x for v in extra for x in _paths_in(v)]


def _run_references(run_dir) -> dict:
    """References a finished run recorded (config.used.yaml): panel-build inherits them, so they are inputs."""
    f = Path(run_dir) / "config.used.yaml"
    try:
        return (yaml.safe_load(f.read_text()) or {}).get("references", {}) if f.exists() else {}
    except yaml.YAMLError:
        return {}


def _guard(stack, command: str, outdir, inputs) -> Path:
    """Before the first write: refuse overlap of the output folder with any input, then lock the folder."""
    outdir = Path(outdir)
    preflight(command, inputs, owned_dirs=[outdir])
    stack.enter_context(folder_lock(outdir, command))
    return outdir


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
            "assembler": getattr(args, "assembler", None), "keep_work": getattr(args, "keep_work", None) or None,
            "min_score": getattr(args, "min_score", None), "min_hitlen": getattr(args, "min_hitlen", None)}
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


def _threshold_args(p):
    g = p.add_argument_group("Centrifuge / Centrifuger thresholds (not for Kraken2)")
    g.add_argument("--min-score", type=int, default=None,
                   help="assignments with a score below this count as unclassified (same rule as centrifuger-quant --min-score)")
    g.add_argument("--min-hitlen", type=int, default=None,
                   help="assignments with a hit length below this count as unclassified (centrifuger-quant --min-length)")


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
    _threshold_args(r)

    e = sub.add_parser("extract", help="only write each sample's unclassified reads (FASTQ.gz)")
    _input_args(e)
    _threshold_args(e)

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
    dl.add_argument("--allow-unverified-silva", action="store_true",
                    help="use SILVA files even when their published MD5 cannot be fetched (recorded as unverified)")

    c = sub.add_parser("config", help="print the default configuration as YAML")
    c.add_argument("-o", "--outdir", default=None)

    pb = sub.add_parser("panel-build", help="build a quantification panel (features + decoys) from a finished run")
    pb.add_argument("--run", required=True, help="output folder of a finished `unclassified-meta run`")
    pb.add_argument("-o", "--outdir", default="panel")
    pb.add_argument("--categories", default="virus",
                    help="comma-separated contig categories that become features (default: virus)")
    pb.add_argument("--contigs", help="file with contig names to use as features instead of --categories")
    pb.add_argument("--refs", help="download-refs folder (default: the one the run used)")
    pb.add_argument("--virus", help="virus genome FASTA to pull reference genomes from (default: from the run)")
    pb.add_argument("--rrna", help="rRNA FASTA used as decoy (default: from the run)")
    pb.add_argument("--host", action="append", metavar="NAME=GENOME.fa",
                    help="host genome(s) used as decoy (default: from the run)")
    pb.add_argument("--extra-targets", help="FASTA of extra target sequences; each sequence becomes one feature")
    pb.add_argument("--same-identity", type=float, default=95.0,
                    help="read-weighted nucleotide identity from which a bin is called X instead of X-like (default 95)")
    pb.add_argument("--no-reference-genomes", action="store_true", help="do not add the reference genomes that were hit")
    pb.add_argument("--no-host-decoy", action="store_true")
    pb.add_argument("--no-rrna-decoy", action="store_true")
    pb.add_argument("--no-other-contig-decoy", action="store_true",
                    help="do not add the run's other contigs as decoys")
    pb.add_argument("-t", "--threads", type=int, default=8)

    pq = sub.add_parser("panel-quant", help="count all reads of every sample against a panel (one fragment policy)")
    pq.add_argument("--panel", required=True, help="folder made by panel-build")
    _input_args(pq)
    pq.add_argument("-k", "--max-alignments", type=int, default=10, help="alignments reported per read (bowtie2 -k, default 10)")
    pq.add_argument("--preset", default="--sensitive", help="bowtie2 preset (default --sensitive)")
    pq.add_argument("--reads-scope", choices=["auto", "all", "unclassified"], default="auto",
                    help="what the input reads are: all reads of the samples, or the unclassified_reads/ written by "
                         "`extract`/`run` (auto: recognised from the folder). Recorded for `taxa-table --merge-panel`.")
    pq.add_argument("--input-summary", help="tables/input_summary.tsv of the extract run (found automatically "
                                            "when --input-dir is its unclassified_reads/ folder)")
    pq.add_argument("--ambiguous", choices=["ledger", "em"], default="ledger",
                    help="fragments whose best alignments fall in several features: 'ledger' books them as ambiguous "
                         "and counts them for none (default); 'em' shares them by the features' estimated abundance")
    pq.add_argument("--imputations", type=int, default=0,
                    help="with --ambiguous em: also write this many tables with every ambiguous fragment drawn to one "
                         "candidate (for `da --imputations`)")
    pq.add_argument("--seed", type=int, default=1, help="random seed of the imputations (default 1)")
    pq.add_argument("-t", "--threads", type=int, default=8)

    tt = sub.add_parser("taxa-table", help="single-rank taxon x sample counts from Kraken-style reports, or from "
                                           "Centrifuge(r) per-read output with --per-read (for `da`)")
    g = tt.add_argument_group("input (choose one)")
    g.add_argument("--input-dir", help="folder with one file per sample (use with --suffix)")
    g.add_argument("--reports", help="TSV with the columns sample and report (the report, or the per-read file)")
    tt.add_argument("--suffix", help="file-name ending of the files in --input-dir, e.g. .quant.kreport or "
                                     ".cfr.tsv; the sample name is the file name without it")
    tt.add_argument("--report-format", choices=["kraken2", "centrifuger-kreport", "centrifuger-quant"],
                    help="required for reports: kraken2 (Kraken2 --report), centrifuger-kreport, or centrifuger-quant "
                         "(--output-format 3, rows rounded one by one); it decides how exactly the tree must add up")
    pr = tt.add_argument_group("Centrifuge / Centrifuger per-read output instead of reports")
    pr.add_argument("--per-read", action="store_true",
                    help="the files are per-read classifier outputs (readID seqID taxID score ...), not reports")
    pr.add_argument("--nodes", help="nodes.dmp of the taxonomy the index was built with "
                                    "(or: centrifuger-inspect -x INDEX --taxonomy-tree > nodes.dmp)")
    pr.add_argument("--names", help="names.dmp of that taxonomy (or: centrifuger-inspect -x INDEX --name-table)")
    pr.add_argument("--min-score", type=int, default=0, help="same value as for `extract` / `run` (default 0)")
    pr.add_argument("--min-hitlen", type=int, default=0, help="same value as for `extract` / `run` (default 0)")
    pr.add_argument("--assign", choices=["lca", "em"], default="lca",
                    help="reads with several accepted taxIDs (Centrifuger -k > 1): 'lca' gives their lowest common "
                         "ancestor (default); 'em' shares them among the candidates by their estimated abundance "
                         "(expected counts)")
    pr.add_argument("--imputations", type=int, default=0,
                    help="with --assign em: also write this many count tables with every ambiguous read drawn to one "
                         "candidate (for `da --imputations`; 10-20 is usual)")
    pr.add_argument("--seed", type=int, default=1, help="random seed of the imputations (default 1)")
    tt.add_argument("--rank", default="S", choices=["D", "K", "P", "C", "O", "F", "G", "S"],
                    help="Kraken rank code of the taxa (default S = species); one rank only")
    tt.add_argument("--drop-lineage", default="", help="comma-separated names or taxIDs; these rows and everything "
                                                       "below them are left out (e.g. Viridiplantae,Metazoa)")
    tt.add_argument("--namespace", "--taxonomy-name", dest="namespace", default=None,
                    help="prefix of the feature IDs; name the index/taxonomy used, e.g. ncbi_core_nt_2025 "
                         "(default taxid; required with --per-read)")
    tt.add_argument("-o", "--outdir", default="taxa")
    m = tt.add_argument_group("merge with a panel table (background features for `da`)")
    m.add_argument("--merge-panel", metavar="QUANT/counts.tsv", help="counts.tsv of panel-quant for the same samples")
    m.add_argument("--merge-panel-features", help="features.tsv of the panel (default: next to the counts)")
    m.add_argument("--panel-taxonomy", metavar="TSV",
                   help="feature_id -> taxid of the panel features (from `panel-taxonomy`): a feature assigned to a "
                        "taxon at or below --rank is added to that taxon's row, so one organism is one feature "
                        "(needs --nodes/--names)")
    m.add_argument("--panel-lineage", help="lineage(s) the panel covers, e.g. Viruses; dropped from the reports. "
                                           "Required when panel-quant counted ALL reads; 'none' to state that the "
                                           "classifier database does not hold the panel's organisms")

    cr = sub.add_parser("control-reads", help="write reversed (not complemented) control reads for a chance-hit test")
    _input_args(cr)
    cr.add_argument("--n-reads", default="100000", help="reads (pairs) per sample, or 'all' (default 100000)")
    cr.add_argument("--seed", type=int, default=11)
    cr.add_argument("--prefix", default="CTRL_", help="prefix put in front of the read IDs (default CTRL_)")

    cs = sub.add_parser("control-summary", help="chance-hit rate of control reads over a grid of thresholds")
    cs.add_argument("per_read", nargs="+", help="per-read classifier output(s) of the control reads; files may also "
                                                "hold real reads (IDs without the prefix)")
    cs.add_argument("--classifier", choices=["centrifuger", "centrifuge", "kraken2"], default="centrifuger")
    cs.add_argument("--min-score-grid", default="0,16,25,36,49,64,81,100,150,200,300")
    cs.add_argument("--min-hitlen-grid", default="0")
    cs.add_argument("--target-pct", type=float, default=1.0, help="accepted chance-hit rate in %% (default 1)")
    cs.add_argument("--prefix", default="CTRL_")
    cs.add_argument("--all-control", action="store_true", help="every read in the files is a control read")
    cs.add_argument("-o", "--outdir", default="control")

    da = sub.add_parser("da", help="differential abundance with ANCOM-BC2 (needs the R environment)")
    da.add_argument("--counts", required=True, help="counts.tsv from panel-quant")
    da.add_argument("--features", required=True, help="features.tsv from panel-quant")
    da.add_argument("--metadata", required=True, help="TSV with a 'sample' column, the group column and covariates")
    da.add_argument("--group", required=True, help="metadata column that holds the groups to compare")
    da.add_argument("--formula", help="fixed-effects formula, e.g. 'group + batch' (default: the group column)")
    da.add_argument("--reference-level", help="group used as the baseline (default: first in alphabetical order)")
    da.add_argument("--rand-formula", metavar="'(1|unit)'",
                    help="random effects for repeated measures, e.g. '(1|unit)' when several samples come from "
                         "the same plant, plot or animal ('unit' is a metadata column)")
    da.add_argument("--lme-relax-checks", action="store_true",
                    help="with --rand-formula: switch off lme4's checks on observations per unit (sparse features)")
    da.add_argument("--prv-cut", type=float, default=0.10, help="drop features present in fewer than this share of samples")
    da.add_argument("--lib-cut", type=float, default=0, help="drop samples with fewer assigned fragments than this")
    da.add_argument("--struc-zero", action="store_true",
                    help="treat a feature absent from a whole group as a structural zero (off by default)")
    da.add_argument("--neg-lb", action="store_true",
                    help="with --struc-zero: also treat features whose abundance in a group has a negative lower "
                         "confidence bound (rare but present) as structural zeros (ANCOM-BC2 neg_lb; off by default)")
    da.add_argument("--imputations", action="store_true",
                    help="also fit every imputation table sealed next to --counts (ambiguous reads drawn by "
                         "taxa-table/panel-quant --imputations) and combine the fits with Rubin's rules "
                         "(mi_results.tsv), so that the test includes the uncertainty of the read assignment")
    da.add_argument("--alpha", type=float, default=0.05)
    da.add_argument("--p-adjust", default="BH",
                    help="multiple-testing adjustment per term (R p.adjust method; default BH, the false discovery rate; holm controls the family-wise error rate as ANCOM-BC2's own default does)")
    da.add_argument("-t", "--threads", type=int, default=1)
    da.add_argument("-o", "--outdir", default="da")
    da.add_argument("--rscript", help="Rscript of the R environment (default: $UNCLASSIFIED_META_RSCRIPT, then PATH)")
    da.add_argument("--unsealed-input", action="store_true",
                    help="use count tables that were not made by panel-quant / taxa-table (no manifest); they are "
                         "used as they are, without the completeness and checksum checks")

    ptx = sub.add_parser("panel-taxonomy", help="taxID of each panel feature from the classifier's own taxonomy "
                                               "(for taxa-table --merge-panel --panel-taxonomy)")
    ptx.add_argument("--panel", required=True, help="folder made by panel-build")
    ptx.add_argument("--write-fasta", metavar="OUT.fa", help="write the panel's target sequences (step 1), then "
                                                            "classify them with the index used for the reads")
    ptx.add_argument("--classified", help="per-read classifier output of those sequences (step 3)")
    ptx.add_argument("--nodes", help="nodes.dmp of the index taxonomy")
    ptx.add_argument("--names", help="names.dmp of the index taxonomy")
    ptx.add_argument("--min-score", type=int, default=0)
    ptx.add_argument("--min-hitlen", type=int, default=0)
    ptx.add_argument("--min-support", type=float, default=0.8,
                     help="share of a feature's classified length its taxon must hold (default 0.8)")
    ptx.add_argument("--min-classified", type=float, default=0.5,
                     help="share of a feature's length that must be classified at all (default 0.5)")
    ptx.add_argument("-o", "--outdir", default="panel_taxonomy")

    sr = sub.add_parser("seal-refs", help="adopt a reference folder prepared by hand or copied from another "
                                          "machine: check its files and write the manifests that `run` requires")
    sr.add_argument("refs_dir", help="the folder with viral/, rrna/, protein/ (the layout of download-refs)")
    sr.add_argument("--version", dest="versions", action="append", default=[], metavar="NAME=VERSION",
                    help="release of a reference, e.g. refseq_viral=230 (repeatable; unrecorded ones are 'unknown')")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        with contextlib.ExitStack() as stack:
            return _dispatch(args, stack)
    except (StepError, ValueError, FileNotFoundError, RuntimeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def _dispatch(args, stack) -> int:
    if True:
        if args.cmd == "config":
            print(yaml.safe_dump(load_config(None), sort_keys=False))
        elif args.cmd == "detect":
            samples = detect_samples(args.input_dir)
            # detect writes one file and owns no folder (its default folder is the current one): check that file
            out = Path(args.outdir) / "samplesheet.detected.tsv"
            preflight("detect", _inputs(args, samples), files=[out])
            write_samplesheet(samples, out)
            print(detection_table(samples))
            print(f"wrote {out}")
        elif args.cmd == "extract":
            from .pipeline import extract_all
            samples = _samples(args)
            outdir = _guard(stack, "extract", args.outdir or "results", _inputs(args, samples))
            setup_logging(outdir / "logs" / "pipeline.log")
            write_samplesheet(samples, outdir / "samplesheet.used.tsv")
            df = extract_all(samples, outdir, args.min_score or 0, args.min_hitlen or 0)
            print(df.to_string(index=False))
            print(f"unclassified reads written to {outdir / 'unclassified_reads'}")
        elif args.cmd == "download-refs":
            from . import refs
            outdir = _guard(stack, "download-refs", args.outdir, _inputs(args))
            setup_logging(outdir / "download-refs.log")
            if args.dry_run:
                rows = refs.plan(outdir)
                for name, url, dest in rows:
                    n = refs.remote_size(url)
                    print(f"{name}\t{'?' if n is None else f'{n / 1e6:.1f} MB'}\t{url}\t{dest}")
                refs.size_check(outdir, rows)
            else:
                vt = refs.download_refs(outdir, args.threads, Path(args.rrna_fasta) if args.rrna_fasta else None,
                                        args.keep_sources, tuple(args.only), args.allow_unverified_silva)
                print(open(vt).read())
        elif args.cmd == "panel-build":
            from . import panel
            _guard(stack, "panel-build", args.outdir, _inputs(args, extra=[_run_references(args.run)]))
            setup_logging(Path(args.outdir) / "panel-build.log")
            out = panel.build_panel(
                Path(args.run), Path(args.outdir), [c.strip() for c in args.categories.split(",") if c.strip()],
                [] if args.no_host_decoy else _hosts(args.host), None if args.no_rrna_decoy else args.rrna,
                args.virus, args.refs, args.extra_targets, args.contigs, args.same_identity,
                not args.no_reference_genomes, not args.no_other_contig_decoy, args.threads,
                use_run_hosts=not args.no_host_decoy, use_run_rrna=not args.no_rrna_decoy)
            print(open(out / "features.tsv").read())
            print(f"panel: {out}")
        elif args.cmd in ("panel-quant", "taxa-table") and getattr(args, "seed", 0) < 0:
            raise SystemExit("--seed must be 0 or a positive whole number; nothing was written")
        elif args.cmd == "panel-quant":
            from . import panel
            samples = _samples(args)
            outdir = _guard(stack, "panel-quant", args.outdir or "quant", _inputs(args, samples))
            setup_logging(outdir / "panel-quant.log")
            write_samplesheet(samples, outdir / "samplesheet.used.tsv")
            out = panel.quantify(Path(args.panel), samples, outdir, args.threads, args.max_alignments, args.preset,
                                 args.reads_scope, args.input_summary, args.ambiguous, args.imputations, args.seed)
            print(open(outdir / "ledger.tsv").read())
            print(f"counts: {out}")
        elif args.cmd == "taxa-table":
            from . import taxa
            if bool(args.input_dir) == bool(args.reports):
                raise SystemExit("give exactly one of --input-dir FOLDER (with --suffix) or --reports SHEET.tsv")
            # every indirect input is resolved before the guard: the reports a sheet lists, the panel's tables and
            # the extract summary its provenance names
            reports = taxa.find_reports(args.input_dir, args.suffix) if args.input_dir else taxa.read_report_sheet(args.reports)
            extra = list(reports.values())
            if args.merge_panel:
                mp = Path(args.merge_panel)
                extra += [mp.parent / "features.tsv", taxa.panel_scope(mp)[1]]
                if args.panel_taxonomy:
                    extra.append(Path(args.panel_taxonomy))
                # the panel's imputation draws are inputs too (a linked folder must not be written over)
                extra += sorted((mp.parent / "imputations").glob("*.tsv")) if (mp.parent / "imputations").is_dir() else []
            elif args.panel_taxonomy:
                raise SystemExit("--panel-taxonomy is used together with --merge-panel")
            if args.panel_taxonomy and not args.per_read:
                raise SystemExit("--panel-taxonomy needs --per-read (the taxa table must record the taxonomy it was "
                                 "made with); nothing was written")
            if args.panel_taxonomy and taxa.panel_scope(args.merge_panel)[0] != "unclassified":
                raise SystemExit("--panel-taxonomy needs a panel counted on the unclassified reads (panel-quant on "
                                 "the unclassified_reads/ of `run`/`extract`); nothing was written")
            outdir = _guard(stack, "taxa-table", args.outdir, _inputs(args, extra=extra))
            setup_logging(outdir / "taxa-table.log")
            drop = [x.strip() for x in args.drop_lineage.split(",") if x.strip()]
            strict = []          # panel lineages of a panel counted on all reads: dropped, and checked against --rank
            if args.merge_panel:
                lineages = taxa.check_merge_mode(args.merge_panel, args.panel_lineage)
                if taxa.panel_scope(args.merge_panel)[0] == "all":
                    strict = list(lineages)             # checked even when also given in --drop-lineage
                    drop = [x for x in drop if x.lower() not in {y.lower() for y in strict}]
                else:
                    drop += [x for x in lineages if x not in drop]
            elif args.panel_lineage:
                raise SystemExit("--panel-lineage is used together with --merge-panel")
            if args.per_read:
                if not (args.nodes and args.names):
                    raise SystemExit("--per-read needs --nodes nodes.dmp and --names names.dmp (the taxonomy of the index)")
                if not args.namespace:
                    raise SystemExit("--per-read needs --taxonomy-name (a short name of the index taxonomy, "
                                     "e.g. ncbi_core_nt_2025); it becomes part of the feature IDs")
                taxa.build_table_from_per_read(reports, outdir, args.nodes, args.names, args.rank, tuple(drop),
                                               args.namespace, args.min_score, args.min_hitlen, tuple(strict),
                                               args.assign, args.imputations, args.seed)
            else:
                if args.assign != "lca" or args.imputations:
                    raise SystemExit("--assign em and --imputations need --per-read: a report has already resolved "
                                     "each read to one taxon")
                if args.min_score or args.min_hitlen or ((args.nodes or args.names) and not args.panel_taxonomy):
                    raise SystemExit("--min-score, --min-hitlen, --nodes and --names belong to --per-read; a report "
                                     "already has its thresholds applied (centrifuger-quant --min-score/--min-length)")
                if not args.report_format:
                    raise SystemExit("give --report-format kraken2, centrifuger-kreport or centrifuger-quant "
                                     "(the format is not guessed)")
                taxa.build_table(reports, outdir, args.rank, tuple(drop), args.namespace or "taxid", tuple(strict),
                                 report_format=args.report_format)
            print(open(outdir / "ledger.tsv").read())
            print(f"taxa table: {outdir / 'counts.tsv'}")
            if args.merge_panel:
                out = taxa.merge_with_panel(outdir, args.merge_panel, args.merge_panel_features, outdir / "merged",
                                            args.panel_taxonomy, args.nodes, args.names, tuple(drop) + tuple(strict))
                print(open(outdir / "merged" / "reconciliation.tsv").read())
                print(f"merged table for `da`: {out}")
        elif args.cmd == "control-reads":
            from . import control
            from .config import parse_subsample
            samples = _samples(args)
            outdir = _guard(stack, "control-reads", args.outdir or "control", _inputs(args, samples))
            setup_logging(outdir / "control-reads.log")
            df = control.make_controls(samples, outdir, parse_subsample(args.n_reads), args.seed, args.prefix)
            print(df.to_string(index=False))
            print(f"control reads: {outdir / 'control_reads'}; classify them with the same index and options as the "
                  "real reads, then run `unclassified-meta control-summary` on the per-read output")
        elif args.cmd == "control-summary":
            from . import control
            _guard(stack, "control-summary", args.outdir, _inputs(args))
            setup_logging(Path(args.outdir) / "control-summary.log")
            grid = lambda t: tuple(int(x) for x in t.split(",") if x.strip())  # noqa: E731
            df = control.summarise(args.per_read, args.outdir, args.classifier, grid(args.min_score_grid),
                                   grid(args.min_hitlen_grid), args.target_pct, args.prefix, args.all_control)
            print(df.to_string(index=False))
            best = df[df["recommended"]]
            if len(best):
                b = best.iloc[0]
                print(f"least strict thresholds with a chance-hit rate <= {args.target_pct} %: "
                      f"--min-score {int(b.min_score)} --min-hitlen {int(b.min_hitlen)} "
                      f"(control reads classified: {b.control_classified_pct} %)")
            else:
                print(f"no threshold in the grid reaches a chance-hit rate <= {args.target_pct} %; extend the grids")
        elif args.cmd == "da":
            from . import da
            ins = _inputs(args)
            if getattr(args, "imputations", False):     # the imputation tables are inputs too
                imp = Path(args.counts).parent / "imputations"
                ins += sorted(imp.glob("*.tsv")) if imp.is_dir() else []
            ins = [x for x in ins if str(x) != str(args.rscript)] if args.rscript else ins
            _guard(stack, "da", args.outdir, ins)
            da.run_ancombc2(args, ins)
        elif args.cmd == "panel-taxonomy":
            from . import panel_taxonomy as ptx
            if bool(args.write_fasta) == bool(args.classified):
                raise SystemExit("give --write-fasta OUT.fa (step 1) or --classified FILE with --nodes/--names (step 3)")
            if args.write_fasta:
                print(f"target sequences: {ptx.write_fasta(args.panel, args.write_fasta)}")
            else:
                if not (args.nodes and args.names):
                    raise SystemExit("--classified needs --nodes and --names (the taxonomy of the index)")
                import os
                panel_files = [Path(r) / f for r, _, fs in os.walk(args.panel, followlinks=True) for f in fs]
                ins = [*panel_files, Path(args.classified), Path(args.nodes), Path(args.names)]
                outdir = _guard(stack, "panel-taxonomy", args.outdir, ins)
                setup_logging(outdir / "panel-taxonomy.log")
                out = ptx.build(args.panel, args.classified, args.nodes, args.names, outdir, args.min_score,
                                args.min_hitlen, args.min_support, args.min_classified)
                print(open(out).read())
                print(f"panel taxonomy: {out}")
        elif args.cmd == "seal-refs":
            from . import refs
            versions = dict(v.split("=", 1) for v in args.versions if "=" in v)
            _guard(stack, "seal-refs", args.refs_dir, [])
            vt = refs.seal_refs(Path(args.refs_dir), versions)
            print(open(vt).read())
        elif args.cmd == "run":
            from .pipeline import run
            cfg = _cfg(args)
            samples = _samples(args)
            _guard(stack, "run", cfg["outdir"], _inputs(args, samples, [cfg.get("references", {})]))
            rep = run(cfg, samples)
            print(f"report: {rep}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
