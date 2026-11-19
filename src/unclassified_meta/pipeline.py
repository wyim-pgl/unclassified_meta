"""The pipeline: extract -> subsample -> pool -> assemble -> classify contigs -> read-back -> tables/figures/report."""
from __future__ import annotations

import json
import logging
import os
import shlex
import shutil
from pathlib import Path

import pandas as pd
import yaml

from . import __version__, manifest
from . import assemble as asm
from . import classify, plots, readback, report, summarize
from .config import Sample, check_references, detection_table, write_samplesheet
from .extract import copy_reads, filter_fastq, filter_fastq_pair, parse_classifier
from .io_utils import count_fastq, refuse_overlap, staged
from .runner import Steps, file_id, params_key, require, setup_logging, sh

log = logging.getLogger("unclassified_meta")



def classification_key(contigs: Path, cfg: dict) -> str:
    """Key of the classification step: the CONTENT of the contigs (not the assembly parameters), the content of
    every reference, the thresholds and the priority order. Changed contigs are therefore classified again,
    and so is everything that depends on them."""
    return params_key("classify", manifest.sha256(contigs), reference_identity(cfg["references"]),
                      cfg["thresholds"], cfg["priority"])


def verify_reference_folders(refs: dict) -> None:
    """Every configured reference that lies in a download-refs folder (refs_dir/<viral|rrna|protein>) must belong
    to a complete, sealed download: an interrupted download or build is never used."""
    rd = refs.get("refs_dir")
    if not rd:
        return
    rd = Path(rd).resolve()
    folders = set()
    for k, v in refs.items():
        if k in ("refs_dir", "host") or not isinstance(v, (str, Path)) or not v:
            continue
        p = Path(v).resolve()
        if p.parent.parent == rd:
            folders.add(p.parent)
    for f in sorted(folders):
        try:
            manifest.verify(f, what="reference files")
        except manifest.ManifestError as e:
            raise manifest.ManifestError(f"{e}. A reference folder prepared by hand or copied from another machine "
                                         f"is adopted with `unclassified-meta seal-refs {rd}`") from None


def _content_id(path):
    """Whole-file sha256 of a reference file (folders such as a Kraken2 database: runner.file_id)."""
    if path is None or not Path(path).exists():
        return [str(path), None]
    if Path(path).is_dir():
        return file_id(path)
    return [Path(path).stat().st_size, classify.source_digest(path)]


def reference_identity(refs: dict) -> dict:
    """References keyed by what they contain, not where they are: a file replaced in place, or a different file
    at the same path, changes the key and the contigs are classified again (see runner.file_id)."""
    out = {}
    for k, v in sorted(refs.items()):
        if k == "host":
            out[k] = [[h.get("name"), _content_id(h.get("fasta"))] for h in (v or [])]
        elif isinstance(v, (str, Path)) and v and k != "refs_dir":
            out[k] = _content_id(v)
        else:
            out[k] = v if k != "refs_dir" else None
    return out

def q(x) -> str:
    return shlex.quote(str(x))


def out_reads(outdir: Path, s: Sample) -> tuple[Path, Path | None]:
    d = outdir / "unclassified_reads"
    if s.paired:
        return d / f"{s.name}_R1.fastq.gz", d / f"{s.name}_R2.fastq.gz"
    return d / f"{s.name}.fastq.gz", None


def extract_sample(s: Sample, outdir: Path, min_score: int = 0, min_hitlen: int = 0) -> dict:
    """Write the sample's unclassified reads and return its input summary.

    `unclassified_reads` is what is written and explained by the pipeline: the reads without any assignment
    (`no_hit_reads`) plus those assigned only below the thresholds (`low_score_reads`).
    """
    o1, o2 = out_reads(outdir, s)
    o1.parent.mkdir(parents=True, exist_ok=True)
    if s.classifier == "none":
        # a copy, never a link: a later run with a classifier rewrites these files, and writing through a link
        # would overwrite the original reads
        if s.r2:
            filter_fastq_pair(s.r1, s.r2, o1, o2, None)          # all reads; the mates are checked record by record
        else:
            copy_reads(s.r1, o1)
        n = count_fastq(o1)
        return {"sample": s.name, "paired": s.paired, "classifier": "none", "total_reads": n,
                "classified_reads": 0, "low_score_reads": 0, "no_hit_reads": n, "unclassified_reads": n}
    ids, total, classified, low, all_ids = parse_classifier(s.classifier_output, s.classifier, min_score, min_hitlen)
    if classified + low + (len(ids) - low) != total:
        raise ValueError(f"{s.name}: classifier read accounting does not close "
                         f"({classified} classified + {low} low-score + {len(ids) - low} without a hit != {total})")
    # every selected read must be found exactly once, the FASTQ must hold the reads the classifier saw, and
    # paired files must list the same reads in the same order
    if s.r2:
        n1 = filter_fastq_pair(s.r1, s.r2, o1, o2, ids, all_ids)
    else:
        n1 = filter_fastq(s.r1, o1, ids, all_ids)
    return {"sample": s.name, "paired": s.paired, "classifier": s.classifier, "total_reads": total,
            "classified_reads": classified, "low_score_reads": low, "no_hit_reads": len(ids) - low,
            "unclassified_reads": n1}


def extract_all(samples: list[Sample], outdir: Path, min_score: int = 0, min_hitlen: int = 0) -> pd.DataFrame:
    min_score, min_hitlen = int(min_score or 0), int(min_hitlen or 0)
    if min_score < 0 or min_hitlen < 0:
        raise ValueError("--min-score and --min-hitlen must be 0 or larger")
    if min_score or min_hitlen:
        bad = [s.name for s in samples if s.classifier == "kraken2"]
        if bad:
            raise ValueError("--min-score / --min-hitlen apply to Centrifuge and Centrifuger per-read output only, "
                             f"but these samples have Kraken2 output: {', '.join(bad)}")
        if not any(s.classifier in ("centrifuge", "centrifuger") for s in samples):
            raise ValueError("--min-score / --min-hitlen were given, but no sample has Centrifuge or Centrifuger output")
    # before anything is written: no output of any sample may be an input of any sample
    inputs = [f for s in samples for f in (s.r1, s.r2, s.classifier_output)]
    outputs = [o for s in samples for o in out_reads(outdir, s)]
    outputs += [outdir / "work" / "inputs" / f"{s.name}.json" for s in samples] + [outdir / "tables" / "input_summary.tsv"]
    refuse_overlap(outputs, inputs, "extract")
    manifest.verify_consumed(inputs, "reads")         # reads written by an earlier extract: only a sealed set
    # the per-sample read files and input_summary.tsv are one generation: marker first, manifest last, so a
    # rerun that stops midway leaves no folder a consumer accepts
    reads_dir = outdir / "unclassified_reads"
    manifest.begin(reads_dir, "extract")
    work = outdir / "work"
    steps = Steps(work)
    rows = []
    for s in samples:
        js = work / "inputs" / f"{s.name}.json"
        js.parent.mkdir(parents=True, exist_ok=True)
        key = params_key("extract", file_id(s.r1), file_id(s.r2), file_id(s.classifier_output), s.classifier,
                         min_score, min_hitlen)
        if not steps.done(f"extract_{s.name}", key):
            log.info("extract unclassified reads: %s (%s)", s.name, s.classifier)
            steps.invalidate(f"extract_{s.name}")
            js.unlink(missing_ok=True)
            summary = extract_sample(s, outdir, min_score, min_hitlen)
            js.write_text(json.dumps(summary))
            steps.mark(f"extract_{s.name}", key, [*[o for o in out_reads(outdir, s) if o is not None], js])
        rows.append(json.loads(js.read_text()))
    df = pd.DataFrame(rows)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    with staged(outdir / "tables" / "input_summary.tsv") as tmp:
        df.to_csv(tmp, sep="\t", index=False)
    seal_extraction(outdir, samples)
    return df


def seal_extraction(outdir: Path, samples: list[Sample]) -> None:
    """Seal unclassified_reads/ together with tables/input_summary.tsv (one generation)."""
    reads = [o.name for s in samples for o in out_reads(outdir, s) if o is not None]
    manifest.seal(outdir / "unclassified_reads", reads + ["../tables/input_summary.tsv"],
                  {"step": "extract", "samples": [s.name for s in samples]})


def subsample(s: Sample, outdir: Path, n: int | None, seed: int) -> tuple[Path, Path | None, int]:
    """Return the reads used for assembly/read-back (subsampled or all) and their count."""
    r1, r2 = out_reads(outdir, s)
    if not n:
        return r1, r2, count_fastq(r1)
    d = outdir / "work" / "subsample"
    d.mkdir(parents=True, exist_ok=True)
    o1 = d / f"{s.name}_1.fq.gz"
    o2 = d / f"{s.name}_2.fq.gz" if r2 else None
    done = d / f"{s.name}.done"
    key = params_key("subsample", n, seed, file_id(r1), file_id(r2))
    state = done.read_text().split("\t") if done.exists() else []
    if len(state) != 2 or state[1] != key:      # missing, written by an older version, or other parameters
        require("seqtk")
        done.unlink(missing_ok=True)
        for src, dst in ((r1, o1), (r2, o2)):
            if src is None:
                continue
            tmp = dst.with_name(dst.name + ".part")
            sh(f"seqtk sample -s {seed} {q(src)} {n} | gzip -1 > {q(tmp)}")
            os.replace(tmp, dst)
        state = [str(count_fastq(o1)), key]
        tmpd = done.with_name(done.name + ".part")
        tmpd.write_text("\t".join(state))
        os.replace(tmpd, done)
    return o1, o2, int(state[0])


def run(cfg: dict, samples: list[Sample]) -> Path:
    outdir = Path(cfg["outdir"]).resolve()
    work = outdir / "work"
    setup_logging(outdir / "logs" / "pipeline.log")
    check_references(cfg)
    require("seqtk", "minimap2", "blastn", "makeblastdb", "diamond", "bowtie2", "bowtie2-build", "samtools")
    verify_reference_folders(cfg["references"])
    # nothing this run writes may be one of its inputs (reads, classifier outputs, references)
    ref_inputs = [v for k, v in cfg["references"].items() if isinstance(v, (str, Path)) and v and k != "refs_dir"]
    ref_inputs += [h["fasta"] for h in cfg["references"].get("host") or []]
    refuse_overlap([outdir / f for f in manifest.RUN_ARTIFACTS] + [outdir / "samplesheet.used.tsv", outdir / "report.html"],
                   [f for s in samples for f in (s.r1, s.r2, s.classifier_output)] + ref_inputs, "run")
    # the run's outputs are about to change: mark the folder as a generation and unpublish the previous one first
    manifest.begin(outdir, "run")
    steps = Steps(work)
    write_samplesheet(samples, outdir / "samplesheet.used.tsv")
    log.info("samples and files:\n%s", detection_table(samples))
    rd = cfg["references"].get("refs_dir")
    vt = Path(rd) / "REFS_VERSION.tsv" if rd else None
    (outdir / "references.used.json").write_text(json.dumps({
        "refs_dir": str(Path(rd).resolve()) if rd else None,
        "refs_version": vt.read_text() if vt and vt.exists() else "",
        "references": reference_identity(cfg["references"])}, indent=1, default=str) + "\n")
    (outdir / "config.used.yaml").write_text(yaml.safe_dump(
        {k: v for k, v in cfg.items() if not k.startswith("_")}, sort_keys=False))
    threads = int(cfg["threads"])

    inputs = extract_all(samples, outdir, cfg.get("min_score", 0), cfg.get("min_hitlen", 0))

    used = {}
    for s in samples:
        used[s.name] = subsample(s, outdir, cfg["subsample"], int(cfg["seed"]))
    inputs["reads_used"] = inputs["sample"].map({k: v[2] for k, v in used.items()})
    # the summary belongs to the sealed extraction: unpublish, rewrite, seal again
    manifest.begin(outdir / "unclassified_reads", "run")
    with staged(outdir / "tables" / "input_summary.tsv") as tmp:
        inputs.to_csv(tmp, sep="\t", index=False)
    seal_extraction(outdir, samples)

    # pool
    pool_dir = work / "pool"
    pool_dir.mkdir(parents=True, exist_ok=True)
    pool = {"r1": None, "r2": None, "s": None}
    paired = [v for v in used.values() if v[1] is not None and v[2] > 0]
    single = [v for v in used.values() if v[1] is None and v[2] > 0]
    if paired:
        pool["r1"], pool["r2"] = pool_dir / "pool_1.fq.gz", pool_dir / "pool_2.fq.gz"
    if single:
        pool["s"] = pool_dir / "pool_s.fq.gz"
    pool_key = params_key("pool", cfg["subsample"], cfg["seed"],
                          [[k, file_id(v[0]), file_id(v[1]), v[2]] for k, v in sorted(used.items())])
    assembler = asm.choose_assembler(cfg)
    asm_key = params_key("assemble", pool_key, assembler, cfg["min_contig_len"])
    pool_files = [f for f in pool.values() if f is not None]
    if not steps.done("assemble", asm_key) and not steps.done("pool", pool_key):   # done() checks every pool file
        steps.invalidate("pool")
        for files, dst in (([v[0] for v in paired], pool["r1"]), ([v[1] for v in paired], pool["r2"]),
                           ([v[0] for v in single], pool["s"])):
            if files:
                with staged(dst) as tmp:
                    sh("cat " + " ".join(q(f) for f in files) + f" > {q(tmp)}")
        steps.mark("pool", pool_key, pool_files)
    n_pool = sum(v[2] for v in used.values())
    if n_pool == 0:
        raise SystemExit("no unclassified reads to assemble")

    # assemble
    cfg["_assembler"] = assembler
    est = asm.estimate_memory_gb(assembler, n_pool)
    log.info("assembler %s, %d pooled reads (pairs counted once), estimated peak memory %.1f GB",
             assembler, n_pool, est)
    warn = asm.memory_warning(assembler, n_pool, float(cfg["memory_gb"]))
    if warn:
        log.warning(warn)
        print(warn, flush=True)
    contigs = outdir / "contigs" / "contigs.fa"
    res_file = work / "resources.json"
    resources = json.loads(res_file.read_text()) if res_file.exists() else []
    if not steps.done("assemble", asm_key):
        steps.invalidate("assemble")
        fasta, res = asm.assemble(assembler, pool, work / "assembly", threads, float(cfg["memory_gb"]))
        contigs.parent.mkdir(parents=True, exist_ok=True)
        n = asm.filter_and_rename(fasta, contigs, outdir / "tables" / "contig_names.tsv", int(cfg["min_contig_len"]))
        resources.append({"step": f"assembly ({assembler})", "input_reads": n_pool, "estimated_gb": est,
                          **res, "contigs_kept": n})
        res_file.write_text(json.dumps(resources))
        steps.mark("assemble", asm_key, [contigs, outdir / "tables" / "contig_names.tsv"])
    log.info("contigs >= %d bp: %s", cfg["min_contig_len"], contigs)

    # classify
    ctab = outdir / "tables" / "contig_category.tsv"
    cls_key = classification_key(contigs, cfg)
    contig_digest = manifest.sha256(contigs)          # read-back depends on the contigs' content, not on asm_key
    if not steps.done("classify", cls_key):
        steps.invalidate("classify")
        result = classify.Classifier(cfg, work / "classify", threads).run(contigs) if contigs.stat().st_size else {}
        tmp = ctab.with_name(ctab.name + ".part")
        classify.write_contig_table(result, tmp)
        os.replace(tmp, ctab)
        steps.mark("classify", cls_key, [ctab])
    contig_df = pd.read_csv(ctab, sep="\t", dtype={"best_hit": str, "best_hit_title": str})
    cat_of = dict(zip(contig_df["contig"], contig_df["category"]))

    # read-back
    rb = work / "readback"
    rb.mkdir(parents=True, exist_ok=True)
    counts, per_sample = {}, {}
    if len(contig_df) and not steps.done("readback_index", contig_digest):
        steps.invalidate("readback_index")
        for old in rb.glob("contigs.*bt2*"):      # index of an earlier assembly
            old.unlink()
        readback.build_index(contigs, rb / "contigs", threads)
        steps.mark("readback_index", contig_digest, sorted(rb.glob("contigs.*bt2*")))
    for s in samples:
        r1, r2, n = used[s.name]
        cf = rb / f"{s.name}.counts.tsv"
        rb_key = params_key("readback", contig_digest, file_id(r1), file_id(r2), n, cfg["readback_preset"])
        if not steps.done(f"readback_{s.name}", rb_key):
            steps.invalidate(f"readback_{s.name}")
            tmp = cf.with_name(cf.name + ".part")
            if len(contig_df) and n:
                readback.count_sample(rb / "contigs", r1, r2, tmp, threads, cfg["readback_preset"])
            else:
                tmp.write_text("")
            os.replace(tmp, cf)
            steps.mark(f"readback_{s.name}", rb_key, [cf])
        counts[s.name] = readback.read_counts(cf)
        per_sample[s.name] = readback.category_counts(counts[s.name], cat_of, n)

    # tables, figures, report
    tdir, fdir = outdir / "tables", outdir / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    table = summarize.read_category_table(inputs, per_sample)
    table.to_csv(tdir / "read_category_by_sample.tsv", sep="\t", index=False)
    summarize.contig_count_matrix(contig_df, counts).to_csv(tdir / "contig_counts.tsv", sep="\t", index=False)
    refs = summarize.reference_table(contig_df, counts, {k: v[2] for k, v in used.items()})
    refs.to_csv(tdir / "reference_by_sample.tsv", sep="\t", index=False)
    res_df = pd.DataFrame(resources)
    res_df.to_csv(tdir / "resources.tsv", sep="\t", index=False)
    bar = plots.stacked_bars(table, inputs, fdir / "read_categories")
    sk = [plots.sankey(inputs, table, refs, fdir / "sankey_all_samples.html")]
    (fdir / "sankey").mkdir(exist_ok=True)
    for s in samples:
        sk.append(plots.sankey(inputs, table, refs, fdir / "sankey" / f"{s.name}.html", sample=s.name))
    rep = report.write_report(outdir, cfg, inputs, table, refs, contig_df, bar[0], sk, res_df)
    # publish: the manifest is written last and seals the files later steps read (panel-build verifies it)
    manifest.seal(outdir, manifest.RUN_ARTIFACTS_PRESENT(outdir), {"step": "run", "version": __version__})
    if not cfg.get("keep_work"):
        shutil.rmtree(work / "pool", ignore_errors=True)
    log.info("done: %s", rep)
    return rep
