"""The pipeline: extract -> subsample -> pool -> assemble -> classify contigs -> read-back -> tables/figures/report."""
from __future__ import annotations

import json
import logging
import os
import shlex
import shutil
from pathlib import Path

import pandas as pd

from . import assemble as asm
from . import classify, plots, readback, report, summarize
from .config import Sample, check_references, write_samplesheet
from .extract import filter_fastq, parse_classifier
from .io_utils import count_fastq
from .runner import Steps, require, setup_logging, sh

log = logging.getLogger("unclassified_meta")


def q(x) -> str:
    return shlex.quote(str(x))


def out_reads(outdir: Path, s: Sample) -> tuple[Path, Path | None]:
    d = outdir / "unclassified_reads"
    if s.paired:
        return d / f"{s.name}_R1.fastq.gz", d / f"{s.name}_R2.fastq.gz"
    return d / f"{s.name}.fastq.gz", None


def extract_sample(s: Sample, outdir: Path) -> dict:
    """Write the sample's unclassified reads and return its input summary."""
    o1, o2 = out_reads(outdir, s)
    o1.parent.mkdir(parents=True, exist_ok=True)
    if s.classifier == "none":
        for src, dst in ((s.r1, o1), (s.r2, o2)):
            if src is None:
                continue
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            if str(src).endswith(".gz"):
                os.symlink(Path(src).resolve(), dst)
            else:
                sh(f"gzip -1 -c {q(src)} > {q(dst)}")
        n = count_fastq(o1)
        return {"sample": s.name, "paired": s.paired, "classifier": "none", "total_reads": n,
                "classified_reads": 0, "unclassified_reads": n}
    ids, total, classified = parse_classifier(s.classifier_output, s.classifier)
    n1 = filter_fastq(s.r1, o1, ids)
    if s.r2:
        n2 = filter_fastq(s.r2, o2, ids)
        if n2 != n1:
            raise ValueError(f"{s.name}: {n1} mate-1 vs {n2} mate-2 unclassified reads; are R1/R2 in the same order?")
    if n1 != len(ids):
        log.warning("%s: classifier lists %d unclassified reads but %d were found in %s",
                    s.name, len(ids), n1, s.r1)
    return {"sample": s.name, "paired": s.paired, "classifier": s.classifier, "total_reads": total,
            "classified_reads": classified, "unclassified_reads": n1}


def extract_all(samples: list[Sample], outdir: Path) -> pd.DataFrame:
    work = outdir / "work"
    steps = Steps(work)
    rows = []
    for s in samples:
        js = work / "inputs" / f"{s.name}.json"
        js.parent.mkdir(parents=True, exist_ok=True)
        if not steps.done(f"extract_{s.name}"):
            log.info("extract unclassified reads: %s (%s)", s.name, s.classifier)
            js.write_text(json.dumps(extract_sample(s, outdir)))
            steps.mark(f"extract_{s.name}")
        rows.append(json.loads(js.read_text()))
    df = pd.DataFrame(rows)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / "tables" / "input_summary.tsv", sep="\t", index=False)
    return df


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
    if not done.exists():
        require("seqtk")
        sh(f"seqtk sample -s {seed} {q(r1)} {n} | gzip -1 > {q(o1)}")
        if r2:
            sh(f"seqtk sample -s {seed} {q(r2)} {n} | gzip -1 > {q(o2)}")
        done.write_text(str(count_fastq(o1)))
    return o1, o2, int(done.read_text())


def run(cfg: dict, samples: list[Sample]) -> Path:
    outdir = Path(cfg["outdir"]).resolve()
    work = outdir / "work"
    setup_logging(outdir / "logs" / "pipeline.log")
    check_references(cfg)
    require("seqtk", "minimap2", "blastn", "makeblastdb", "diamond", "bowtie2", "bowtie2-build", "samtools")
    steps = Steps(work)
    write_samplesheet(samples, outdir / "samplesheet.used.tsv")
    threads = int(cfg["threads"])

    inputs = extract_all(samples, outdir)

    used = {}
    for s in samples:
        used[s.name] = subsample(s, outdir, cfg["subsample"], int(cfg["seed"]))
    inputs["reads_used"] = inputs["sample"].map({k: v[2] for k, v in used.items()})
    inputs.to_csv(outdir / "tables" / "input_summary.tsv", sep="\t", index=False)

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
    if not steps.done("pool"):
        if paired:
            sh("cat " + " ".join(q(v[0]) for v in paired) + f" > {q(pool['r1'])}")
            sh("cat " + " ".join(q(v[1]) for v in paired) + f" > {q(pool['r2'])}")
        if single:
            sh("cat " + " ".join(q(v[0]) for v in single) + f" > {q(pool['s'])}")
        steps.mark("pool")
    n_pool = sum(v[2] for v in used.values())
    if n_pool == 0:
        raise SystemExit("no unclassified reads to assemble")

    # assemble
    assembler = asm.choose_assembler(cfg)
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
    if not steps.done("assemble"):
        fasta, res = asm.assemble(assembler, pool, work / "assembly", threads, float(cfg["memory_gb"]))
        contigs.parent.mkdir(parents=True, exist_ok=True)
        n = asm.filter_and_rename(fasta, contigs, outdir / "tables" / "contig_names.tsv", int(cfg["min_contig_len"]))
        resources.append({"step": f"assembly ({assembler})", "input_reads": n_pool, "estimated_gb": est,
                          **res, "contigs_kept": n})
        res_file.write_text(json.dumps(resources))
        steps.mark("assemble")
    log.info("contigs >= %d bp: %s", cfg["min_contig_len"], contigs)

    # classify
    ctab = outdir / "tables" / "contig_category.tsv"
    if not steps.done("classify"):
        result = classify.Classifier(cfg, work / "classify", threads).run(contigs) if contigs.stat().st_size else {}
        classify.write_contig_table(result, ctab)
        steps.mark("classify")
    contig_df = pd.read_csv(ctab, sep="\t", dtype={"best_hit": str, "best_hit_title": str})
    cat_of = dict(zip(contig_df["contig"], contig_df["category"]))

    # read-back
    rb = work / "readback"
    rb.mkdir(parents=True, exist_ok=True)
    counts, per_sample = {}, {}
    if len(contig_df):
        readback.build_index(contigs, rb / "contigs", threads)
    for s in samples:
        r1, r2, n = used[s.name]
        cf = rb / f"{s.name}.counts.tsv"
        if not steps.done(f"readback_{s.name}"):
            if len(contig_df) and n:
                readback.count_sample(rb / "contigs", r1, r2, cf, threads, cfg["readback_preset"])
            else:
                cf.write_text("")
            steps.mark(f"readback_{s.name}")
        counts[s.name] = readback.read_counts(cf)
        per_sample[s.name] = readback.category_counts(counts[s.name], cat_of, n)

    # tables, figures, report
    tdir, fdir = outdir / "tables", outdir / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    table = summarize.read_category_table(inputs, per_sample)
    table.to_csv(tdir / "read_category_by_sample.tsv", sep="\t", index=False)
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
    if not cfg.get("keep_work"):
        shutil.rmtree(work / "pool", ignore_errors=True)
    log.info("done: %s", rep)
    return rep
