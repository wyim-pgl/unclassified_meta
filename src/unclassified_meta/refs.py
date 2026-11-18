"""Download and prepare the required references (`unclassified-meta download-refs`).

Layout written under --outdir (default refs/):
  viral/refseq_viral.fna               RefSeq viral genomes (+ BLAST database with the same prefix)
  viral/refseq_viral_protein.dmnd      RefSeq viral proteins, DIAMOND database
  rrna/plant_rrna.fasta                SILVA 138.2 plant/plastid/mitochondrial SSU+LSU + Rfam 5.8S (+ BLAST db)
  protein/sprot_viridiplantae.dmnd     Swiss-Prot reviewed Viridiplantae proteins, DIAMOND database
  REFS_VERSION.tsv                     what was fetched, from where, when, sizes and sequence counts
Every step is skipped when its output exists, so an interrupted run resumes.
"""
from __future__ import annotations

import datetime as dt
import gzip
import logging
import re
import shutil
import urllib.request
from pathlib import Path

from .runner import require, sh

log = logging.getLogger("unclassified_meta")

REFSEQ_VIRAL_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
REFSEQ_RELEASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
REFSEQ_GENOMIC_RE = re.compile(r'href="(viral\.\d+\.1\.genomic\.fna\.gz)"')
REFSEQ_PROTEIN_RE = re.compile(r'href="(viral\.\d+\.protein\.faa\.gz)"')
SILVA_URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/"
SILVA_FILES = ("SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz", "SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz")
SILVA_KEEP = (";Chloroplastida;", ";Chloroplast;", ";Mitochondria;")
RFAM_58S_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF00002.fa.gz"
SPROT_PLANT_URL = ("https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=true"
                   "&query=%28reviewed%3Atrue%29%20AND%20%28taxonomy_id%3A33090%29")
KRAKEN2_INDEX_PAGE = "https://benlangmead.github.io/aws-indexes/k2"
USER_AGENT = "unclassified-meta/0.1 (reference downloader)"


def fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def download(url: str, dest: Path, is_gzip: bool = True) -> Path:
    """Download to dest (via dest.part). Skip when dest exists and, for gzip files, decompresses cleanly."""
    if dest.exists() and dest.stat().st_size > 0 and (not is_gzip or gzip_ok(dest)):
        log.info("have %s", dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    log.info("GET %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as r, open(part, "wb") as fo:
        shutil.copyfileobj(r, fo, 1 << 22)
    if is_gzip and not gzip_ok(part):
        raise RuntimeError(f"downloaded file is not valid gzip: {url}")
    part.rename(dest)
    return dest


def list_refseq_viral() -> tuple[list[str], list[str]]:
    html = fetch_text(REFSEQ_VIRAL_URL)
    genomic = sorted(set(REFSEQ_GENOMIC_RE.findall(html)))
    protein = sorted(set(REFSEQ_PROTEIN_RE.findall(html)))
    if not genomic or not protein:
        raise RuntimeError(f"no viral.*.1.genomic.fna.gz / viral.*.protein.faa.gz listed at {REFSEQ_VIRAL_URL}")
    return genomic, protein


def concat_gz(parts: list[Path], out: Path, transform=None) -> int:
    """Decompress and concatenate FASTA parts into out; return the number of sequences."""
    n = 0
    tmp = out.with_name(out.name + ".part")
    with open(tmp, "w") as fo:
        for p in parts:
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    elif transform:
                        line = transform(line)
                    fo.write(line)
    tmp.rename(out)
    return n


def filter_silva(src: Path, keep=SILVA_KEEP):
    """Yield FASTA lines of SILVA records whose taxonomy path contains one of `keep`, U -> T."""
    keep_rec = False
    with gzip.open(src, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                keep_rec = any(k in line for k in keep)
                if keep_rec:
                    yield line
            elif keep_rec:
                yield line.replace("U", "T").replace("u", "t")


def count_fasta(path: Path) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def plan(outdir: Path) -> list[tuple[str, str, str]]:
    """(reference, url, destination) rows without downloading the data (dry run)."""
    genomic, protein = list_refseq_viral()
    rows = [("refseq_viral_genomic", REFSEQ_VIRAL_URL + f, str(outdir / "viral/src" / f)) for f in genomic]
    rows += [("refseq_viral_protein", REFSEQ_VIRAL_URL + f, str(outdir / "viral/src" / f)) for f in protein]
    rows += [("silva", SILVA_URL + f, str(outdir / "rrna/src" / f)) for f in SILVA_FILES]
    rows += [("rfam_5.8S", RFAM_58S_URL, str(outdir / "rrna/src/RF00002.fa.gz")),
             ("sprot_viridiplantae", SPROT_PLANT_URL, str(outdir / "protein/src/sprot_viridiplantae.faa.gz"))]
    return rows


# Prepared disk use (uncompressed FASTA + BLAST/DIAMOND index) relative to the compressed download,
# measured on a full run (see README "References"); used only for the free-space check.
PREPARED_PER_DOWNLOAD = 3.5


def remote_size(url: str) -> int | None:
    """Content-Length from a HEAD request (None when the server streams without a length)."""
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n and int(n) > 0 else None   # streamed responses report 0 or nothing
    except OSError:
        return None


def size_check(outdir: Path, rows: list[tuple[str, str, str]]) -> tuple[int, int, int]:
    """Print download and prepared-size estimates and warn when free disk space is short.

    Returns (known download bytes, estimated prepared bytes, free bytes)."""
    total, unknown = 0, []
    for name, url, _ in rows:
        n = remote_size(url)
        if n is None:
            unknown.append(name)
        else:
            total += n
    prepared = int(total * PREPARED_PER_DOWNLOAD)
    outdir.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(outdir).free
    gb = 1024 ** 3
    msg = (f"download about {total / gb:.2f} GB"
           + (f" (+ {', '.join(sorted(set(unknown)))} without a size header)" if unknown else "")
           + f"; prepared references about {prepared / gb:.1f} GB; free in {outdir}: {free / gb:.1f} GB")
    log.info(msg)
    print(msg)
    if free < total + prepared:
        warn = (f"WARNING: free space in {outdir} ({free / gb:.1f} GB) is below the expected need "
                f"({(total + prepared) / gb:.1f} GB including temporary downloads).")
        log.warning(warn)
        print(warn)
    return total, prepared, free


def download_refs(outdir: Path, threads: int = 4, rrna_fasta: Path | None = None,
                  keep_sources: bool = False, only: tuple[str, ...] = ("viral", "rrna", "protein")) -> Path:
    require("makeblastdb", "diamond")
    outdir.mkdir(parents=True, exist_ok=True)
    prefixes = {"viral": ("refseq",), "rrna": () if rrna_fasta else ("silva", "rfam"), "protein": ("sprot",)}
    wanted = tuple(p for k in only for p in prefixes[k])
    try:
        size_check(outdir, [r for r in plan(outdir) if r[0].startswith(wanted)] if wanted else [])
    except OSError as e:
        log.warning("could not check remote sizes: %s", e)
    today = dt.date.today().isoformat()
    rows = []

    if "viral" in only:
        vdir = outdir / "viral"
        fna, faa, dmnd = vdir / "refseq_viral.fna", vdir / "refseq_viral_protein.faa", vdir / "refseq_viral_protein.dmnd"
        try:
            release = fetch_text(REFSEQ_RELEASE_URL).strip()
        except OSError:
            release = "unknown"
        if not fna.exists() or not dmnd.exists():
            genomic, protein = list_refseq_viral()
            if not fna.exists():
                parts = [download(REFSEQ_VIRAL_URL + f, vdir / "src" / f) for f in genomic]
                concat_gz(parts, fna)
            if not dmnd.exists():
                parts = [download(REFSEQ_VIRAL_URL + f, vdir / "src" / f) for f in protein]
                concat_gz(parts, faa)
                sh(["diamond", "makedb", "--in", str(faa), "-d", str(dmnd.with_suffix("")), "-p", str(threads)])
        if not Path(str(fna) + ".ndone").exists():
            sh(["makeblastdb", "-in", str(fna), "-dbtype", "nucl", "-out", str(fna)])
            Path(str(fna) + ".ndone").touch()
        rows.append(("refseq_viral_genomic", f"RefSeq release {release}", REFSEQ_VIRAL_URL + "viral.*.1.genomic.fna.gz",
                     str(fna), count_fasta(fna), fna.stat().st_size))
        rows.append(("refseq_viral_protein", f"RefSeq release {release}", REFSEQ_VIRAL_URL + "viral.*.protein.faa.gz",
                     str(dmnd), count_fasta(faa) if faa.exists() else "", dmnd.stat().st_size))

    if "rrna" in only:
        rdir = outdir / "rrna"
        rdir.mkdir(parents=True, exist_ok=True)
        out = rdir / "plant_rrna.fasta"
        if rrna_fasta:
            if not out.exists():
                shutil.copyfile(rrna_fasta, out)
            source = f"copied from {rrna_fasta}"
        else:
            source = SILVA_URL + " + " + RFAM_58S_URL
            if not out.exists():
                srcs = [download(SILVA_URL + f, rdir / "src" / f) for f in SILVA_FILES]
                rfam = download(RFAM_58S_URL, rdir / "src" / "RF00002.fa.gz")
                tmp = out.with_name(out.name + ".part")
                with open(tmp, "w") as fo:
                    for s in srcs:
                        fo.writelines(filter_silva(s))
                    with gzip.open(rfam, "rt") as fh:
                        for line in fh:
                            fo.write(line if line.startswith(">") else line.replace("U", "T").replace("u", "t"))
                tmp.rename(out)
        if not Path(str(out) + ".ndone").exists():
            sh(["makeblastdb", "-in", str(out), "-dbtype", "nucl", "-out", str(out)])
            Path(str(out) + ".ndone").touch()
        rows.append(("plant_rrna", "SILVA 138.2 NR99 (Chloroplastida, Chloroplast, Mitochondria) + Rfam RF00002",
                     source, str(out), count_fasta(out), out.stat().st_size))

    if "protein" in only:
        pdir = outdir / "protein"
        faa, dmnd = pdir / "sprot_viridiplantae.faa", pdir / "sprot_viridiplantae.dmnd"
        if not dmnd.exists():
            gz = download(SPROT_PLANT_URL, pdir / "src" / "sprot_viridiplantae.faa.gz")
            concat_gz([gz], faa)
            sh(["diamond", "makedb", "--in", str(faa), "-d", str(dmnd.with_suffix("")), "-p", str(threads)])
        rows.append(("sprot_viridiplantae", f"UniProtKB reviewed, taxonomy 33090, fetched {today}", SPROT_PLANT_URL,
                     str(dmnd), count_fasta(faa) if faa.exists() else "", dmnd.stat().st_size))

    if not keep_sources:
        for sub in ("viral/src", "rrna/src", "protein/src"):
            shutil.rmtree(outdir / sub, ignore_errors=True)

    vt = outdir / "REFS_VERSION.tsv"
    old = {}
    if vt.exists():
        for line in open(vt).read().splitlines()[1:]:
            old[line.split("\t")[0]] = line
    for r in rows:
        old[r[0]] = "\t".join(map(str, (r[0], r[1], r[2], r[3], r[4], r[5], today)))
    with open(vt, "w") as fo:
        fo.write("reference\tversion\tsource\tfile\tn_sequences\tbytes\tdate\n")
        fo.write("\n".join(old.values()) + "\n")
    return vt
