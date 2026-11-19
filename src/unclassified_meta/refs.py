"""Download and prepare the required references (`unclassified-meta download-refs`).

Layout written under --outdir (default refs/):
  viral/refseq_viral.fna               RefSeq viral genomes (+ BLAST database with the same prefix)
  viral/refseq_viral_protein.dmnd      RefSeq viral proteins, DIAMOND database
  rrna/plant_rrna.fasta                SILVA 138.2 plant/plastid/mitochondrial SSU+LSU + Rfam 5.8S (+ BLAST db)
  protein/sprot_viridiplantae.dmnd     Swiss-Prot reviewed Viridiplantae proteins, DIAMOND database
  REFS_VERSION.tsv                     what was fetched, from where, when, sizes and sequence counts
A folder is reused only as a sealed generation (MANIFEST.json), and a downloaded part only when its sidecar
(<part>.release.json) names the same URL and release and its checksum still matches, so an interrupted run
resumes without ever mixing releases. A folder prepared by hand or copied from elsewhere is adopted with
`unclassified-meta seal-refs <dir>` (see seal_refs).
"""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import urllib.request
from pathlib import Path

from . import manifest
from .io_utils import refuse_overlap, staged
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


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def upstream_md5(url: str) -> str | None:
    """The md5 the provider publishes next to a file (<url>.md5, as SILVA does), or None when there is none."""
    try:
        m = re.search(r"\b([0-9a-f]{32})\b", fetch_text(url + ".md5"))
    except OSError:
        m = None
    if m is None:
        log.info("no upstream checksum for %s", url)
    return m.group(1) if m else None


def download(url: str, dest: Path, is_gzip: bool = True, release: str | None = None,
             md5: str | None = None) -> Path:
    """Download url to dest through a unique temporary file. A cached dest is reused only when its sidecar
    (dest.release.json) records the same URL and release and dest still has the recorded sha256; a source
    without a release number is identified by the day it was fetched. With `md5` (the provider's checksum)
    the download is verified against it."""
    dest = Path(dest)
    identity = {"url": url, "release": release or f"fetched {dt.date.today().isoformat()}"}
    if md5:
        identity["md5"] = md5           # a cached file is reused only if it was verified against this checksum
    side = dest.with_name(dest.name + ".release.json")
    if dest.exists() and side.exists():
        try:
            rec = json.loads(side.read_text())
        except ValueError:
            rec = {}
        if {k: rec.get(k) for k in identity} == identity and rec.get("sha256") == _file_sha256(dest) and \
                (not is_gzip or gzip_ok(dest)):
            log.info("have %s (%s)", dest, identity["release"])
            return dest
        log.info("cached %s belongs to another release or changed; fetching again", dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    log.info("GET %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    side.unlink(missing_ok=True)                       # the cache record goes first, the file follows
    with staged(dest) as part:
        with urllib.request.urlopen(req, timeout=600) as r, open(part, "wb") as fo:
            shutil.copyfileobj(r, fo, 1 << 22)
        if is_gzip and not gzip_ok(part):
            raise RuntimeError(f"downloaded file is not valid gzip: {url}")
        if md5:
            h = hashlib.md5()
            with open(part, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(chunk)
            if h.hexdigest() != md5:
                raise RuntimeError(f"checksum mismatch for {url}: md5 {h.hexdigest()}, provider says {md5}")
    with staged(side) as tmp:
        tmp.write_text(json.dumps({**identity, "sha256": _file_sha256(dest)}) + "\n")
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
    with staged(out) as tmp, open(tmp, "w") as fo:
        for p in parts:
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    elif transform:
                        line = transform(line)
                    fo.write(line)
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
PREPARED_PER_DOWNLOAD = 2.5   # 1.3 GB prepared / 0.53 GB downloaded, 2026-11-17


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


def bind_blast_db(fasta: Path) -> None:
    """The BLAST database next to a FASTA, bound to its content by a sidecar (see classify.build_blast_db)."""
    from .classify import blast_db_complete, build_blast_db, sidecar_matches
    if not (blast_db_complete(fasta) and sidecar_matches(fasta, [fasta])):
        build_blast_db(fasta, fasta)


def db_files(fasta: Path) -> list[str]:
    return sorted(f.name for f in fasta.parent.glob(fasta.name + ".n*")) + [fasta.name + ".source.json"]


def _makedb_diamond(faa: Path, dmnd: Path, threads: int) -> None:
    """DIAMOND database built under a unique temporary stem and renamed into place."""
    import tempfile
    stem = Path(tempfile.mkdtemp(prefix=".diamond-", dir=dmnd.parent)) / dmnd.stem
    try:
        sh(["diamond", "makedb", "--in", str(faa), "-d", str(stem), "-p", str(threads)])
        os.replace(str(stem) + ".dmnd", dmnd)
    finally:
        shutil.rmtree(stem.parent, ignore_errors=True)


def write_versions(outdir: Path, rows, fresh: set, today: str, adopted: dict | None = None) -> Path:
    """REFS_VERSION.tsv: fresh rows get today's date; reused ones keep their record, or 'unknown' without one.
    `adopted` (seal-refs) maps a reference name to the version the user gave, or 'unknown'."""
    vt = outdir / "REFS_VERSION.tsv"
    old = {}
    if vt.exists():
        for line in open(vt).read().splitlines()[1:]:
            old[line.split("\t")[0]] = line
    for r in rows:
        if adopted is not None:
            old[r[0]] = "\t".join(map(str, (r[0], adopted.get(r[0], "unknown"), "adopted with seal-refs", r[3], r[4],
                                            r[5], "unknown")))
        elif r[0] in fresh:
            old[r[0]] = "\t".join(map(str, (r[0], r[1], r[2], r[3], r[4], r[5], today)))
        elif r[0] not in old:
            # reused, but nothing records where and when it was fetched: say so instead of guessing
            old[r[0]] = "\t".join(map(str, (r[0], "unknown", r[2], r[3], r[4], r[5], "unknown")))
        # else: reused database with a record: keep the version and date it was downloaded with
    with staged(vt) as tmp:
        tmp.write_text("reference\tversion\tsource\tfile\tn_sequences\tbytes\tdate\n" + "\n".join(old.values()) + "\n")
    return vt


SEAL_LAYOUT = {   # folder -> (FASTA files that must hold sequences, DIAMOND db and the FASTA it is built from)
    "viral": (("refseq_viral.fna",), ("refseq_viral_protein.dmnd", "refseq_viral_protein.faa")),
    "rrna": (("plant_rrna.fasta",), None),
    "protein": ((), ("sprot_viridiplantae.dmnd", "sprot_viridiplantae.faa")),
}
SEAL_NAMES = {"viral": ("refseq_viral_genomic", "refseq_viral_protein"), "rrna": ("plant_rrna",),
              "protein": ("sprot_viridiplantae",)}


def seal_refs(refs_dir: Path, versions: dict | None = None, threads: int = 4) -> Path:
    """Adopt a reference folder made by hand or copied from another machine (the download-refs layout): check
    that every expected FASTA holds sequences, build missing BLAST / DIAMOND databases, and write the manifest
    of each folder last. Versions are 'unknown' unless given (NAME=VERSION)."""
    refs_dir, versions = Path(refs_dir), dict(versions or {})
    unknown = set(versions) - {n for v in SEAL_NAMES.values() for n in v}
    if unknown:
        raise ValueError(f"unknown reference name(s) {', '.join(sorted(unknown))}; use "
                         f"{', '.join(n for v in SEAL_NAMES.values() for n in v)}")
    rows = []
    for sub, (fastas, dia) in SEAL_LAYOUT.items():
        d = refs_dir / sub
        if not d.is_dir():
            raise FileNotFoundError(f"{d} not found: seal-refs expects the download-refs layout (viral/, rrna/, protein/)")
        manifest.begin(d, f"seal-refs {sub}")
        files = []
        for name in fastas:
            f = d / name
            if not f.is_file() or count_fasta(f) == 0:
                raise ValueError(f"{f} is missing or holds no FASTA sequence")
            bind_blast_db(f)
            files += [name] + db_files(f)
        if dia:
            dmnd, faa = d / dia[0], d / dia[1]
            if not dmnd.exists():
                if not faa.is_file() or count_fasta(faa) == 0:
                    raise ValueError(f"{dmnd} is missing and {faa} holds no protein sequence to build it from")
                require("diamond")
                _makedb_diamond(faa, dmnd, threads)
            files.append(dmnd.name)
        manifest.seal(d, files, {"step": "seal-refs", "versions": {n: versions.get(n, "unknown")
                                                                   for n in SEAL_NAMES[sub]}})
        # one row per reference: (file it is used from, FASTA its sequences are counted in)
        used = [(d / n, d / n) for n in fastas] + ([(d / dia[0], d / dia[1])] if dia else [])
        for name, (f, src) in zip(SEAL_NAMES[sub], used):
            rows.append((name, "", "", str(f), count_fasta(src) if src.exists() else "", f.stat().st_size))
    return write_versions(refs_dir, rows, set(), "", adopted=versions)


def download_refs(outdir: Path, threads: int = 4, rrna_fasta: Path | None = None,
                  keep_sources: bool = False, only: tuple[str, ...] = ("viral", "rrna", "protein"),
                  allow_unverified_silva: bool = False) -> Path:
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
    fresh: set[str] = set()          # references (re)built in this call; the others keep their recorded version

    def complete(folder: Path, files: list[str]) -> bool:
        """A folder is reused only as a complete, sealed generation (manifest written after its last file)."""
        try:
            manifest.verify(folder, required=files, what="reference files")
            return True
        except manifest.ManifestError as e:
            if (folder / manifest.NAME).exists() or any((folder / f).exists() for f in files):
                log.warning("%s; rebuilding", e)
            return False

    def makedb_diamond(faa: Path, dmnd: Path) -> None:
        _makedb_diamond(faa, dmnd, threads)

    if "viral" in only:
        vdir = outdir / "viral"
        vdir.mkdir(parents=True, exist_ok=True)
        fna, faa, dmnd = vdir / "refseq_viral.fna", vdir / "refseq_viral_protein.faa", vdir / "refseq_viral_protein.dmnd"
        release = None
        if not complete(vdir, [fna.name, dmnd.name]):
            manifest.begin(vdir, "download-refs viral")
            try:
                release = fetch_text(REFSEQ_RELEASE_URL).strip()
            except OSError:
                release = "unknown"
            # cached parts are reused only for the same release (an unknown release: only on the same day)
            rel = None if release == "unknown" else f"RefSeq {release}"
            genomic, protein = list_refseq_viral()
            concat_gz([download(REFSEQ_VIRAL_URL + f, vdir / "src" / f, release=rel) for f in genomic], fna)
            concat_gz([download(REFSEQ_VIRAL_URL + f, vdir / "src" / f, release=rel) for f in protein], faa)
            makedb_diamond(faa, dmnd)
            fresh |= {"refseq_viral_genomic", "refseq_viral_protein"}
        bind_blast_db(fna)
        manifest.seal(vdir, [fna.name, dmnd.name] + db_files(fna), {"step": "download-refs viral"})
        rows.append(("refseq_viral_genomic", f"RefSeq release {release}", REFSEQ_VIRAL_URL + "viral.*.1.genomic.fna.gz",
                     str(fna), count_fasta(fna), fna.stat().st_size))
        rows.append(("refseq_viral_protein", f"RefSeq release {release}", REFSEQ_VIRAL_URL + "viral.*.protein.faa.gz",
                     str(dmnd), count_fasta(faa) if faa.exists() else "", dmnd.stat().st_size))

    if "rrna" in only:
        rdir = outdir / "rrna"
        rdir.mkdir(parents=True, exist_ok=True)
        out = rdir / "plant_rrna.fasta"
        source = f"copied from {rrna_fasta}" if rrna_fasta else SILVA_URL + " + " + RFAM_58S_URL
        if rrna_fasta:
            refuse_overlap([out], [rrna_fasta], "download-refs")
        # SILVA verification status travels with the generation: a cached folder keeps the status its manifest
        # records; it becomes "verified" only by a rebuild from downloads checked against the published MD5
        cached = complete(rdir, [out.name])
        status = None
        if cached:
            try:
                rec = json.loads((rdir / manifest.NAME).read_text())
            except (OSError, ValueError):
                rec = {}
            status = rec.get("silva_md5", "unknown")
            if status not in ("verified", "not applicable") and not rrna_fasta and not allow_unverified_silva:
                missing = [f for f in SILVA_FILES if upstream_md5(SILVA_URL + f) is None]
                if missing:
                    raise RuntimeError(f"the SILVA files in {rdir} are recorded as {status} (not checked against "
                                       "the checksum SILVA publishes), and the checksum still cannot be fetched "
                                       f"({SILVA_URL}<file>.md5). Retry when arb-silva.de is reachable (the files "
                                       "are then downloaded again and checked), or pass --allow-unverified-silva "
                                       "to keep using them as unverified")
                log.warning("SILVA files in %s are %s: downloading them again and checking the MD5", rdir, status)
                cached = False
        if not cached:
            status = "not applicable" if rrna_fasta else "verified"
            if not rrna_fasta:
                # SILVA publishes an MD5 next to every file; without it a download cannot be verified
                silva_md5 = {f: upstream_md5(SILVA_URL + f) for f in SILVA_FILES}
                missing = [f for f, m in silva_md5.items() if m is None]
                if missing and not allow_unverified_silva:
                    raise RuntimeError(f"could not fetch the SILVA checksum(s) for {', '.join(missing)} "
                                       f"({SILVA_URL}<file>.md5), so the download cannot be verified. Retry when "
                                       "arb-silva.de is reachable, or pass --allow-unverified-silva to use the files "
                                       "without a checksum (recorded as unverified in REFS_VERSION.tsv)")
                if missing:
                    status = "unverified"
            manifest.begin(rdir, "download-refs rrna")
            with staged(out) as tmp:
                if rrna_fasta:
                    shutil.copyfile(rrna_fasta, tmp)
                else:
                    srcs = [download(SILVA_URL + f, rdir / "src" / f, release="SILVA 138.2", md5=silva_md5[f])
                            for f in SILVA_FILES]
                    rfam = download(RFAM_58S_URL, rdir / "src" / "RF00002.fa.gz")
                    with open(tmp, "w") as fo:
                        for s_ in srcs:
                            fo.writelines(filter_silva(s_))
                        with gzip.open(rfam, "rt") as fh:
                            for line in fh:
                                fo.write(line if line.startswith(">") else line.replace("U", "T").replace("u", "t"))
            fresh.add("plant_rrna")
        bind_blast_db(out)
        manifest.seal(rdir, [out.name] + db_files(out), {"step": "download-refs rrna", "silva_md5": status})
        rows.append(("plant_rrna", "SILVA 138.2 NR99 (Chloroplastida, Chloroplast, Mitochondria) + Rfam RF00002"
                     + (" (SILVA MD5 unverified: --allow-unverified-silva)" if status == "unverified" else ""),
                     source, str(out), count_fasta(out), out.stat().st_size))

    if "protein" in only:
        pdir = outdir / "protein"
        pdir.mkdir(parents=True, exist_ok=True)
        faa, dmnd = pdir / "sprot_viridiplantae.faa", pdir / "sprot_viridiplantae.dmnd"
        if not complete(pdir, [dmnd.name]):
            manifest.begin(pdir, "download-refs protein")
            gz = download(SPROT_PLANT_URL, pdir / "src" / "sprot_viridiplantae.faa.gz")
            concat_gz([gz], faa)
            makedb_diamond(faa, dmnd)
            fresh.add("sprot_viridiplantae")
        manifest.seal(pdir, [dmnd.name], {"step": "download-refs protein"})
        rows.append(("sprot_viridiplantae", f"UniProtKB reviewed, taxonomy 33090, fetched {today}", SPROT_PLANT_URL,
                     str(dmnd), count_fasta(faa) if faa.exists() else "", dmnd.stat().st_size))

    if not keep_sources:
        for sub in ("viral/src", "rrna/src", "protein/src"):
            shutil.rmtree(outdir / sub, ignore_errors=True)

    return write_versions(outdir, rows, fresh, today)
