"""Read evidence behind presence and structural zeros, per feature and group level (`unclassified-meta detection`).

ANCOM-BC2 calls a feature a structural zero in a group when every sample of the group has count == 0 (presence is
count != 0). With EM expected counts (taxa-table --assign em) a fraction of a shared read is a nonzero count, so the
zero pattern depends on how the shared reads were divided. This module sets the evidence side by side, with the
same rule for every source and no threshold on any value:

  lca           counts.tsv of a taxa-table --assign lca folder (optional, --lca): whole reads
  unique        unique_counts.tsv of the --assign em folder: whole reads whose candidate set is that taxon alone
  em_expected   counts.tsv of the --assign em folder: expected reads, fractions as written
  imputations   imputations/counts_NNN.tsv of the --assign em folder: whole reads drawn from the posterior

Outputs (sealed with MANIFEST.json, written last):

  detection_by_sample.tsv  feature_id, label, sample, level, [lca_count], unique_count, em_expected,
                           [imputation_mean, imputation_nonzero_frac (share of imputations with count > 0)]
  detection_by_group.tsv   per feature x level: n_samples; per source <source>_present_samples (count != 0) and
                           <source>_reads (sum); [imputation_all_zero_frac (share of imputations in which every
                           sample of the level is 0)]; struc_zero_<source> (every sample of the level is 0)
  detection_summary.tsv    per feature: struc_zero_levels_<source> (levels that are structural zeros under that
                           source), [imputation_all_zero_frac (<group> = L) per level, inconsistent (some level is a
                           structural zero in some imputations and not in others)], interpretation

`level` holds the values of the metadata column given with --group. Columns in [] are written only when their
source exists (--lca given; imputations sealed in the EM folder).

interpretation: each level gets the first matching class below (u = unique reads, e = EM expected count,
l = LCA reads, each summed over the samples of the level; l only with --lca):

  zero          u == 0, e == 0 and (no --lca or l == 0)   "zero in every source in <levels>"
  shared only   u == 0 and e != 0                         "no unique-read support in <levels>; EM expected count
                                                           there comes only from shared reads"
  lca only      u == 0, e == 0 and l != 0                 "LCA reads but no unique or EM expected reads in <levels>"
  lca zero      u != 0 and l == 0 (with --lca)            "LCA zero but unique reads present in <levels>"
  unique        u != 0                                    "unique-read support in <levels>"

When every level falls in the same class the text says "in every level" ("zero in every source" for the zero class);
otherwise the statements of the classes present are joined with " | " in the order above.
"""
from __future__ import annotations

import hashlib
import io
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import manifest
from .io_utils import staged

log = logging.getLogger("unclassified_meta")
OUTPUTS = ("detection_by_sample.tsv", "detection_by_group.tsv", "detection_summary.tsv")

_TEXT = {
    "zero": "zero in every source in {}",
    "shared": "no unique-read support in {}; EM expected count there comes only from shared reads",
    "lca_only": "LCA reads but no unique or EM expected reads in {}",
    "lca_zero": "LCA zero but unique reads present in {}",
    "unique": "unique-read support in {}",
}


def _sealed(folder: Path, name: str, rec: dict) -> bytes:
    """The bytes of a sealed file, checked against the sha256 its manifest recorded (read once, then parsed)."""
    want = (rec.get("artifacts") or {}).get(name)
    if want is None:
        raise manifest.ManifestError(f"{folder / name} is not covered by {folder / manifest.NAME}")
    data = (folder / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != want:
        raise manifest.ManifestError(f"{folder / name} changed after it was sealed; run taxa-table again")
    return data


def _counts(data: bytes, what: str) -> pd.DataFrame:
    df = pd.read_csv(io.BytesIO(data), sep="\t", dtype=str, keep_default_na=False).set_index("feature_id")
    try:
        return df.apply(pd.to_numeric)
    except ValueError as e:
        raise ValueError(f"{what}: a count is not a number ({e})") from None


def _labels(data: bytes) -> dict:
    df = pd.read_csv(io.BytesIO(data), sep="\t", dtype=str, keep_default_na=False)
    return dict(zip(df["feature_id"], df["label"])) if "label" in df.columns else {}


def load_em(em_dir: str | Path) -> dict:
    """counts, unique counts, labels and imputation draws of a sealed taxa-table --assign em folder."""
    from .taxa import read_draw
    em_dir = Path(em_dir)
    rec = manifest.verify(em_dir, what="taxa tables")
    if rec.get("step") != "taxa-table" or rec.get("assignment") != "em":
        raise ValueError(f"--em {em_dir}: not a taxa-table --per-read --assign em folder (its manifest records "
                         f"step {rec.get('step')!r}, assignment {rec.get('assignment')!r})")
    if "unique_counts.tsv" not in rec["artifacts"]:
        raise ValueError(f"--em {em_dir} holds no unique_counts.tsv (made by an older version); run taxa-table "
                         "--assign em again")
    counts = _counts(_sealed(em_dir, "counts.tsv", rec), f"{em_dir}/counts.tsv")
    unique = _counts(_sealed(em_dir, "unique_counts.tsv", rec), f"{em_dir}/unique_counts.tsv")
    if list(unique.index) != list(counts.index) or list(unique.columns) != list(counts.columns):
        raise ValueError(f"{em_dir}: unique_counts.tsv and counts.tsv differ in rows or samples")
    draws = []
    for m in range(int(rec.get("imputations") or 0)):
        name = f"imputations/counts_{m + 1:03d}.tsv"
        d = read_draw(io.BytesIO(_sealed(em_dir, name, rec)))
        if list(d.index) != list(counts.index) or list(d.columns) != list(counts.columns):
            raise ValueError(f"{em_dir}/{name}: rows or samples differ from counts.tsv")
        draws.append(d)
    labels = _labels(_sealed(em_dir, "features.tsv", rec)) if "features.tsv" in rec["artifacts"] else {}
    return {"rec": rec, "counts": counts, "unique": unique, "draws": draws, "labels": labels}


def load_lca(lca_dir: str | Path, em_rec: dict) -> dict:
    """counts and labels of a sealed taxa-table folder made without EM, at the rank and namespace of the EM one."""
    lca_dir = Path(lca_dir)
    rec = manifest.verify(lca_dir, what="taxa tables")
    if rec.get("step") != "taxa-table" or rec.get("assignment") == "em":
        raise ValueError(f"--lca {lca_dir}: not a taxa-table folder made with --assign lca (its manifest records "
                         f"step {rec.get('step')!r}, assignment {rec.get('assignment')!r})")
    for key in ("rank", "namespace"):
        if rec.get(key) != em_rec.get(key):
            raise ValueError(f"--lca and --em were made with different {key}s ({rec.get(key)!r} and "
                             f"{em_rec.get(key)!r}); their feature IDs cannot be matched")
    counts = _counts(_sealed(lca_dir, "counts.tsv", rec), f"{lca_dir}/counts.tsv")
    labels = _labels(_sealed(lca_dir, "features.tsv", rec)) if "features.tsv" in rec["artifacts"] else {}
    return {"rec": rec, "counts": counts, "labels": labels}


def _samples(metadata: str | Path, group: str, have: list, lca: list | None) -> tuple[list, list]:
    meta = pd.read_csv(metadata, sep="\t", dtype=str, keep_default_na=False)
    for col in ("sample", group):
        if col not in meta.columns:
            raise ValueError(f"{metadata} has no column {col!r}")
    if meta["sample"].duplicated().any():
        raise ValueError(f"{metadata}: duplicate sample names")
    empty = meta[group].str.strip().isin(["", "NA"])
    if empty.any():
        log.warning("detection: %d sample(s) without a value in %s were left out: %s", int(empty.sum()), group,
                    ", ".join(meta.loc[empty, "sample"]))
        meta = meta[~empty]
    missing = [s for s in meta["sample"] if s not in have]
    if missing:
        raise ValueError(f"samples in the metadata without counts in --em: {', '.join(missing[:5])}")
    if lca is not None:
        missing = [s for s in meta["sample"] if s not in lca]
        if missing:
            raise ValueError(f"samples in the metadata without counts in --lca: {', '.join(missing[:5])}")
    extra = [s for s in have if s not in set(meta["sample"])]
    if extra:
        log.warning("detection: samples with counts but not in the metadata were left out: %s", ", ".join(extra))
    if meta.empty:
        raise ValueError("no sample of the metadata has a group value")
    return list(meta["sample"]), list(meta[group])


def _pick(features: list, labels: dict, wanted: list) -> list:
    if not wanted:
        return features
    out = []
    for w in wanted:
        hit = [w] if w in features else [f for f in features if labels.get(f) == w]
        if not hit:
            raise ValueError(f"--feature {w}: no feature with this ID or label in the tables")
        out += [f for f in hit if f not in out]
    return out


def _classify(u: float, e: float, lca: float | None) -> str:
    if u == 0 and e == 0:
        return "lca_only" if lca is not None and lca != 0 else "zero"
    if u == 0:
        return "shared"
    if lca is not None and lca == 0:
        return "lca_zero"
    return "unique"


def interpretation(classes: dict) -> str:
    """The summary text of one feature from {level: class} (classes of `_classify`)."""
    if len(set(classes.values())) == 1:
        c = next(iter(classes.values()))
        return "zero in every source" if c == "zero" else _TEXT[c].format("every level")
    parts = []
    for c in _TEXT:
        lv = [str(k) for k, v in classes.items() if v == c]
        if lv:
            parts.append(_TEXT[c].format(", ".join(lv)))
    return " | ".join(parts)


def detection_tables(em: dict, lca: dict | None, samples: list, groups: list, features: list | None = None,
                     group: str = "level") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list]:
    """The three tables (by sample, by group, summary) and the group levels."""
    ids = list(em["counts"].index)
    if lca is not None:
        ids += [f for f in lca["counts"].index if f not in set(ids)]
    labels = {**(lca["labels"] if lca else {}), **em["labels"]}
    feats = _pick(ids, labels, list(features or []))
    take = lambda df: df.reindex(index=feats, fill_value=0)[samples].to_numpy(dtype=float)  # noqa: E731
    src = {"unique": take(em["unique"]), "em_expected": take(em["counts"])}
    if lca is not None:
        src = {"lca": take(lca["counts"]), **src}
    draws = np.stack([take(d) for d in em["draws"]]) if em["draws"] else None
    levels = sorted(set(groups))
    g = np.asarray(groups)
    nf, ns = len(feats), len(samples)
    lab = [labels.get(f, "") for f in feats]

    by_sample = pd.DataFrame({"feature_id": np.repeat(feats, ns), "label": np.repeat(lab, ns),
                              "sample": np.tile(samples, nf), "level": np.tile(groups, nf)})
    names = {"lca": "lca_count", "unique": "unique_count", "em_expected": "em_expected"}
    for k, x in src.items():
        by_sample[names[k]] = x.ravel()
    if draws is not None:
        by_sample["imputation_mean"] = draws.mean(axis=0).ravel()
        by_sample["imputation_nonzero_frac"] = (draws > 0).mean(axis=0).ravel()

    rows, cls = [], {f: {} for f in feats}
    zero_frac = {}
    for lv in levels:
        idx = np.flatnonzero(g == lv)
        part = pd.DataFrame({"feature_id": feats, "label": lab, "level": lv, "n_samples": len(idx)})
        sums = {}
        for k, x in src.items():
            part[f"{k}_present_samples"] = (x[:, idx] != 0).sum(axis=1)
            part[f"{k}_reads"] = sums[k] = x[:, idx].sum(axis=1)
        if draws is not None:
            part["imputation_all_zero_frac"] = zero_frac[lv] = (draws[:, :, idx] == 0).all(axis=2).mean(axis=0)
        for k in src:
            part[f"struc_zero_{k}"] = part[f"{k}_present_samples"] == 0
        rows.append(part)
        for i, f in enumerate(feats):
            cls[f][lv] = _classify(sums["unique"][i], sums["em_expected"][i],
                                   sums["lca"][i] if "lca" in sums else None)
    by_group = pd.concat(rows, ignore_index=True)

    summary = pd.DataFrame({"feature_id": feats, "label": lab})
    for k in src:
        z = by_group[by_group[f"struc_zero_{k}"]]
        hit = z.groupby("feature_id", sort=False)["level"].agg(lambda s: ",".join(map(str, s)))
        summary[f"struc_zero_levels_{k}"] = [hit.get(f, "") for f in feats]
    if draws is not None:
        for lv in levels:
            summary[f"imputation_all_zero_frac ({group} = {lv})"] = zero_frac[lv]
        fr = np.column_stack([zero_frac[lv] for lv in levels])
        summary["inconsistent"] = ((fr > 0) & (fr < 1)).any(axis=1)
    summary["interpretation"] = [interpretation(cls[f]) for f in feats]
    return by_sample, by_group, summary, levels


def run(em_dir: str | Path, metadata: str | Path, group: str, outdir: str | Path, lca_dir: str | Path | None = None,
        features: list | None = None) -> Path:
    """Read and verify the sealed inputs, write the three tables and seal the output folder."""
    em = load_em(em_dir)
    lca = load_lca(lca_dir, em["rec"]) if lca_dir else None
    samples, groups = _samples(metadata, group, list(em["counts"].columns),
                               list(lca["counts"].columns) if lca else None)
    by_sample, by_group, summary, levels = detection_tables(em, lca, samples, groups, features, group)
    outdir = Path(outdir)
    manifest.begin(outdir, "detection")
    sources = "unique reads, EM expected counts" + (", LCA counts" if lca else "") + \
              (f", {len(em['draws'])} imputations" if em["draws"] else "")
    note = (f"# detection evidence ({sources}); level = values of metadata column '{group}' "
            f"({', '.join(map(str, levels))}); presence is count != 0, ANCOM-BC2's structural-zero rule, for every "
            "source; EM expected counts are reported as written (a fraction of a shared read is a nonzero count); "
            "no threshold is applied\n")
    for name, df in zip(OUTPUTS, (by_sample, by_group, summary)):
        with staged(outdir / name) as tmp:
            with open(tmp, "w") as fo:
                fo.write(note)
                df.to_csv(fo, sep="\t", index=False)
    rec = {"step": "detection", "group": group, "em": str(Path(em_dir).resolve()),
           "lca": str(Path(lca_dir).resolve()) if lca_dir else None, "imputations": len(em["draws"]),
           "inputs_sha256": {"em/counts.tsv": em["rec"]["artifacts"]["counts.tsv"],
                             "em/unique_counts.tsv": em["rec"]["artifacts"]["unique_counts.tsv"],
                             **({"lca/counts.tsv": lca["rec"]["artifacts"]["counts.tsv"]} if lca else {}),
                             "metadata": manifest.sha256(metadata)}}
    manifest.seal(outdir, list(OUTPUTS), rec)
    log.info("detection: %d features x %d levels -> %s", len(summary), len(levels), outdir)
    return outdir / "detection_summary.tsv"
