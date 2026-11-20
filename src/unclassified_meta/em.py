"""Reads that fit several targets equally well: EM expected counts, multiple imputation and Rubin's rules.

A read (or read pair) that a classifier reports for several taxa (Centrifuger `-k` > 1), or a fragment whose best
alignments fall in several panel features, is ambiguous. Three common treatments each bias the counts:
  - lowest common ancestor: the read leaves the rank, so the most ambiguous taxa lose reads;
  - an even split (what `centrifuger-quant` writes in its read-count column): rare taxa gain reads from common
    relatives;
  - discarding it (the panel ledger 'ambiguous'): close relatives lose reads.

Here the reads are grouped into equivalence classes (reads with the same set of candidate targets) and the
maximum-likelihood target proportions are found with the EM algorithm of a mixture of multinomials: each
ambiguous read is shared among its candidates in proportion to their estimated abundance, which is learnt mainly
from the reads that fit one target only (the same E-step as centrifuger-quant's abundance EM, applied to read
counts and without genome-length scaling, because the counts are what the differential-abundance test needs).
Independent groups of targets are fitted separately. A group whose classes are laminar (disjoint or nested) is
solved in closed form; any other group by EM accelerated with SQUAREM (Varadhan & Roland 2008,
doi:10.1111/j.1467-9469.2007.00585.x), which stops when that group's log-likelihood is certified to be within 0.5 of
its own maximum (the tangent-plane bound max_i grad_i - N of the concave likelihood on the simplex). In directions
where the likelihood is that flat, the counts are not identified by the data to better than this, and the
imputations below show their spread.

The expected counts are one point estimate. The uncertainty of the assignment is carried into the test by
multiple imputation (Rubin 1987; Little & Rubin, Statistical Analysis with Missing Data, 3rd ed., 2019): each
imputation is the last state of an independent data-augmentation chain (Tanner & Wong 1987), a Gibbs sampler that
alternates the assignment of every ambiguous read given the proportions and the proportions given the
assignments (Dirichlet with a Jeffreys prior of 0.5 per candidate), started from proportions drawn from that
prior. Where the reads inform the proportions the chain reaches their posterior within a few sweeps (one read of
its own is enough for a candidate to keep its posterior share); where the likelihood is flat over millions of
shared reads a chain barely moves. Groups whose classes are laminar are therefore drawn exactly from their
Dirichlet-tree posterior; only groups with crossing classes use a chain, started by sampling-importance-resampling
from the exact posterior of their largest-read laminar part weighted by the likelihood of the crossing classes,
so that a few crossing reads that carry real information move the start where they should. When the crossing
classes hold very many reads the weights concentrate and the start, and with it the draws, are approximate.
Each imputed table is tested and the M estimates are combined with Rubin's rules, with the Barnard & Rubin (1999)
degrees of freedom; the same idea as the inferential replicates of Salmon used by fishpond/swish (Zhu et al.
2019, doi:10.1093/nar/gkz622).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

PRIOR = 0.5            # Dirichlet prior (Jeffreys) of the imputation chains
SWEEPS = 500           # Gibbs sweeps per imputation chain


@dataclass
class Fit:
    counts: dict          # target -> expected reads (float); sums to the number of reads
    theta: dict           # target -> estimated proportion
    iterations: int
    converged: bool


class _Classes:
    """Equivalence classes as flat index arrays for vectorised E-steps."""

    def __init__(self, classes: dict):
        classes = {tuple(k): float(n) for k, n in classes.items() if n}
        self.targets = sorted({t for k in classes for t in k})
        index = {t: i for i, t in enumerate(self.targets)}
        self.T = len(self.targets)
        self.total = sum(classes.values())
        self.single = np.zeros(self.T)
        multi = [(k, n) for k, n in classes.items() if len(k) > 1]
        for k, n in classes.items():
            if len(k) == 1:
                self.single[index[k[0]]] += n
        self.n = np.array([n for _, n in multi], dtype=float)
        self.cls = np.concatenate([np.full(len(k), i) for i, (k, _) in enumerate(multi)]).astype(np.int64) \
            if multi else np.zeros(0, dtype=np.int64)
        self.tgt = np.array([index[t] for k, _ in multi for t in k], dtype=np.int64)
        self.K = len(multi)
        self.width = max((len(k) for k, _ in multi), default=0)
        # padded (classes x width) layout for vectorised multinomial draws
        self.pad = np.full((self.K, self.width), -1, dtype=np.int64)
        for i, (k, _) in enumerate(multi):
            self.pad[i, :len(k)] = [index[t] for t in k]

    def estep(self, theta: np.ndarray) -> np.ndarray:
        """Expected reads per target given proportions theta."""
        out = self.single.copy()
        if self.K:
            w = theta[self.tgt]
            s = np.bincount(self.cls, weights=w, minlength=self.K)
            even = s[self.cls] <= 0
            share = np.where(even, 0.0, w / np.where(s[self.cls] > 0, s[self.cls], 1.0))
            if even.any():                       # every candidate at zero: split evenly
                size = np.bincount(self.cls, minlength=self.K)[self.cls]
                share = np.where(even, 1.0 / size, share)
            out += np.bincount(self.tgt, weights=share * self.n[self.cls], minlength=self.T)
        return out

    def loglik(self, theta: np.ndarray) -> float:
        with np.errstate(divide="ignore"):
            ll = float(np.dot(self.single, np.log(np.where(self.single > 0, theta, 1.0))))
            if self.K:
                s = np.bincount(self.cls, weights=theta[self.tgt], minlength=self.K)
                ll += float(np.dot(self.n, np.log(s)))
        return ll


def _gap(c: "_Classes", theta: np.ndarray) -> float:
    """Certified bound on how far the log-likelihood is below its maximum. The log-likelihood is concave on the
    simplex and theta . grad = N (the reads), so by the tangent plane l* - l(theta) <= max_i grad_i - N."""
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(c.single > 0, c.single / np.where(theta > 0, theta, 1.0), 0.0)
        g = np.where((c.single > 0) & (theta <= 0), np.inf, g)
        if c.K:
            s = np.bincount(c.cls, weights=theta[c.tgt], minlength=c.K)
            w = np.where(s > 0, c.n / np.where(s > 0, s, 1.0), np.inf)
            g = g + np.bincount(c.tgt, weights=w[c.cls], minlength=c.T)
    return float(np.max(g) - c.total)


def _squarem(c: "_Classes", theta: np.ndarray, gap_tol: float, max_iter: int) -> tuple[np.ndarray, int, bool]:
    """SQUAREM (scheme S3) with step-length backtracking towards the plain EM step: a proposal is kept strictly
    positive where the current estimate is (projected to at least a thousandth of the current value, never to
    zero, so a later iteration can still give a target reads) and must not lower the likelihood below two plain
    EM steps. Stops when this group's certified likelihood gap (_gap) is at most gap_tol nats."""
    F = lambda th: c.estep(th) / c.total
    it = 0
    support = theta > 0
    if _gap(c, theta) <= gap_tol:
        return theta, 0, True
    while it < max_iter:
        it += 1
        t1 = F(theta)
        t2 = F(t1)
        r, v = t1 - theta, t2 - 2 * t1 + theta
        nv = float(np.sqrt(np.dot(v, v)))
        alpha = min(-float(np.sqrt(np.dot(r, r))) / nv, -1.0) if nv > 0 else -1.0
        ll2 = c.loglik(t2)
        new = t2
        while alpha < -1.0:
            prop = theta - 2 * alpha * r + alpha * alpha * v
            # projected onto the simplex interior: a coordinate may shrink a thousandfold per step but never reach
            # zero, so a target the extrapolation overshoots can still regain reads later
            prop = np.where(support, np.maximum(prop, np.maximum(theta * 1e-3, 1e-300)), 0.0)
            prop = F(prop / prop.sum())
            if c.loglik(prop) >= ll2 - 1e-9:
                new = prop
                break
            alpha = (alpha - 1.0) / 2.0
        theta = new
        if _gap(c, theta) <= gap_tol:
            return theta, it, True
    return theta, it, False


def dominated(classes: dict) -> set:
    """Targets whose maximum-likelihood share is exactly zero: every class that holds B also holds some A, and A
    is also in a class without B. Moving B's share to A then leaves every class containing B unchanged and raises
    the others, so the likelihood never falls. Removing such targets first makes the boundary MLE exact (plain EM
    would approach it only slowly). Targets with identical class sets are not dominated: their split is not
    identified by the data and is left to the even start."""
    members: dict = {}
    for k, n in classes.items():
        if not n:
            continue
        for t in k:
            members.setdefault(t, []).append(frozenset(k))
    out = set()
    for b, ks in members.items():
        cand = set(ks[0])
        for k in ks[1:]:
            cand &= k
            if len(cand) <= 1:
                break
        cand.discard(b)
        if any(len(members[a]) > len(ks) for a in cand):
            out.add(b)
    return out


def round_preserving(values: dict, decimals: int = 6) -> dict:
    """Round to `decimals` while keeping the rounded values' sum equal to the rounded total (largest remainder),
    so that expected counts written with a fixed number of decimals still add up to the reads they share."""
    if not values:
        return {}
    scale = 10 ** decimals
    keys = list(values)
    exact = np.array([values[k] * scale for k in keys], dtype=float)
    floor = np.floor(exact)
    short = int(round(exact.sum())) - int(floor.sum())
    order = np.argsort(-(exact - floor))
    floor[order[:max(short, 0)]] += 1
    return {k: float(f) / scale for k, f in zip(keys, floor)}


def _class_components(classes: dict) -> list[dict]:
    """Split the classes into independent groups (targets joined through shared reads)."""
    parent: dict = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for k in classes:
        for t in k[1:]:
            a, b = find(k[0]), find(t)
            if a != b:
                parent[b] = a
        find(k[0])
    groups: dict = {}
    for k, n in classes.items():
        groups.setdefault(find(k[0]), {})[k] = n
    return list(groups.values())


class _Tree:
    """A laminar group of classes (any two are disjoint or nested) as a tree. The likelihood is then
    prod over edges phi_edge ** N_subtree (N = reads of every class inside the child's subtree), so both the MLE
    (phi = N_c / sum N) and the posterior under the flat Dirichlet(prior) (Dirichlet(prior * targets_c + N_c) at
    every node, by the aggregation property of the Dirichlet) are exact."""

    def __init__(self, classes: dict):
        self.classes = classes
        multi = sorted((k for k in classes if len(k) > 1), key=len, reverse=True)
        self.single = {k[0]: float(n) for k, n in classes.items() if len(k) == 1}
        targets = sorted({t for k in classes for t in k})
        ROOT = ("#root",)
        where = {t: ROOT for t in targets}           # smallest node seen so far that holds the target
        self.parent: dict = {}
        self.ok = True
        for k in multi:
            ps = {where[t] for t in k}
            if len(ps) != 1:                         # crossing classes: not laminar
                self.ok = False
                return
            self.parent[k] = ps.pop()
            for t in k:
                where[t] = k
        self.root = ROOT
        self.children: dict = {ROOT: []}
        for k in multi:
            self.children.setdefault(k, [])
        for k, par in self.parent.items():
            self.children[par].append(k)
        for t in targets:
            self.children[where[t]].append(t)       # a target is a leaf below its smallest class
        self.ntargets, self.nreads = {}, {}
        self._sizes(ROOT)

    def _sizes(self, v):
        if isinstance(v, str):                       # a target
            self.ntargets[v], self.nreads[v] = 1, self.single.get(v, 0.0)
            return
        nt, nr = 0, float(self.classes.get(v, 0)) if v != self.root else 0.0
        for c in self.children[v]:
            self._sizes(c)
            nt, nr = nt + self.ntargets[c], nr + self.nreads[c]
        self.ntargets[v], self.nreads[v] = nt, nr

    def shares(self, v, rng=None, prior=0.5) -> np.ndarray:
        """Child shares at node v: the MLE (rng None; even per target when no child has a read) or a posterior
        draw."""
        ch = self.children[v]
        n = np.array([self.nreads[c] for c in ch])
        if rng is not None:
            return rng.dirichlet(np.array([prior * self.ntargets[c] for c in ch]) + n)
        if n.sum() > 0:
            return n / n.sum()
        size = np.array([self.ntargets[c] for c in ch], dtype=float)
        return size / size.sum()

    def flow(self, rng=None, prior=0.5, draw=False) -> dict:
        """Reads per target: expected (MLE shares) or, with rng, one exact posterior imputation."""
        out = dict(self.single)

        def down(v, arriving):
            ch = self.children[v]
            if not ch:
                return
            own = self.classes.get(v, 0) if v != self.root else 0
            if draw:                                   # whole reads stay whole (no float round trip)
                k = int(arriving) + int(own)
            else:
                k = arriving + float(own)
            phi = self.shares(v, rng, prior)
            parts = (rng.multinomial(k, phi) if draw else k * phi) if k > 0 else np.zeros(len(ch))
            for c, x in zip(ch, parts):
                if isinstance(c, str):
                    out[c] = out.get(c, 0.0) + float(x)
                else:
                    down(c, int(x) if draw else float(x))
        down(self.root, 0 if draw else 0.0)
        return out


def _fit_component(classes: dict, gap_tol: float, max_iter: int) -> tuple[dict, int, bool]:
    """Expected reads of one independent group of targets."""
    multi = [k for k in classes if len(k) > 1]
    if not multi:
        return {k[0]: float(n) for k, n in classes.items()}, 0, True
    tree = _Tree(classes)
    if tree.ok:                                  # laminar (one shared class, or nested ones): closed form
        return tree.flow(), 0, True
    c = _Classes(classes)
    dom = dominated(classes)
    zero = np.array([t in dom for t in c.targets])
    counts = c.single.copy()                      # start from the even split among the non-dominated candidates
    live = (~zero[c.tgt]).astype(float)
    size = np.bincount(c.cls, weights=live, minlength=c.K)
    counts += np.bincount(c.tgt, weights=live * (c.n / np.where(size > 0, size, 1))[c.cls], minlength=c.T)
    theta = counts / c.total
    theta, it, converged = _squarem(c, theta, gap_tol, max_iter)
    counts = theta * c.total                     # the certified proportions themselves (no further step)
    return {t: float(counts[i]) for i, t in enumerate(c.targets)}, it, converged


GAP_TOL = 0.5          # certified log-likelihood distance to the maximum at which the EM stops (LR statistic <= 1)


def expected_counts(classes: dict, gap_tol: float = GAP_TOL, max_iter: int = 200000) -> Fit:
    """Maximum-likelihood expected reads per target from equivalence classes {tuple of candidates: reads}.
    Independent groups of targets are fitted separately: in closed form when the group's classes are laminar
    (disjoint or nested; e.g. one shared class), by SQUAREM otherwise."""
    classes = {tuple(k): float(n) for k, n in classes.items() if n}
    total = sum(classes.values())
    if not total:
        return Fit({}, {}, 0, True)
    counts, its, ok = {}, 0, True
    for comp in _class_components(classes):
        c, it, conv = _fit_component(comp, gap_tol, max_iter)
        counts.update(c)
        its, ok = max(its, it), ok and conv
    return Fit(counts, {t: v / total for t, v in counts.items()}, its, ok)


def _draw_assignments(c: _Classes, theta: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Reads per target with every ambiguous read drawn to one candidate (vectorised conditional binomials)."""
    out = c.single.copy()
    if not c.K:
        return out
    w = np.where(c.pad >= 0, theta[np.clip(c.pad, 0, None)], 0.0)
    rest = np.rint(c.n).astype(np.int64)
    left = w.sum(axis=1)
    zero = left <= 0
    if zero.any():                                # all candidates at zero: uniform
        w[zero] = (c.pad[zero] >= 0).astype(float)
        left = w.sum(axis=1)
    for j in range(c.width):
        p = np.where(left > 0, w[:, j] / np.where(left > 0, left, 1.0), 0.0)
        p = np.clip(p, 0.0, 1.0)
        last = (j == c.width - 1) | (c.pad[:, j + 1] < 0 if j + 1 < c.width else True)
        x = np.where(last, rest, rng.binomial(rest, p))
        x = np.where(c.pad[:, j] >= 0, x, 0)
        np.add.at(out, np.clip(c.pad[:, j], 0, None), np.where(c.pad[:, j] >= 0, x, 0))
        rest = rest - x
        left = left - w[:, j]
    return out


def _components(c: "_Classes") -> tuple[np.ndarray, np.ndarray]:
    """(component of every target, number of ambiguous classes in each component): targets joined by sharing reads."""
    parent = np.arange(c.T)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for row in c.pad:
        idx = row[row >= 0]
        for t in idx[1:]:
            a, b = find(idx[0]), find(t)
            if a != b:
                parent[b] = a
    comp = np.array([find(t) for t in range(c.T)])
    nclass = np.bincount(comp[c.pad[:, 0]], minlength=c.T) if c.K else np.zeros(c.T, dtype=int)
    return comp, nclass


def _laminar_part(classes: dict) -> dict:
    """The largest-read sub-family of laminar classes (classes taken by decreasing reads, each kept when it is
    disjoint from or nested with every class kept so far); singletons are always kept."""
    kept: list = []
    by_target: dict = {}
    for k in sorted((k for k in classes if len(k) > 1), key=lambda k: -classes[k]):
        s = frozenset(k)
        others = {o for t in s for o in by_target.get(t, ())}
        if all(not (s & o) or s <= o or o <= s for o in others):
            kept.append(k)
            for t in s:
                by_target.setdefault(t, set()).add(s)
    return {**{k: n for k, n in classes.items() if len(k) == 1}, **{k: classes[k] for k in kept}}


def _tree_theta(tree: "_Tree", rng, prior: float) -> dict:
    """Target proportions drawn from the Dirichlet-tree posterior of a laminar family."""
    out = {}

    def down(v, mass):
        ch = tree.children[v]
        if not ch:
            return
        for c, x in zip(ch, tree.shares(v, rng, prior)):
            if isinstance(c, str):
                out[c] = out.get(c, 0.0) + mass * float(x)
            else:
                down(c, mass * float(x))
    down(tree.root, 1.0)
    return out


def _chain(classes: dict, m: int, rng: np.random.Generator, prior: float, sweeps: int) -> list[dict]:
    """Independent data-augmentation chains over crossing classes. Each starts from an exact posterior draw of
    the laminar part of the classes holding the most reads (a crossing class with few reads then barely moves the
    chain from a correct start), and the chain adds the crossing classes."""
    c = _Classes(classes)
    part = _laminar_part(classes)
    for t in c.targets:                           # targets only in crossing classes: prior leaves of the tree
        part.setdefault((t,), 0)
    lam = _Tree(part)
    crossing = [(np.array([c.targets.index(t) for t in k]), float(n)) for k, n in classes.items()
                if len(k) > 1 and k not in part]
    out = []
    for _ in range(m):
        if lam.ok:
            # sampling-importance-resampling start: candidates from the exact posterior of the laminar part,
            # weighted by the likelihood of the crossing classes it leaves out
            cands = []
            for _ in range(256):
                d = _tree_theta(lam, rng, prior)
                cands.append(np.array([d.get(t, 0.0) for t in c.targets]))
            cands = np.array(cands)
            with np.errstate(divide="ignore"):
                logw = np.zeros(len(cands))
                for idx, n in crossing:
                    logw += n * np.log(cands[:, idx].sum(axis=1))
            logw = np.where(np.isfinite(logw), logw, -np.inf)
            w = np.exp(logw - logw.max()) if np.isfinite(logw.max()) else np.ones(len(cands))
            theta = cands[rng.choice(len(cands), p=w / w.sum())] + 1e-300
            theta /= theta.sum()
        else:
            theta = rng.dirichlet(c.single + prior)
        draw = c.single
        for _ in range(max(1, sweeps)):
            draw = _draw_assignments(c, theta, rng)
            theta = rng.dirichlet(draw + prior)
        out.append({t: float(draw[i]) for i, t in enumerate(c.targets)})
    return out


def imputations(classes: dict, fit: Fit, m: int, rng: np.random.Generator, prior: float = PRIOR,
                sweeps: int = SWEEPS) -> list[dict]:
    """M imputations: integer reads per target (dict), each summing to the reads.

    Each independent group of targets is drawn on its own. A group whose classes are laminar (one shared class,
    or nested classes) is drawn exactly from its posterior: Dirichlet shares at every node of the class tree, then
    multinomial draws down the tree. A group with crossing classes is drawn by an independent data-augmentation
    chain (see the module notes for its limits)."""
    classes = {tuple(k): int(n) for k, n in classes.items() if n}
    if not classes or m <= 0:
        return [dict() for _ in range(max(m, 0))]
    out = [dict() for _ in range(m)]
    crossing = {}
    for comp in _class_components(classes):
        if all(len(k) == 1 for k in comp):
            for d in out:
                for k, n in comp.items():
                    d[k[0]] = d.get(k[0], 0) + n
            continue
        tree = _Tree(comp)
        if tree.ok:
            for d in out:
                for t, v in tree.flow(rng, prior, draw=True).items():
                    d[t] = d.get(t, 0) + v
        else:
            crossing.update(comp)
    for comp in _class_components(crossing) if crossing else []:      # each group its own chains and start
        for d, x in zip(out, _chain(comp, m, rng, prior, sweeps)):
            for t, v in x.items():
                d[t] = d.get(t, 0) + v
    return [{t: int(round(v)) for t, v in d.items() if round(v)} for d in out]


# ---------------------------------------------------------------- tail probabilities (no scipy needed)

def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction of the regularized incomplete beta function (Lentz's method)."""
    tiny, qab, qap, qam = 1e-300, a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h


def _lgamma_ratio(a: float, b: float) -> float:
    """log Gamma(a + b) - log Gamma(a), without cancellation when a is huge (asymptotic series for b = 1/2)."""
    if a > 1e5 and b == 0.5:
        return 0.5 * math.log(a) - 1 / (8 * a) + 1 / (192 * a ** 3)
    return math.lgamma(a + b) - math.lgamma(a)


def _betainc(a: float, b: float, x: float, x1: float | None = None) -> float:
    """Regularized incomplete beta I_x(a, b); x1 = 1 - x when the caller knows it more exactly than 1 - x."""
    if x1 is None:
        x1 = 1.0 - x
    if x <= 0:
        return 0.0
    if x1 <= 0:
        return 1.0
    la = math.log(x) if x < 0.5 else math.log1p(-x1)
    lb = math.log(x1) if x1 < 0.5 else math.log1p(-x)
    lbt = _lgamma_ratio(a, b) - math.lgamma(b) + a * la + b * lb
    if x < (a + 1) / (a + b + 2):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, x1) / b


def two_sided_p(z: float, df: float = float("inf")) -> float:
    """Two-sided tail probability of z under the normal (df = inf) or Student t distribution."""
    if math.isnan(z) or math.isnan(df):
        return float("nan")
    if math.isinf(z):
        return 0.0
    if math.isinf(df):
        return math.erfc(abs(z) / math.sqrt(2.0))
    if df > 1e7:
        # t tail as the normal tail plus the first two terms of its 1/df expansion (Hill 1970); accurate to far
        # below 1e-9 relative here, where the incomplete beta function loses its precision
        a = abs(z)
        phi = math.exp(-0.5 * a * a) / math.sqrt(2 * math.pi)
        g = (a ** 3 + a) / (4 * df) + (5 * a ** 7 + 16 * a ** 5 + 3 * a ** 3 - 3 * a) / (96 * df * df)
        return min(1.0, math.erfc(a / math.sqrt(2.0)) + 2 * phi * g)
    # p = I_y(df/2, 1/2) with y = df/(df+z^2); for y >= 1/2 use 1 - I_{1-y}(1/2, df/2) so that neither a small
    # p nor a p near 1 is formed by cancellation
    if z == 0:
        return 1.0
    z2 = z * z
    return _betainc(df / 2.0, 0.5, df / (df + z2), z2 / (df + z2))   # I_y(df/2, 1/2), with 1 - y given exactly


def df_from_p(stat: float, p: float, lo: float = 0.5, hi: float = 1e7) -> float:
    """Degrees of freedom of a two-sided t test with statistic `stat` and p-value `p` (inf when normal)."""
    if not (math.isfinite(stat) and math.isfinite(p)) or p <= 0 or p >= 1 or stat == 0:
        return float("nan")
    if two_sided_p(stat, hi) >= p * (1 - 1e-9):
        return float("inf")
    if two_sided_p(stat, lo) <= p:
        return lo
    a, b = math.log(lo), math.log(hi)
    for _ in range(100):                          # p decreases as df grows
        mid = 0.5 * (a + b)
        if two_sided_p(stat, math.exp(mid)) > p:
            a = mid
        else:
            b = mid
    return math.exp(0.5 * (a + b))


def rubin(estimates: np.ndarray, ses: np.ndarray, df_complete=None) -> dict:
    """Combine M fits (features x M arrays of estimates and standard errors) with Rubin's rules.

    T = W + (1 + 1/M) B, W the mean within-fit variance and B the between-fit variance of the estimates.
    Reference: t with Barnard & Rubin's (1999) degrees of freedom from the complete-data degrees of freedom
    `df_complete` (one per feature; inf = normal, the default), the formula of mice's `barnard.rubin`:
    nu = 1 / (1 / nu_m + 1 / nu_obs), with nu_m = (M - 1) / lambda^2 and
    nu_obs = nu_com (nu_com + 1) / (nu_com + 3) (1 - lambda).
    The same formula holds when the imputations agree (lambda = 0): nu = nu_obs = nu_com (nu_com + 1) / (nu_com + 3),
    smaller than nu_com (8.46 for nu_com = 10), so the pooled test is a little more conservative than one fit's.
    A feature missing (NaN) from any fit is not combined (NaN).
    """
    est = np.asarray(estimates, dtype=float)
    se = np.asarray(ses, dtype=float)
    nf, m = est.shape
    nu_com = np.full(nf, np.inf) if df_complete is None else np.asarray(df_complete, dtype=float)
    ok = np.isfinite(est).all(axis=1) & np.isfinite(se).all(axis=1) & (se > 0).all(axis=1) & ~np.isnan(nu_com)
    q = np.full(nf, np.nan)
    t_se, df, p, fmi = q.copy(), q.copy(), q.copy(), q.copy()
    for i in np.flatnonzero(ok):
        e, s = est[i], se[i]
        qbar = e.mean()
        w = float((s ** 2).mean())
        b = float(e.var(ddof=1)) if m > 1 else 0.0
        if b <= 1e-12 * w:                       # identical fits up to floating-point noise
            b = 0.0
        t = w + (1 + 1 / m) * b
        lam = (1 + 1 / m) * b / t if t > 0 else 0.0
        nc = nu_com[i]
        nu_m = (m - 1) / lam ** 2 if lam > 0 else math.inf
        if math.isinf(nc):
            nu = nu_m
        else:
            nu_obs = (nc + 1) / (nc + 3) * nc * (1 - lam)
            nu = (1 / (1 / nu_m + 1 / nu_obs) if math.isfinite(nu_m) else nu_obs) if nu_obs > 0 else 0.0
        z = qbar / math.sqrt(t) if t > 0 else (math.inf if qbar else float("nan"))
        q[i], t_se[i], df[i], fmi[i] = qbar, math.sqrt(t), nu, lam
        p[i] = two_sided_p(z, nu) if (not math.isnan(z) and nu > 1e-8) else 1.0
    return {"estimate": q, "se": t_se, "df": df, "p": p, "fmi_between": fmi}
