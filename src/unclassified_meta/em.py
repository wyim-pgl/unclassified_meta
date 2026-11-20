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
doi:10.1111/j.1467-9469.2007.00585.x), which stops when that group's log-likelihood is certified to be within
DUAL_TOL (1e-4) of its own maximum, by the smaller of the tangent-plane bound max_i grad_i - N of the concave
likelihood on the simplex and a Lagrangian-dual bound that stays tight where a candidate with a tiny share holds
the tangent bound open (_dual_gap). In directions where the likelihood is flat the counts are not identified by
the data, and the imputations below show their spread.

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
The chains of all crossing groups and all M imputations run side by side as arrays, and the reads of a class are
drawn together as one multinomial; these are the same distributions as one chain per group with conditional
binomials per class, with the random numbers used in a different order.
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
        self.ncls = self.n[self.cls]                  # reads of the class of every (class, candidate) entry
        self.has_single = self.single > 0

    def sums(self, theta: np.ndarray) -> np.ndarray:
        """Proportion of every ambiguous class (sum of its candidates' theta). The E-step, the log-likelihood and
        the gap at one theta all start from these, so callers compute them once and pass them on."""
        return np.bincount(self.cls, weights=theta[self.tgt], minlength=self.K)

    def estep(self, theta: np.ndarray, s: np.ndarray | None = None) -> np.ndarray:
        """Expected reads per target given proportions theta (s: self.sums(theta), when already known)."""
        out = self.single.copy()
        if self.K:
            w = theta[self.tgt]
            if s is None:
                s = np.bincount(self.cls, weights=w, minlength=self.K)
            sc = s[self.cls]
            if sc.min() > 0:                         # the usual case; the same arithmetic as the general branch
                share = w / sc
            else:
                even = sc <= 0
                share = np.where(even, 0.0, w / np.where(sc > 0, sc, 1.0))
                size = np.bincount(self.cls, minlength=self.K)[self.cls]     # every candidate at zero: even split
                share = np.where(even, 1.0 / size, share)
            out += np.bincount(self.tgt, weights=share * self.ncls, minlength=self.T)
        return out

    def loglik(self, theta: np.ndarray, s: np.ndarray | None = None) -> float:
        with np.errstate(divide="ignore"):
            ll = float(np.dot(self.single, np.log(np.where(self.has_single, theta, 1.0))))
            if self.K:
                if s is None:
                    s = self.sums(theta)
                ll += float(np.dot(self.n, np.log(s)))
        return ll


def _grad(c: "_Classes", theta: np.ndarray, s: np.ndarray | None = None) -> np.ndarray:
    """Gradient of the log-likelihood in theta (s: c.sums(theta), when already known)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(c.has_single, c.single / np.where(theta > 0, theta, 1.0), 0.0)
        g = np.where(c.has_single & (theta <= 0), np.inf, g)
        if c.K:
            if s is None:
                s = c.sums(theta)
            w = np.where(s > 0, c.n / np.where(s > 0, s, 1.0), np.inf)
            g = g + np.bincount(c.tgt, weights=w[c.cls], minlength=c.T)
    return g


def _gap(c: "_Classes", theta: np.ndarray, s: np.ndarray | None = None) -> float:
    """Certified bound on how far the log-likelihood is below its maximum. The log-likelihood is concave on the
    simplex and theta . grad = N (the reads), so by the tangent plane l* - l(theta) <= max_i grad_i - N.
    (s: c.sums(theta), when already known.)"""
    return float(np.max(_grad(c, theta, s)) - c.total)


def _dual_gap(c: "_Classes", theta: np.ndarray, s: np.ndarray, g: np.ndarray, gap_tol: float) -> float:
    """A certified bound on l* - l(theta) at most the tangent-plane bound max_i g_i - N, often far smaller.

    Every term of the log-likelihood is n_j log(y_j), with y_j = a_j . theta the share of a class (or of one
    target, for its own reads). From log(x) <= x - 1, for any weights w_j > 0 and every theta' on the simplex,
        l(theta') <= sum_j n_j log(n_j / w_j) - N + max_i (A^T w)_i,
    the Lagrangian dual of the mixture likelihood. With w_j = n_j / y_j(theta) the right side is
    l(theta) + max_i g_i - N, the tangent-plane bound. A target with a tiny share whose gradient still exceeds
    N holds that bound open long after the likelihood has stopped changing, because EM grows its share only by
    the factor g_i / N per iteration and N counts every read of the group. Shrinking the weights of that
    target's terms by a factor (1 - eps_i) lowers its (A^T w)_i to a level lam at a cost of only
    -sum n_j log(1 - eps_i) over its own terms; the bound is evaluated exactly for the shrunken weights, so it
    stays a certificate. Returns the smaller of the two bounds."""
    tot = c.total
    plain = float(np.max(g)) - tot
    if plain <= gap_tol or not np.isfinite(plain):
        return plain
    best = plain
    with np.errstate(divide="ignore", invalid="ignore"):
        own = np.where(c.has_single, c.single / np.where(theta > 0, theta, 1.0), 0.0)   # weight of own reads
        w = c.n / s
    for frac in (0.1, 0.3, 0.6):
        lam = tot + frac * gap_tol
        eps = np.where(g > lam, (g - lam) / g, 0.0)
        if not eps.any() or eps.max() >= 1:
            continue
        log_keep = np.log1p(-eps)
        cls_keep = np.bincount(c.cls, weights=log_keep[c.tgt], minlength=c.K) if c.K else np.zeros(0)
        cost = -float(np.dot(c.n, cls_keep)) - float(np.dot(c.single, log_keep))
        g2 = own * np.exp(log_keep)
        if c.K:
            g2 = g2 + np.bincount(c.tgt, weights=(w * np.exp(cls_keep))[c.cls], minlength=c.T)
        best = min(best, cost + float(np.max(g2)) - tot)
    return best


def _squarem(c: "_Classes", theta: np.ndarray, gap_tol: float, max_iter: int,
             dual_tol: float | None = None) -> tuple[np.ndarray, int, bool]:
    """SQUAREM (scheme S3) with step-length backtracking towards the plain EM step: a proposal is kept strictly
    positive where the current estimate is (projected to at least a thousandth of the current value, never to
    zero, so a later iteration can still give a target reads) and must not lower the likelihood below two plain
    EM steps. Stops when this group's log-likelihood is certified to be within gap_tol nats of its maximum by
    the tangent-plane bound (_gap); with dual_tol given, when it is certified to be within min(gap_tol,
    dual_tol) by the smaller of the tangent-plane bound and the dual bound of _dual_gap.

    The class proportions (c.sums) of a point are computed once and shared by its log-likelihood, its gap and
    the next EM step from it; with dual_tol None the iterates and the stop are those of the original code."""
    total = c.total
    it = 0
    support = theta > 0

    tol = gap_tol if dual_tol is None else min(gap_tol, dual_tol)

    def done(th, sums):
        g = _grad(c, th, sums)
        if float(np.max(g)) - total <= tol:
            return True
        return dual_tol is not None and _dual_gap(c, th, sums, g, tol) <= tol
    s = c.sums(theta)
    if done(theta, s):
        return theta, 0, True
    while it < max_iter:
        it += 1
        t1 = c.estep(theta, s) / total
        t2 = c.estep(t1) / total
        r, v = t1 - theta, t2 - 2 * t1 + theta
        nv = float(np.sqrt(np.dot(v, v)))
        alpha = min(-float(np.sqrt(np.dot(r, r))) / nv, -1.0) if nv > 0 else -1.0
        s2 = c.sums(t2)
        ll2 = c.loglik(t2, s2)
        new, s_new = t2, s2
        while alpha < -1.0:
            prop = theta - 2 * alpha * r + alpha * alpha * v
            # projected onto the simplex interior: a coordinate may shrink a thousandfold per step but never reach
            # zero, so a target the extrapolation overshoots can still regain reads later
            prop = np.where(support, np.maximum(prop, np.maximum(theta * 1e-3, 1e-300)), 0.0)
            prop = c.estep(prop / prop.sum()) / total
            s_prop = c.sums(prop)
            if c.loglik(prop, s_prop) >= ll2 - 1e-9:
                new, s_new = prop, s_prop
                break
            alpha = (alpha - 1.0) / 2.0
        theta, s = new, s_new
        if done(theta, s):
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


def _fit_component(classes: dict, gap_tol: float, max_iter: int,
                   dual_tol: float | None = None) -> tuple[dict, int, bool]:
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
    theta, it, converged = _squarem(c, theta, gap_tol, max_iter, dual_tol)
    counts = theta * c.total                     # the certified proportions themselves (no further step)
    return {t: float(counts[i]) for i, t in enumerate(c.targets)}, it, converged


GAP_TOL = 0.5          # certified log-likelihood distance to the maximum at which the EM stops (LR statistic <= 1)
DUAL_TOL = 1e-4        # the same, certified by the dual bound (_dual_gap): far closer to the maximum than GAP_TOL
                       # and reached in a few hundred iterations where the tangent bound needed up to 1.5e5


def expected_counts(classes: dict, gap_tol: float = GAP_TOL, max_iter: int = 200000,
                    dual_tol: float | None = DUAL_TOL) -> Fit:
    """Maximum-likelihood expected reads per target from equivalence classes {tuple of candidates: reads}.
    Independent groups of targets are fitted separately: in closed form when the group's classes are laminar
    (disjoint or nested; e.g. one shared class), by SQUAREM otherwise."""
    classes = {tuple(k): float(n) for k, n in classes.items() if n}
    total = sum(classes.values())
    if not total:
        return Fit({}, {}, 0, True)
    counts, its, ok = {}, 0, True
    for comp in _class_components(classes):
        c, it, conv = _fit_component(comp, gap_tol, max_iter, dual_tol)
        counts.update(c)
        its, ok = max(its, it), ok and conv
    return Fit(counts, {t: v / total for t, v in counts.items()}, its, ok)


EXPAND = 16            # classes with at most this many reads are drawn read by read (categorical), larger ones
                       # by conditional binomials; both are exact multinomial draws


def _draw_layout(c: _Classes) -> tuple:
    """Index arrays for _draw_assignments, built once per set of classes."""
    lay = getattr(c, "_layout", None)
    if lay is None:
        reads = np.rint(c.n).astype(np.int64)
        small = reads <= EXPAND
        read_cls = np.repeat(np.flatnonzero(small), reads[small])            # one entry per read of a small class
        big = np.flatnonzero(~small)
        nvalid = (c.pad >= 0).sum(axis=1)
        lay = c._layout = (reads, read_cls, big, nvalid)
    return lay


def _draw_assignments(c: _Classes, theta: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Reads per target with every ambiguous read drawn to one candidate, for M chains at once.

    theta: (M, T) non-negative weights (only their ratios within a class matter, so unnormalised gamma draws
    serve). All reads of a class share one categorical distribution over its candidates (uniform when every
    candidate weighs zero), so the class's reads are one multinomial draw: read by read (one uniform each,
    against the class's cumulative weights) for classes with up to EXPAND reads, by conditional binomials over
    the candidate columns for larger ones. Returns (M, T) reads."""
    theta = np.atleast_2d(theta)
    m = theta.shape[0]
    out = np.broadcast_to(c.single, (m, c.T)).copy()
    if not c.K:
        return out
    reads, read_cls, big, nvalid = _draw_layout(c)
    valid = c.pad >= 0
    w = np.where(valid, theta[:, np.clip(c.pad, 0, None)], 0.0)            # (M, K, width)
    left = w.sum(axis=2)
    zero = left <= 0
    if zero.any():                                # all candidates at zero: uniform
        w[zero] = np.broadcast_to(valid, w.shape)[zero].astype(float)
        left = w.sum(axis=2)
    flat = np.arange(m)[:, None] * c.T            # row offsets into the (M*T) output
    acc = np.zeros(m * c.T)
    if len(read_cls):
        cum = np.cumsum(w[:, read_cls, :], axis=2)                          # (M, R, width)
        u = rng.random((m, len(read_cls))) * cum[:, :, -1]
        # candidate j is taken when cum[j-1] <= u < cum[j]: never one of weight zero
        pick = np.minimum((cum <= u[:, :, None]).sum(axis=2), nvalid[read_cls] - 1)
        acc += np.bincount((flat + c.pad[read_cls, pick]).ravel(), minlength=m * c.T)
    if len(big):
        pad, wb, lb = c.pad[big], w[:, big, :], left[:, big]
        rest = np.broadcast_to(reads[big], (m, len(big))).copy()
        for j in range(c.width):
            col = pad[:, j]
            ok = col >= 0
            last = ~(pad[:, j + 1] >= 0) if j + 1 < c.width else np.ones(len(big), dtype=bool)
            p = np.clip(np.where(lb > 0, wb[:, :, j] / np.where(lb > 0, lb, 1.0), 0.0), 0.0, 1.0)
            x = np.where(last, rest, rng.binomial(rest, p))
            x = np.where(ok, x, 0)
            acc += np.bincount((flat + np.clip(col, 0, None)).ravel(), weights=x.ravel(), minlength=m * c.T)
            rest = rest - x
            lb = lb - wb[:, :, j]
    return out + acc.reshape(m, c.T)


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


class _TreeDraw:
    """A laminar _Tree as flat edge arrays, to draw many target-proportion vectors from its Dirichlet-tree
    posterior at once: the shares below every node are Dirichlet(prior * targets + reads of each child), drawn
    as independent gamma variates normalised within the node; a target's proportion is the product of the
    shares on its path from the root."""

    def __init__(self, tree: "_Tree", column: dict):
        node_id = {tree.root: 0}
        parent, child, alpha_t, reads, starts, levels = [], [], [], [], [], []
        frontier, leaf_edge, leaf_col = [tree.root], [], []
        while frontier:                               # breadth first: one level of edges at a time
            level, nxt = [], []
            for v in frontier:
                ch = tree.children[v]
                if not ch:
                    continue
                starts.append(len(child))
                for c in ch:
                    node_id.setdefault(c, len(node_id))
                    level.append(len(child))
                    parent.append(node_id[v])
                    child.append(node_id[c])
                    alpha_t.append(tree.ntargets[c])
                    reads.append(tree.nreads[c])
                    if isinstance(c, str):
                        leaf_edge.append(len(child) - 1)
                        leaf_col.append(column[c])
                    else:
                        nxt.append(c)
            if level:
                levels.append(np.array(level, dtype=np.int64))
            frontier = nxt
        self.E, self.nodes = len(child), len(node_id)
        self.parent, self.child = np.array(parent, dtype=np.int64), np.array(child, dtype=np.int64)
        self.ntargets, self.reads = np.array(alpha_t, dtype=float), np.array(reads, dtype=float)
        self.starts = np.array(starts, dtype=np.int64)          # edges are contiguous per parent node
        self.group = np.repeat(np.arange(len(starts)), np.diff(np.r_[self.starts, self.E]))
        self.levels = levels
        self.leaf_node = self.child[np.array(leaf_edge, dtype=np.int64)]
        self.leaf_col = np.array(leaf_col, dtype=np.int64)

    def thetas(self, rng: np.random.Generator, prior: float, size: int, ncol: int) -> np.ndarray:
        """(size, ncol) proportions; columns of targets outside the tree stay zero."""
        out = np.zeros((size, ncol))
        if not self.E:
            return out
        g = rng.standard_gamma(prior * self.ntargets + self.reads, size=(size, self.E))
        tot = np.add.reduceat(g, self.starts, axis=1)[:, self.group]
        phi = g / np.where(tot > 0, tot, 1.0)
        mass = np.zeros((size, self.nodes))
        mass[:, 0] = 1.0
        for lv in self.levels:
            mass[:, self.child[lv]] = mass[:, self.parent[lv]] * phi[:, lv]
        out[:, self.leaf_col] = mass[:, self.leaf_node]
        return out


SIR_CANDIDATES = 256   # laminar-posterior candidates per chain start (sampling-importance-resampling)


def _sir_starts(classes: dict, targets: list, column: dict, m: int, rng: np.random.Generator,
                prior: float) -> np.ndarray:
    """(m, len(targets)) chain starts for one group with crossing classes (zero outside the group): each is
    resampled from SIR_CANDIDATES exact posterior draws of the laminar part of the group's classes holding the
    most reads, weighted by the likelihood of the crossing classes that part leaves out (so that a few
    crossing reads that carry real information move the start where they should)."""
    part = _laminar_part(classes)
    own = sorted({t for k in classes for t in k})
    for t in own:                                 # targets only in crossing classes: prior leaves of the tree
        part.setdefault((t,), 0)
    lam = _Tree(part)
    cols = np.array([column[t] for t in own], dtype=np.int64)
    out = np.zeros((m, len(targets)))
    if not lam.ok:                                # cannot happen for _laminar_part; kept as a safe fallback
        single = np.array([float(classes.get((t,), 0)) for t in own])
        out[:, cols] = rng.dirichlet(single + prior, size=m)
        return out
    local = {t: i for i, t in enumerate(own)}
    draw = _TreeDraw(lam, local)
    cross = [(k, float(n)) for k, n in classes.items() if len(k) > 1 and k not in part]
    ctgt = np.array([local[t] for k, _ in cross for t in k], dtype=np.int64)
    cstart = np.cumsum([0] + [len(k) for k, _ in cross[:-1]]).astype(np.int64)
    cn = np.array([n for _, n in cross])
    # bound the (candidates x crossing entries) block that the weights need at a time
    step = max(1, min(SIR_CANDIDATES, 4_000_000 // max(1, len(ctgt), draw.E)))
    for i in range(m):
        cands = np.empty((SIR_CANDIDATES, len(own)))
        logw = np.zeros(SIR_CANDIDATES)
        for a in range(0, SIR_CANDIDATES, step):
            b = min(SIR_CANDIDATES, a + step)
            cands[a:b] = draw.thetas(rng, prior, b - a, len(own))
            if len(cross):
                with np.errstate(divide="ignore"):
                    logw[a:b] = np.log(np.add.reduceat(cands[a:b][:, ctgt], cstart, axis=1)) @ cn
        logw = np.where(np.isfinite(logw), logw, -np.inf)
        w = np.exp(logw - logw.max()) if np.isfinite(logw.max()) else np.ones(SIR_CANDIDATES)
        theta = cands[rng.choice(SIR_CANDIDATES, p=w / w.sum())] + 1e-300
        out[i, cols] = theta / theta.sum()
    return out


def _chains(classes: dict, m: int, rng: np.random.Generator, prior: float, sweeps: int) -> list[dict]:
    """M independent data-augmentation chains over the groups with crossing classes, all groups and chains
    advanced together. Each chain of each group starts from its own SIR draw (_sir_starts) and alternates the
    assignment of every ambiguous read given the proportions and the proportions given the assignments
    (Dirichlet(reads + prior) within the group, drawn as gamma variates: only the ratios within a group's classes
    enter the next assignment, so one gamma vector over all groups is the same as one Dirichlet per group)."""
    c = _Classes(classes)
    column = {t: i for i, t in enumerate(c.targets)}
    theta = np.zeros((m, c.T))
    for comp in _class_components(classes):
        theta += _sir_starts(comp, c.targets, column, m, rng, prior)
    draw = np.broadcast_to(c.single, (m, c.T))
    for s in range(max(1, sweeps)):
        draw = _draw_assignments(c, theta, rng)
        if s + 1 < max(1, sweeps):                # the proportions after the last assignment are not used
            theta = rng.standard_gamma(draw + prior)
    return [{t: float(draw[i, j]) for j, t in enumerate(c.targets)} for i in range(m)]


def imputations(classes: dict, fit: Fit, m: int, rng: np.random.Generator, prior: float = PRIOR,
                sweeps: int = SWEEPS) -> list[dict]:
    """M imputations: integer reads per target (dict), each summing to the reads.

    Each independent group of targets is drawn on its own. A group whose classes are laminar (one shared class,
    or nested classes) is drawn exactly from its posterior: Dirichlet shares at every node of the class tree, then
    multinomial draws down the tree. Groups with crossing classes are drawn by independent data-augmentation
    chains (see the module notes for their limits), run side by side."""
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
    if crossing:
        for d, x in zip(out, _chains(crossing, m, rng, prior, sweeps)):
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
    Reference: t with Barnard & Rubin's (1999) degrees of freedom, which reduce to the complete-data degrees of
    freedom `df_complete` (one per feature; inf = normal, the default) when the imputations agree.
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
        if lam == 0:
            nu = nc                              # the imputations agree: the complete-data test itself
        elif math.isinf(nc):
            nu = nu_m
        else:
            nu_obs = (nc + 1) / (nc + 3) * nc * (1 - lam)
            nu = 1 / (1 / nu_m + 1 / nu_obs) if nu_obs > 0 else 0.0
        z = qbar / math.sqrt(t) if t > 0 else (math.inf if qbar else float("nan"))
        q[i], t_se[i], df[i], fmi[i] = qbar, math.sqrt(t), nu, lam
        p[i] = two_sided_p(z, nu) if (not math.isnan(z) and nu > 1e-8) else 1.0
    return {"estimate": q, "se": t_se, "df": df, "p": p, "fmi_between": fmi}
