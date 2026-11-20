"""The faster EM and imputation (em.py): the same estimator and the same sampling distributions as before.

Reference implementations of the earlier code paths (one Dirichlet per tree node, conditional binomials per
class, one chain at a time) are kept here so that the vectorised versions are checked against them."""
import numpy as np
import pytest

from unclassified_meta import em

CROSSING = {("A",): 40, ("B",): 3, ("A", "B"): 200, ("B", "C"): 30, ("A", "C", "D"): 12, ("C", "D"): 5,
            ("D", "E"): 2, ("E",): 1, ("A", "E"): 7}


def _mle(classes, iters=200_000):
    """Plain EM to far below any tolerance used here (the reference maximum)."""
    c = em._Classes(classes)
    theta = np.full(c.T, 1.0 / c.T)
    for _ in range(iters):
        new = c.estep(theta) / c.total
        if np.abs(new - theta).max() < 1e-15:
            break
        theta = new
    return c, theta


# ------------------------------------------------------------------ EM: certificate and estimates

def test_the_dual_bound_never_falls_below_the_true_gap():
    c, best = _mle(CROSSING)
    lstar = c.loglik(best)
    rng = np.random.default_rng(3)
    for _ in range(300):
        theta = rng.dirichlet(np.full(c.T, 0.3)) * 0.2 + 0.8 * best
        s = c.sums(theta)
        g = em._grad(c, theta, s)
        true = lstar - c.loglik(theta, s)
        assert em._dual_gap(c, theta, s, g, 1e-3) >= true - 1e-9
        assert em._gap(c, theta, s) >= true - 1e-9


def test_the_dual_bound_is_never_looser_than_the_tangent_bound():
    c, best = _mle(CROSSING)
    theta = 0.9 * best + 0.1 / c.T
    s = c.sums(theta)
    g = em._grad(c, theta, s)
    assert em._dual_gap(c, theta, s, g, 1e-3) <= em._gap(c, theta, s) + 1e-12


def test_a_tiny_share_that_holds_the_tangent_bound_open_is_certified_by_the_dual_bound():
    # B has no reads of its own and shares a few with a huge A: its gradient exceeds N for a very long time
    classes = {("A",): 10_000_000, ("A", "B"): 40, ("B", "C"): 3, ("C",): 2}
    fit = em.expected_counts(classes, dual_tol=1e-3)
    c, best = _mle(classes, iters=2_000_000)
    theta = np.array([fit.theta[t] for t in c.targets])
    assert fit.converged and c.loglik(best) - c.loglik(theta) <= 1e-3 + 1e-9


@pytest.mark.parametrize("classes", [CROSSING,
                                     {("A", "B"): 1_000, ("C", "D"): 500, ("A", "C"): 20, ("B", "D"): 1,
                                      ("A",): 3},
                                     {("A",): 0.5e6, ("A", "B", "C"): 1e5, ("B", "D"): 10, ("C", "D"): 15,
                                      ("D",): 1}])
def test_expected_counts_reach_the_maximum_likelihood_within_the_certified_tolerance(classes):
    fit = em.expected_counts(classes, dual_tol=1e-4)
    c, best = _mle(classes)
    theta = np.array([fit.theta[t] for t in c.targets])
    assert fit.converged
    assert c.loglik(best) - c.loglik(theta) <= 1e-4 + 1e-9
    assert sum(fit.counts.values()) == pytest.approx(sum(classes.values()), rel=1e-12)


def test_without_the_dual_bound_the_iterates_are_those_of_plain_squarem():
    """dual_tol None reproduces the earlier stop rule (tangent-plane gap only)."""
    c = em._Classes(CROSSING)
    theta0 = np.full(c.T, 1.0 / c.T)
    a, ia, _ = em._squarem(c, theta0, 1e-6, 10_000, None)
    s = c.sums(a)
    assert em._gap(c, a, s) <= 1e-6 and ia > 0


def test_a_target_with_zero_maximum_likelihood_share_gets_no_reads():
    # B is dominated by A (every class with B also has A, and A has reads of its own)
    out = em.expected_counts({("A",): 50, ("A", "B"): 20, ("A", "B", "C"): 4, ("C",): 9})
    assert out.counts["B"] == pytest.approx(0, abs=1e-9)


def test_a_dominant_taxon_takes_almost_all_shared_reads():
    out = em.expected_counts({("A",): 1_000_000, ("B",): 1, ("A", "B"): 10_000, ("B", "C"): 1, ("C",): 1})
    assert out.counts["A"] == pytest.approx(1_010_000, rel=1e-4) and out.counts["B"] < 15


def test_reads_with_the_same_candidates_form_one_class():
    """A class is the set of candidates: reads with the same set add up (Centrifuger reports every accepted
    candidate of a read with the same weight, so the set carries all of a read's information here)."""
    a = em.expected_counts({("A", "B"): 10, ("A",): 5, ("B", "C"): 3})
    b = em.expected_counts({("B", "A"): 4, ("A", "B"): 6, ("A",): 5, ("B", "C"): 3})
    assert all(b.counts[t] == pytest.approx(v, abs=1e-6) for t, v in a.counts.items())


# ------------------------------------------------------------------ imputation: class-level draws

def _old_draw(c, theta, rng):
    """The earlier per-chain draw: conditional binomials over the candidate columns for every class."""
    out = c.single.copy()
    w = np.where(c.pad >= 0, theta[np.clip(c.pad, 0, None)], 0.0)
    rest = np.rint(c.n).astype(np.int64)
    left = w.sum(axis=1)
    for j in range(c.width):
        p = np.clip(np.where(left > 0, w[:, j] / np.where(left > 0, left, 1.0), 0.0), 0.0, 1.0)
        last = (j == c.width - 1) | (c.pad[:, j + 1] < 0 if j + 1 < c.width else True)
        x = np.where(last, rest, rng.binomial(rest, p))
        x = np.where(c.pad[:, j] >= 0, x, 0)
        np.add.at(out, np.clip(c.pad[:, j], 0, None), x)
        rest = rest - x
        left = left - w[:, j]
    return out


@pytest.mark.parametrize("reads", [1, 7, em.EXPAND, em.EXPAND + 1, 500])
def test_class_level_draws_are_multinomial_with_the_class_probabilities(reads):
    classes = {("A", "B", "C"): reads, ("C", "D"): reads, ("A",): 2}
    c = em._Classes(classes)
    theta = np.array([0.5, 0.3, 0.2, 0.0])                 # D at zero: never drawn in ("C", "D")
    rng = np.random.default_rng(11)
    draws = em._draw_assignments(c, np.tile(theta, (4000, 1)), rng)
    assert (draws.sum(axis=1) == c.total).all()               # every read assigned once
    assert (draws[:, 3] == 0).all()
    a_share = (draws[:, 0] - 2) / reads                       # A's share of the first class: 0.5
    assert a_share.mean() == pytest.approx(0.5, abs=4 * np.sqrt(0.25 / reads / 4000))
    assert a_share.var() * reads == pytest.approx(0.25, rel=0.1)
    old = np.array([_old_draw(c, theta, rng) for _ in range(4000)])
    assert np.abs(old.mean(axis=0) - draws.mean(axis=0)).max() < 5 * np.sqrt(reads / 4000) + 1e-9


def test_a_class_whose_candidates_all_weigh_zero_is_split_uniformly():
    c = em._Classes({("A", "B"): 1000, ("C",): 1})
    draws = em._draw_assignments(c, np.array([[0.0, 0.0, 1.0]]), np.random.default_rng(1))
    assert draws[0, 0] + draws[0, 1] == 1000 and 400 < draws[0, 0] < 600


# ------------------------------------------------------------------ imputation: Dirichlet-tree starts

def _old_tree_theta(tree, rng, prior):
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


def test_batched_tree_draws_follow_the_dirichlet_tree_posterior():
    classes = {("A", "B", "C", "D"): 30, ("A", "B"): 12, ("C",): 4, ("D",): 1, ("E",): 3, ("A",): 2}
    tree = em._Tree(classes)
    assert tree.ok
    targets = sorted({t for k in classes for t in k})
    col = {t: i for i, t in enumerate(targets)}
    rng = np.random.default_rng(5)
    new = em._TreeDraw(tree, col).thetas(rng, em.PRIOR, 20000, len(targets))
    old = np.array([[d[t] for t in targets] for d in (_old_tree_theta(tree, rng, em.PRIOR) for _ in range(20000))])
    assert np.allclose(new.sum(axis=1), 1.0)
    assert np.abs(new.mean(axis=0) - old.mean(axis=0)).max() < 0.01
    assert np.abs(new.std(axis=0) - old.std(axis=0)).max() < 0.01


# ------------------------------------------------------------------ imputation: whole pipeline

def test_imputations_keep_every_read_and_are_reproducible_with_a_seed():
    fit = em.expected_counts(CROSSING)
    a = em.imputations(CROSSING, fit, 5, np.random.default_rng([11, 1, 0]))
    b = em.imputations(CROSSING, fit, 5, np.random.default_rng([11, 1, 0]))
    assert a == b
    assert all(sum(d.values()) == sum(CROSSING.values()) for d in a)
    assert all(isinstance(v, int) for d in a for v in d.values())


def test_independent_groups_drawn_together_keep_their_own_distribution():
    """Chains of several crossing groups run side by side (one gamma vector over all groups); each group's
    draws must match the same group imputed alone."""
    g1 = {("A", "B"): 300, ("B", "C"): 20, ("A", "C"): 50, ("A",): 10}
    g2 = {("X", "Y"): 5000, ("Y", "Z"): 3, ("X", "Z"): 900, ("Z",): 40}
    both = {**g1, **g2}
    alone = em.imputations(g1, em.expected_counts(g1), 300, np.random.default_rng(1), sweeps=50)
    joint = em.imputations(both, em.expected_counts(both), 300, np.random.default_rng(2), sweeps=50)
    for t in "ABC":
        x, y = np.array([d.get(t, 0) for d in alone]), np.array([d.get(t, 0) for d in joint])
        assert x.mean() == pytest.approx(y.mean(), abs=5 * np.sqrt((x.var() + y.var()) / 300) + 1)
    assert all(sum(d.get(t, 0) for t in "XYZ") == sum(g2.values()) for d in joint)
