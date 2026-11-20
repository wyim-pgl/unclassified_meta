"""EM read assignment and multiple imputation (em.py)."""
import numpy as np
import pytest

from unclassified_meta import em


def test_unique_reads_are_counted_as_they_are():
    out = em.expected_counts({("A",): 5, ("B",): 3})
    assert out.counts == {"A": 5.0, "B": 3.0} and out.converged


def test_ambiguous_reads_follow_the_abundance_of_the_unique_reads():
    # 90 unique A, 10 unique B, 100 shared: the MLE puts 90 of the shared reads on A
    out = em.expected_counts({("A",): 90, ("B",): 10, ("A", "B"): 100})
    assert out.counts["A"] == pytest.approx(180, abs=1e-3)
    assert out.counts["B"] == pytest.approx(20, abs=1e-3)


def test_expected_counts_keep_every_read():
    classes = {("A",): 7, ("A", "B"): 11, ("B", "C"): 13, ("A", "B", "C"): 17, ("C",): 2}
    out = em.expected_counts(classes)
    assert sum(out.counts.values()) == pytest.approx(sum(classes.values()), abs=1e-9)


def test_a_target_never_seen_alone_gets_its_share_from_the_others():
    out = em.expected_counts({("A", "B"): 10})
    assert out.counts["A"] == pytest.approx(5) and out.counts["B"] == pytest.approx(5)


def test_imputations_are_integers_that_keep_every_read_and_vary_only_for_ambiguous_reads():
    classes = {("A",): 50, ("B",): 50, ("A", "B"): 40}
    fit = em.expected_counts(classes)
    draws = em.imputations(classes, fit, 200, np.random.default_rng(1))
    assert len(draws) == 200
    for d in draws:
        assert all(isinstance(v, int) for v in d.values()) and sum(d.values()) == 140
        assert 50 <= d.get("A", 0) <= 90
    assert np.mean([d.get("A", 0) for d in draws]) == pytest.approx(fit.counts["A"], abs=2.5)


def test_imputation_is_reproducible_with_a_seed():
    classes = {("A", "B"): 25, ("B", "C"): 25, ("C",): 3}
    fit = em.expected_counts(classes)
    a = em.imputations(classes, fit, 3, np.random.default_rng(7))
    b = em.imputations(classes, fit, 3, np.random.default_rng(7))
    assert a == b


def test_em_reaches_the_mle_of_a_small_class_next_to_a_huge_one():
    # one A-only read and 1,000 A/B reads: the MLE gives B nothing, however many other reads there are
    out = em.expected_counts({("A",): 1, ("A", "B"): 1000, ("C",): 600_000_000})
    assert out.converged and out.counts["B"] < 0.01 and out.counts["A"] > 1000.99


def test_em_never_clips_a_target_that_the_mle_keeps():
    out = em.expected_counts({("A", "C"): 13, ("A", "B"): 10, ("A",): 1, ("B", "C"): 1})
    assert out.counts["A"] == pytest.approx(22.917, abs=0.05)               # within the certified tolerance
    assert out.counts["C"] == pytest.approx(2.083, abs=0.05) and out.counts["B"] < 0.05


@pytest.mark.parametrize("shared", [72_000_000, 200_000_000])
def test_em_drives_a_flat_direction_to_its_boundary_mle(shared):
    out = em.expected_counts({("A",): 1, ("A", "B"): shared, ("C",): 600_000_000})
    assert (not out.converged) or out.counts["B"] < 1.0                # never a false "converged"


def test_one_read_of_its_own_keeps_the_posterior_share_of_the_other_candidate():
    classes = {("A",): 1, ("A", "B"): 100}
    fit = em.expected_counts(classes)
    draws = [d.get("B", 0) for d in em.imputations(classes, fit, 300, np.random.default_rng(5))]
    assert np.mean(draws) == pytest.approx(25, abs=6)                  # Beta(0.5, 1.5) posterior share of B x 100


def test_overlapping_unidentified_ratios_match_exact_enumeration():
    classes = {("A", "B"): 30, ("C", "D"): 30, ("A", "C"): 30, ("B", "D"): 30}
    fit = em.expected_counts(classes)
    draws = np.array([d.get("A", 0) for d in em.imputations(classes, fit, 1200, np.random.default_rng(6))])
    assert np.mean((draws >= 12) & (draws <= 48)) == pytest.approx(0.408, abs=0.06)


def test_unidentified_proportions_spread_like_the_posterior_for_millions_of_reads():
    classes = {("A", "B"): 1_000_000}
    fit = em.expected_counts(classes)
    draws = [d.get("A", 0) for d in em.imputations(classes, fit, 60, np.random.default_rng(4))]
    assert np.std(draws) > 200_000                       # Beta(0.5, 0.5) x 1e6 has sd 353,553


def test_unidentified_proportions_are_not_fixed_at_the_even_split():
    # every read of A and B is shared: the data cannot tell their ratio, so the imputations must spread far
    # beyond the binomial spread of a fixed 50:50 split (sd 5 for 100 reads)
    classes = {("A", "B"): 100}
    fit = em.expected_counts(classes)
    draws = [d.get("A", 0) for d in em.imputations(classes, fit, 60, np.random.default_rng(3))]
    assert np.std(draws) > 10


def test_identical_imputations_use_the_barnard_rubin_observed_degrees_of_freedom():
    # nu_obs = nu (nu + 1) / (nu + 3) (1 - lambda) with lambda = 0: 4 * 5 / 7, as in mice's barnard.rubin
    r = em.rubin(np.array([[0.25, 0.25, 0.25]]), np.array([[0.1, 0.1, 0.1]]), df_complete=[4])
    assert r["df"][0] == pytest.approx(20 / 7)
    assert r["p"][0] == pytest.approx(0.09192389, abs=1e-6)              # 2 * pt(-2.5, 20/7)


def test_identical_decimal_estimates_count_as_identical():
    r = em.rubin(np.array([[0.1, 0.1, 0.1]]), np.array([[0.05, 0.05, 0.05]]), df_complete=[4])
    assert r["fmi_between"][0] == 0 and r["df"][0] == pytest.approx(20 / 7)
    assert r["p"][0] == pytest.approx(0.14388810, abs=1e-6)              # 2 * pt(-2, 20/7)


def test_degrees_of_freedom_do_not_jump_when_the_imputations_barely_disagree():
    same = em.rubin(np.array([[0.25, 0.25, 0.25]]), np.array([[0.1, 0.1, 0.1]]), df_complete=[10])
    near = em.rubin(np.array([[0.25, 0.25, 0.2500001]]), np.array([[0.1, 0.1, 0.1]]), df_complete=[10])
    assert same["df"][0] == pytest.approx(10 * 11 / 13) and near["df"][0] == pytest.approx(same["df"][0], rel=1e-6)


def test_a_fit_with_zero_standard_error_is_not_combined():
    r = em.rubin(np.array([[0.1, 0.3, 0.2]]), np.array([[0.0, 0.0, 0.0]]), df_complete=[4])
    assert np.isnan(r["p"][0])


def test_t_tail_is_accurate_far_in_the_tail():
    assert em.two_sided_p(10.0, 1000) == pytest.approx(1.66707e-22, rel=1e-3)
    assert em.two_sided_p(30.0, 1.01e6) > em.two_sided_p(30.0)            # t tails are heavier than the normal
    assert em.two_sided_p(1.96, 1e15) == pytest.approx(0.0499958, rel=1e-5)


def test_between_variance_that_swamps_a_tiny_within_variance_gives_p_near_one():
    r = em.rubin(np.array([[9.0, 10.0, 11.0]]), np.array([[1e-10] * 3]), df_complete=[4])
    assert r["p"][0] > 0.5


def test_t_tail_is_accurate_for_large_and_huge_degrees_of_freedom_and_small_p():
    assert em.two_sided_p(3.0, 1e12) == pytest.approx(0.0026998, rel=1e-4)
    assert em.two_sided_p(3.0, 1e5) == pytest.approx(0.0026998, rel=1e-3)
    assert em.two_sided_p(10.0, 5) == pytest.approx(1.7105e-4, rel=1e-3)
    assert em.two_sided_p(40.0, 30) < 1e-25
    est = np.full((1, 10), 0.3) + np.arange(10) * 1e-15
    assert em.rubin(est, np.full((1, 10), 0.1))["p"][0] == pytest.approx(0.0026998, rel=1e-3)


def test_df_is_recovered_from_a_t_test():
    p = em.two_sided_p(2.5, 7.0)
    assert em.df_from_p(2.5, p) == pytest.approx(7.0, rel=1e-4)
    assert em.df_from_p(2.5, em.two_sided_p(2.5)) == float("inf")


def test_rubin_rules_without_between_imputation_variance_return_the_single_fit():
    r = em.rubin(np.array([[0.5, 0.5, 0.5]]), np.array([[0.2, 0.2, 0.2]]))
    assert r["estimate"][0] == pytest.approx(0.5) and r["se"][0] == pytest.approx(0.2)
    assert r["p"][0] == pytest.approx(0.0124193, rel=1e-5)          # two-sided normal p of z = 2.5


def test_rubin_rules_add_the_between_imputation_variance():
    est = np.array([[1.0, 2.0, 3.0]]); se = np.array([[1.0, 1.0, 1.0]])
    r = em.rubin(est, se)
    # T = W + (1 + 1/M) B = 1 + (4/3) * 1
    assert r["estimate"][0] == pytest.approx(2.0)
    assert r["se"][0] == pytest.approx(np.sqrt(1 + 4 / 3))
    assert r["fmi_between"][0] == pytest.approx((4 / 3) / (1 + 4 / 3))
    assert np.isfinite(r["df"][0]) and r["df"][0] > 0


def test_rubin_rules_leave_a_feature_missing_from_any_fit_as_nan():
    r = em.rubin(np.array([[1.0, np.nan]]), np.array([[0.5, 0.5]]))
    assert np.isnan(r["estimate"][0]) and np.isnan(r["p"][0])


@pytest.mark.parametrize("z, df, expected", [(1.959964, float("inf"), 0.05), (2.570582, 5, 0.05),
                                              (12.7062, 1, 0.05), (2.228139, 10, 0.05), (0.0, 3, 1.0)])
def test_two_sided_p_matches_the_normal_and_t_tables(z, df, expected):
    assert em.two_sided_p(z, df) == pytest.approx(expected, rel=1e-4, abs=1e-12)


def test_each_component_must_converge_on_its_own():
    out = em.expected_counts({("A",): 1, ("A", "B"): 1_000_000_000, ("C",): 500_000_008, ("D",): 500_000_000,
                              ("C", "D"): 500_000_000})
    assert (not out.converged) or out.counts["B"] < 1.0


def test_one_read_set_shared_by_ten_candidates_is_drawn_from_the_exact_posterior():
    cls = {("A",): 1, tuple("ABCDEFGHIJ"): 1_000_000}
    share = np.array([d.get("A", 0) - 1 for d in em.imputations(cls, em.expected_counts(cls), 400,
                                                                  np.random.default_rng(1))]) / 1e6
    assert share.mean() == pytest.approx(0.25, abs=0.03) and share.std() == pytest.approx(0.164, abs=0.03)


def test_cauchy_tail_when_the_complement_rounds_to_one():
    assert em.two_sided_p(1e9, 1) == pytest.approx(6.3662e-10, rel=1e-4)


@pytest.mark.parametrize("classes", [
    {("A",): 1, ("A", "B"): 2_000_000_000, ("C",): 1_000_000_001, ("A", "B", "C"): 2_000_000_000},
    {("A",): 1, ("A", "B"): 200_000_000, ("C",): 100_000_001, ("A", "B", "C"): 200_000_000}])
def test_a_dominated_candidate_gets_exactly_its_mle_of_zero(classes):
    out = em.expected_counts(classes)
    assert out.converged and out.counts["B"] == 0.0
    assert sum(out.counts.values()) == pytest.approx(sum(classes.values()), rel=1e-12)


def test_candidates_with_identical_classes_are_not_dominated():
    assert em.dominated({("A", "B"): 10, ("A", "B", "C"): 5, ("C",): 2}) == set()
    assert em.dominated({("A", "B"): 10, ("A",): 1}) == {"B"}


def test_nested_classes_give_the_exact_mle():
    out = em.expected_counts({("A", "B"): 100_000_000, ("A", "B", "C"): 1, ("A",): 2, ("B",): 1, ("C",): 1})
    assert out.counts["A"] == pytest.approx(66_666_669.333, abs=1e-2)


def test_nested_classes_are_imputed_from_the_exact_posterior():
    cls = {("A", "B", "C"): 1_000_000, ("A", "B"): 100}
    c = np.array([d.get("C", 0) for d in em.imputations(cls, em.expected_counts(cls), 400, np.random.default_rng(7))])
    assert c.mean() == pytest.approx(4926, abs=1500)              # Beta(0.5, 101) x 1e6; MC se ~340
    cls = {("A",): 100, ("A", "B"): 1_000_000, ("A", "B", "C"): 1}
    b = np.array([d.get("B", 0) for d in em.imputations(cls, em.expected_counts(cls), 400, np.random.default_rng(2))])
    assert b.mean() == pytest.approx(4951, abs=1500)


def test_a_crossing_class_with_one_read_does_not_bias_the_imputations():
    cls = {("A", "B", "C"): 1_000_000, ("A", "B"): 100, ("B", "C"): 1}
    c = np.array([d.get("C", 0) for d in em.imputations(cls, em.expected_counts(cls), 200, np.random.default_rng(2))])
    assert c.mean() == pytest.approx(4974, abs=1500)


@pytest.mark.parametrize("classes", [
    {("A", "B"): 1_000_000_000, ("C", "D"): 500_000_000, ("A", "C"): 20, ("B", "D"): 1},
    {("A", "B", "E", "F"): 5_000_000, ("B", "E"): 50_000_000, ("A", "C"): 110, ("B", "F"): 15_000_000,
     ("C", "E", "F"): 10}])
def test_a_converged_fit_is_certified_close_to_the_maximum(classes):
    f = em.expected_counts(classes)
    c = em._Classes(classes)
    theta = np.array([f.theta[t] for t in c.targets])
    s = c.sums(theta)
    bound = min(em._gap(c, theta, s), em._dual_gap(c, theta, s, em._grad(c, theta, s), em.DUAL_TOL))
    assert f.converged and bound <= em.DUAL_TOL


def test_crossing_imputations_follow_the_information_of_a_single_crossing_read():
    cls = {("A", "B", "C"): 1_000_000, ("A", "B"): 100, ("B", "C"): 1}
    d = em.imputations(cls, em.expected_counts(cls), 200, np.random.default_rng(2))
    a = np.mean([x.get("A", 0) for x in d])
    assert a == pytest.approx(251_197, abs=60_000)                      # Beta(0.5, 1.5) x ~1e6; MC se ~18k


def test_rounded_shares_keep_their_sum():
    out = em.round_preserving({f"t{i}": 1 / 6 for i in range(6)}, 6)
    assert round(sum(out.values()), 6) == 1.0 and all(abs(v - 1 / 6) < 1e-6 for v in out.values())
