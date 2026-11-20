"""combine_imputations: which ANCOM-BC2 fits count, and Rubin / Barnard-Rubin pooling checked against R.

Reference numbers come from R 4.5.3 (pt, p.adjust); the ANCOM-BC2 behaviour mimicked here is ANCOMBC 2.14.0
.ancombc2_core: p_hat = 2 * pt(abs(W), df = dof, lower.tail = FALSE); p_hat[is.na(p_hat)] = 1, with dof the lm
residual df (0 when the feature is nonzero in no more samples than the design has columns, giving p = NaN -> 1).
"""
import numpy as np
import pandas as pd
import pytest

from unclassified_meta import da
from unclassified_meta.runner import StepError

T = "groupB"


def fit(lfc, se, p=None, w="auto", ids=None, robust="TRUE"):
    """One ANCOM-BC2 result table for term groupB; w='auto' gives W = lfc / se, None leaves the W column out."""
    lfc, se = np.asarray(lfc, float), np.asarray(se, float)
    ids = ids or [f"f{i + 1}" for i in range(len(lfc))]
    d = {"feature_id": ids, "label": ids, "lfc_(Intercept)": 0.0, "se_(Intercept)": 1.0,
         f"lfc_{T}": lfc, f"se_{T}": se}
    if w is not None:
        d[f"W_{T}"] = lfc / se if isinstance(w, str) else np.asarray(w, float)
    if p is not None:
        d[f"p_{T}"] = np.asarray(p, float)
    d.update({f"diff_{T}": "TRUE", f"passed_ss_{T}": "TRUE", f"diff_robust_{T}": robust})
    return pd.DataFrame(d)


def pool(fits, method="holm", alpha=0.05):
    out, _ = da.combine_imputations(fits, alpha, method)
    return out.set_index("feature_id")


# --- which fits count ---------------------------------------------------------------------------------------------

def test_non_testable_p1_with_nonzero_w_in_every_fit_is_reported_untestable_and_never_significant():
    fits = [fit([4.3, 0.5], [0.2, 0.5], p=[1.0, 0.33]) for _ in range(5)]
    row = pool(fits, method="none")
    assert row.loc["f1", "status_groupB"] == "untestable"
    assert row.loc["f1", "n_fits_groupB"] == 0 and row.loc["f1", "n_untestable_groupB"] == 5
    assert row.loc["f1", "p_groupB"] == 1.0 and row.loc["f1", "q_groupB"] == 1.0
    assert np.isnan(row.loc["f1", "lfc_groupB"]) and not row.loc["f1", "diff_groupB"]


def test_valid_p1_with_zero_w_is_an_ordinary_test_and_is_pooled():
    fits = [fit([0.0, 1e-14], [0.5, 0.5], p=[1.0, 1.0]) for _ in range(4)]   # |W| = 2e-14: p rounds to 1
    row = pool(fits)
    assert (row["status_groupB"] == "pooled").all() and (row["n_fits_groupB"] == 4).all()
    assert (row["n_untestable_groupB"] == 0).all() and np.allclose(row["p_groupB"], 1.0)


def test_p1_with_missing_w_is_judged_from_lfc_over_se():
    a = [fit([3.0], [0.5], p=[1.0], w=[np.nan]) for _ in range(3)]          # W written as NA
    b = [fit([3.0], [0.5], p=[1.0], w=None) for _ in range(3)]              # no W column at all
    c = [fit([0.0], [0.5], p=[1.0], w=None) for _ in range(3)]              # lfc 0: a real p = 1
    assert pool(a).loc["f1", "status_groupB"] == "untestable"
    assert pool(b).loc["f1", "status_groupB"] == "untestable"
    assert pool(c).loc["f1", "status_groupB"] == "pooled"


def test_partially_untestable_feature_is_not_pooled_over_the_testable_subset():
    # testable in 2 of 5 fits with a huge effect; pooling those 2 would give p ~ 0
    fits = [fit([5.0], [0.1], p=[1e-30]) for _ in range(2)] + [fit([5.0], [0.1], p=[1.0]) for _ in range(3)]
    row = pool(fits, method="none")
    assert row.loc["f1", "status_groupB"] == "partly_untestable"
    assert row.loc["f1", "n_fits_groupB"] == 2 and row.loc["f1", "n_untestable_groupB"] == 3
    assert row.loc["f1", "p_groupB"] == 1.0 and not row.loc["f1", "diff_groupB"]


def test_missing_or_failed_p_makes_a_fit_not_count():
    fits = [fit([2.0], [0.5], p=[1e-3]) for _ in range(3)] + [fit([2.0], [0.5], p=[np.nan])]
    row = pool(fits, method="none")
    assert row.loc["f1", "status_groupB"] == "incomplete" and row.loc["f1", "p_groupB"] == 1.0


@pytest.mark.parametrize("bad_se", [0.0, np.inf, np.nan, -0.5])
def test_zero_infinite_or_missing_variance_in_one_fit_keeps_the_feature_unpooled(bad_se):
    fits = [fit([2.0, 2.0], [0.5, 0.5], p=[1e-3, 1e-3]) for _ in range(3)]
    fits.append(fit([2.0, 2.0], [0.5, bad_se], p=[1e-3, 1e-3]))
    row = pool(fits, method="none")
    assert row.loc["f1", "status_groupB"] == "pooled"
    assert row.loc["f2", "status_groupB"] == "incomplete" and row.loc["f2", "p_groupB"] == 1.0


def test_zero_between_imputation_variance_reproduces_the_complete_data_t_test():
    # identical fits, W = 4, p from a t test with 3.5 df (R: 2 * pt(4, 3.5, lower.tail = FALSE))
    p = 0.0209206163654037
    row = pool([fit([2.0], [0.5], p=[p]) for _ in range(4)], method="none")
    assert row.loc["f1", "fmi_groupB"] == 0 and row.loc["f1", "se_groupB"] == pytest.approx(0.5)
    assert row.loc["f1", "df_groupB"] == pytest.approx(3.5, rel=1e-6)
    assert row.loc["f1", "p_groupB"] == pytest.approx(p, rel=1e-8)


def test_feature_absent_from_every_fit_is_not_estimated_and_absent_from_some_is_incomplete():
    fits = [fit([1.0, 1.0, np.nan], [0.5, 0.5, np.nan], p=[0.05, 0.05, np.nan]) for _ in range(2)]
    fits.append(fit([1.0, np.nan], [0.5, np.nan], p=[0.05, np.nan], ids=["f1", "f2"]))
    st = pool(fits)["status_groupB"]
    assert st["f1"] == "pooled" and st["f2"] == "incomplete" and st["f3"] == "not_estimated"


# --- malformed fit sets ------------------------------------------------------------------------------------------

def test_a_fit_missing_a_coefficient_column_is_refused():
    fits = [fit([1.0], [0.5], p=[0.05]) for _ in range(3)]
    fits[2] = fits[2].drop(columns=[f"se_{T}"])
    with pytest.raises(StepError, match="fit 3 does not have the coefficients"):
        da.combine_imputations(fits, 0.05, "holm")


def test_a_fit_with_an_extra_term_is_refused():
    fits = [fit([1.0], [0.5], p=[0.05]) for _ in range(3)]
    fits[1] = fits[1].assign(lfc_groupC=1.0, se_groupC=0.5)
    with pytest.raises(StepError, match="extra lfc_groupC"):
        da.combine_imputations(fits, 0.05, "holm")


@pytest.mark.parametrize("n", [0, 1])
def test_fewer_than_two_fits_are_refused(n):
    with pytest.raises(StepError, match="at least 2 fits"):
        da.combine_imputations([fit([1.0], [0.5], p=[0.05]) for _ in range(n)], 0.05, "holm")


# --- numbers against R ---------------------------------------------------------------------------------------------

def test_rubin_and_barnard_rubin_match_r():
    # est 1, 2, 3, se 1: W = 1, B = 1, T = 7/3, lambda = 4/7, nu_old = 6.125
    est, se = [1.0, 2.0, 3.0], [1.0, 1.0, 1.0]
    p10 = [0.34089313230206, 0.0733880347707404, 0.0133436550225696]   # 2 * pt(est / se, 10, lower = FALSE)
    nc = pool([fit([e], [s], p=[p]) for e, s, p in zip(est, se, p10)], method="none").loc["f1"]
    assert nc["lfc_groupB"] == pytest.approx(2.0)
    assert nc["se_groupB"] == pytest.approx(np.sqrt(7 / 3))
    assert nc["fmi_groupB"] == pytest.approx(4 / 7)
    assert nc["df_groupB"] == pytest.approx(2.27778560360614, rel=1e-6)          # Barnard-Rubin, nu_com = 10
    assert nc["p_groupB"] == pytest.approx(0.306871794179897, rel=1e-6)
    # p far below any t(W = 1..3) tail: df_from_p gives the normal complete-data reference, nu = (M - 1) / lambda^2
    inf = pool([fit([e], [s], p=[1e-300]) for e, s in zip(est, se)], method="none").loc["f1"]
    assert inf["df_groupB"] == pytest.approx(6.125)
    assert inf["p_groupB"] == pytest.approx(0.237400350431343, rel=1e-6)


def test_bh_and_holm_after_pooling_match_r_p_adjust():
    p = [2.148356e-11, 1.077986e-05, 1.694492e-05, 6.859258e-05, 7.169939e-05, 0.02, 0.02, 1, 1, 1]
    bh = [2.148356e-10, 5.38993e-05, 5.64830666666667e-05, 1.4339878e-04, 1.4339878e-04, 0.0285714285714286,
          0.0285714285714286, 1, 1, 1]
    holm = [2.148356e-10, 9.701874e-05, 1.3555936e-04, 4.8014806e-04, 4.8014806e-04, 0.1, 0.1, 1, 1, 1]
    assert np.allclose(da.p_adjust(p, "BH"), bh, rtol=1e-9)
    assert np.allclose(da.p_adjust(p, "holm"), holm, rtol=1e-7)


def test_untestable_features_stay_in_the_bh_denominator():
    # one real signal and 9 untestable features: BH over 10 p values, not over the single pooled one
    fits = [fit([3.0] + [5.0] * 9, [1.0] * 10, p=[1e-3] + [1.0] * 9) for _ in range(3)]
    row = pool(fits, method="BH")
    pooled_p = row.loc["f1", "p_groupB"]
    assert row.loc["f1", "q_groupB"] == pytest.approx(min(1, 10 * pooled_p))
    assert (row.drop("f1")["q_groupB"] == 1).all()


def test_pooled_table_is_internally_consistent():
    rng = np.random.default_rng(7)
    n, m = 60, 6
    fits = []
    for _ in range(m):
        lfc = rng.normal(0, 2, n)
        se = rng.uniform(0.2, 1.5, n)
        p = rng.uniform(1e-6, 0.9, n)
        p[:10] = 1.0                                          # untestable in every fit
        p[10:20] = np.where(rng.random(10) < 0.5, 1.0, 0.01)  # untestable in some
        se[20:23] = np.where(rng.random(3) < 0.3, 0.0, se[20:23])
        fits.append(fit(lfc, se, p=p))
    out, terms = da.combine_imputations(fits, 0.05, "BH")
    assert terms == [T]
    assert set(out[f"status_{T}"]) <= set(da.MI_STATUS)
    pooled = out[f"status_{T}"].eq("pooled")
    assert (pooled == (out[f"n_fits_{T}"] == m)).all()
    assert (out.loc[~pooled, f"p_{T}"] == 1).all() and out.loc[~pooled, f"lfc_{T}"].isna().all()
    assert (out[f"status_{T}"].iloc[:10] == "untestable").all()
    assert np.allclose(out[f"q_{T}"], da.p_adjust(out[f"p_{T}"], "BH"))
    assert not out.loc[~pooled, f"diff_{T}"].any() and not out.loc[~pooled, f"diff_robust_{T}"].any()
    assert ((out.loc[pooled, f"fmi_{T}"] >= 0) & (out.loc[pooled, f"fmi_{T}"] < 1)).all()
    assert (out[f"n_fits_{T}"] + out[f"n_untestable_{T}"] <= m).all()


def test_da_adjusts_with_bh_by_default():
    from unclassified_meta.cli import build_parser
    args = build_parser().parse_args(["da", "--counts", "c.tsv", "--features", "f.tsv", "--metadata", "m.tsv",
                                      "--group", "g"])
    assert args.p_adjust == "BH"
    script = da.SCRIPT.read_text()
    assert '`p-adjust` = "BH"' in script                     # the R script's own default agrees
