import numpy as np
import pandas as pd

from ltcue.analysis import bootstrap_mean_ci, leakage_slopes, slope_contrast, slope_summary
from ltcue.metrics import convergence_diagnostics
from ltcue import ExperimentConfig, build_community, parameters_at


def _table():
    rows = []
    for rep in range(6):
        for temp, slope in ((288.0, -1.0), (318.0, -3.0)):
            for leak in (0.0, 0.4, 0.8):
                rows.append({"rep": rep, "T": temp, "leak": leak, "y": 10 + rep + slope * leak})
    return pd.DataFrame(rows)


def test_bootstrap_ci_brackets_mean_and_handles_nan():
    mean, lo, hi = bootstrap_mean_ci(np.array([1.0, 2.0, 3.0, np.nan]))
    assert lo <= mean <= hi and np.isclose(mean, 2.0)
    assert np.isnan(bootstrap_mean_ci(np.array([np.nan]))[0])


def test_leakage_slopes_recover_known_effect_per_community():
    slopes = leakage_slopes(_table(), "y")
    assert len(slopes) == 12
    np.testing.assert_allclose(slopes[slopes["T"] == 288.0]["slope"], -1.0)
    np.testing.assert_allclose(slopes[slopes["T"] == 318.0]["slope"], -3.0)
    summary = slope_summary(slopes).set_index("T")
    assert np.isclose(summary.loc[318.0, "mean"], -3.0)


def test_slope_contrast_is_paired_difference():
    mean, lo, hi = slope_contrast(leakage_slopes(_table(), "y"), 288.0, 318.0)
    assert np.isclose(mean, -2.0) and lo <= mean <= hi


def test_convergence_diagnostics_zero_at_exact_equilibrium():
    cfg = ExperimentConfig(n_communities=1)
    p = parameters_at(build_community(cfg, 0), cfg, 0.0, 303.0)
    # with no consumers the resource equilibrium is s/d and every derivative is zero
    state = np.concatenate([np.zeros(cfg.n_consumers), p.resource_supply / p.resource_decay])
    out = convergence_diagnostics(state, p, np.zeros(cfg.n_consumers, dtype=bool))
    assert out["resid_abs"] < 1e-12 and out["rate_rel"] == 0.0


def test_claims_table_denominator_excludes_communities_without_survivors():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "scripts" / "reference_sensitivity.py"
    spec = importlib.util.spec_from_file_location("reference_sensitivity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    cfg = ExperimentConfig(leakage_fractions=(0.0, 0.8))
    rows = []
    for rep in range(4):
        for leak in cfg.leakage_fractions:
            for T in cfg.temperatures_k:
                value = np.nan if (rep == 0 and T == 318.0 and leak == 0.8) else -0.1 - 0.01 * rep
                rows.append({"R0": 1.0, "leak": leak, "T": T, "rep": rep,
                             "proxy": value, "proxy_flux": value})
    claims = module.claims_table(pd.DataFrame(rows), cfg)
    assert claims["n_with_survivors_at_max_leakage_45C"].iloc[0] == 3
    assert claims["share_negative_at_max_leakage_45C"].iloc[0] == 1.0


def test_reference_consistency_counts_sign_changes_and_reversals():
    from ltcue.analysis import reference_consistency

    # proxy = G - K / R0 for one community: A (G=1.0, K=0.4), B (G=0.8, K=0.1)
    # at different leakage, and C (G=1.0, K=0.45) at the same leakage as A.
    rows = []
    for r0 in (0.5, 1.0, 2.0):
        rows += [
            {"R0": r0, "rep": 0, "leak": 0.0, "T": 300.0, "proxy": 1.0 - 0.4 / r0},
            {"R0": r0, "rep": 0, "leak": 0.2, "T": 300.0, "proxy": 0.8 - 0.1 / r0},
            {"R0": r0, "rep": 0, "leak": 0.0, "T": 310.0, "proxy": 1.0 - 0.45 / r0},
        ]
    out = reference_consistency(pd.DataFrame(rows), "proxy")
    assert out["n_runs"] == 3
    assert out["n_sign_changes"] == 0 and out["n_negative_R0_0.5"] == 0
    # A and C share G, so their order cannot depend on R0
    assert out["n_pairs_equal_leakage"] == 1 and out["n_reversed_equal_leakage"] == 0
    # A and C are below B at R0 = 0.5 and 1 and above it at R0 = 2
    assert out["n_pairs_different_leakage"] == 2 and out["n_reversed_different_leakage"] == 2
