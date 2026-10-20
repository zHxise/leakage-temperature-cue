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
