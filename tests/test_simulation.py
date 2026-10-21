import numpy as np
import pandas as pd
import pytest

from ltcue import ExperimentConfig, run_experiment

DATA = "tests/data/pilot_results.csv"


@pytest.fixture(scope="module")
def small():
    cfg = ExperimentConfig(
        n_communities=2,
        temperatures_k=(294.0, 312.0),
        leakage_fractions=(0.0, 0.4),
        t_end=2000.0,
    )
    return cfg, *run_experiment(cfg)


def test_shape_and_columns(small):
    cfg, table, states = small
    assert len(table) == 2 * 2 * 2
    assert states.shape == (len(table), cfg.n_consumers + cfg.n_resources)
    assert {"leak", "T", "rep", "richness", "biomass", "retention_proxy", "retention_proxy_flux", "gross_retained", "gross_retained_flux", "rate_rel"} <= set(table.columns)


def test_states_nonnegative_and_integrator_succeeds(small):
    _, table, states = small
    assert table["ok"].all()
    assert (states > -1e-9).all()


def test_deterministic(small):
    cfg, table, _ = small
    again, _ = run_experiment(cfg)
    pd.testing.assert_frame_equal(table, again)


def test_parallel_matches_serial(small):
    cfg, table, _ = small
    parallel, _ = run_experiment(cfg, n_jobs=2)
    pd.testing.assert_frame_equal(table, parallel)


@pytest.mark.slow
def test_legacy_mode_reproduces_pilot_results():
    cfg = ExperimentConfig(rng_mode="legacy", n_communities=2, t_end=5000.0)
    table, _ = run_experiment(cfg)
    pilot = pd.read_csv(DATA)
    keys = ["leak", "T", "rep"]
    merged = table.merge(pilot, on=keys, suffixes=("", "_pilot"))
    assert len(merged) == len(table)
    assert (merged["richness"] == merged["richness_pilot"]).all()
    np.testing.assert_allclose(merged["biomass"], merged["biomass_pilot"], rtol=1e-5, atol=1e-8)
    # the pilot stored its net-retention proxy (then called "cue")
    np.testing.assert_allclose(merged["retention_proxy"], merged["cue"], atol=1e-8)


def test_steady_state_mass_balance():
    # at steady state sum(m C) equals total supply minus resource decay flux
    from ltcue import build_community, simulate_point

    cfg = ExperimentConfig()
    row, state = simulate_point(cfg, build_community(cfg, 0), 0.4, 300.0)
    n = cfg.n_consumers
    from ltcue import parameters_at

    p = parameters_at(build_community(cfg, 0), cfg, 0.4, 300.0)
    lhs = float((p.mortality * state[:n]).sum())
    rhs = float(p.resource_supply.sum() - (p.resource_decay * state[n:]).sum())
    assert abs(lhs - rhs) / p.resource_supply.sum() < 1e-6
