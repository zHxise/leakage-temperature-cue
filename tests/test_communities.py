import numpy as np
import pytest

from ltcue import ExperimentConfig, build_community, parameters_at


@pytest.fixture
def cfg():
    return ExperimentConfig(n_communities=3)


def test_leakage_rows_sum_to_treatment_value(cfg):
    community = build_community(cfg, 0)
    for leak in cfg.leakage_fractions:
        p = parameters_at(community, cfg, leak, 303.0)
        np.testing.assert_allclose(p.leakage.sum(axis=2), leak, atol=1e-12)
        np.testing.assert_allclose(p.leakage_fraction, leak)


def test_traits_shared_across_leakage(cfg):
    community = build_community(cfg, 1)
    base = parameters_at(community, cfg, 0.0, 306.0)
    for leak in cfg.leakage_fractions[1:]:
        p = parameters_at(community, cfg, leak, 306.0)
        np.testing.assert_array_equal(p.uptake, base.uptake)
        np.testing.assert_array_equal(p.mortality, base.mortality)
        np.testing.assert_array_equal(p.resource_supply, base.resource_supply)
        np.testing.assert_array_equal(p.resource_decay, base.resource_decay)
        np.testing.assert_allclose(p.leakage, community.leakage_unit * leak, atol=1e-15)


def test_replicates_differ(cfg):
    a, b = build_community(cfg, 0), build_community(cfg, 1)
    assert not np.array_equal(a.uptake, b.uptake)


def test_build_is_deterministic(cfg):
    a, b = build_community(cfg, 2), build_community(cfg, 2)
    np.testing.assert_array_equal(a.uptake, b.uptake)
    np.testing.assert_array_equal(a.leakage_unit, b.leakage_unit)
    np.testing.assert_array_equal(a.uptake_topt, b.uptake_topt)


def test_independent_streams_decouple_traits(cfg):
    # Changing the number of resources changes how many
    # draws the uptake and leakage streams consume; mortality and optima must not move.
    a = build_community(cfg, 0)
    b = build_community(cfg.with_(n_resources=8), 0)
    np.testing.assert_array_equal(a.mortality, b.mortality)
    np.testing.assert_array_equal(a.uptake_topt, b.uptake_topt)


def test_legacy_stream_couples_traits(cfg):
    legacy = cfg.with_(rng_mode="legacy")
    a = build_community(legacy, 0)
    b = build_community(legacy.with_(n_resources=8), 0)
    assert not np.array_equal(a.mortality, b.mortality)


def test_traits_shared_across_temperature(cfg):
    community = build_community(cfg, 0)
    reference = parameters_at(community, cfg, 0.4, 300.0)
    for temperature in cfg.temperatures_k:
        p = parameters_at(community, cfg, 0.4, temperature)
        ratio = p.uptake / community.uptake
        # every row of uptake is rescaled by a single per-consumer factor
        np.testing.assert_allclose(ratio, ratio[:, [0]] * np.ones_like(ratio), rtol=1e-12)
        # mortality is rescaled by the same factor for every consumer
        mort_ratio = p.mortality / community.mortality
        np.testing.assert_allclose(mort_ratio, mort_ratio[0], rtol=1e-12)
        np.testing.assert_array_equal(p.leakage, reference.leakage)
        np.testing.assert_array_equal(p.resource_supply, reference.resource_supply)
        np.testing.assert_array_equal(p.resource_decay, reference.resource_decay)
