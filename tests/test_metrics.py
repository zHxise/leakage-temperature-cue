import numpy as np

from ltcue import ExperimentConfig, build_community, parameters_at
from ltcue.metrics import community_cue, reference_cue


def test_reference_cue_matches_tensor_formulation():
    cfg = ExperimentConfig()
    community = build_community(cfg, 0)
    p = parameters_at(community, cfg, 0.4, 300.0)
    eta = 1.0 - p.leakage.sum(axis=2)
    r0 = np.ones(cfg.n_resources)
    expected = ((p.uptake * eta * r0).sum(axis=1) - p.mortality) / (p.uptake * r0).sum(axis=1)
    np.testing.assert_allclose(reference_cue(p), expected, rtol=1e-12)


def test_cue_without_mortality_is_one_minus_leakage():
    cfg = ExperimentConfig()
    p = parameters_at(build_community(cfg, 0), cfg, 0.3, 303.0)
    gross = reference_cue(p) + p.mortality / (p.uptake @ np.ones(cfg.n_resources))
    np.testing.assert_allclose(gross, 0.7, rtol=1e-12)


def test_community_cue_weighting_and_empty():
    values = np.array([0.2, 0.6, 0.9])
    biomass = np.array([1.0, 3.0, 0.0])
    alive = np.array([True, True, False])
    assert np.isclose(community_cue(values, biomass, alive, "biomass"), 0.5)
    assert np.isclose(community_cue(values, biomass, alive, "equal"), 0.4)
    assert np.isnan(community_cue(values, biomass, np.zeros(3, bool), "biomass"))
