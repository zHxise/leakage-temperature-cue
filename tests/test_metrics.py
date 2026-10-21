import numpy as np

from ltcue import ExperimentConfig, build_community, parameters_at
from ltcue.metrics import (
    community_average,
    gross_retained_fraction,
    net_retention_proxy,
    uptake_at_reference,
)


def _params(leak=0.4, temperature=300.0):
    cfg = ExperimentConfig()
    return cfg, parameters_at(build_community(cfg, 0), cfg, leak, temperature)


def test_proxy_matches_tensor_formulation():
    cfg, p = _params()
    eta = 1.0 - p.leakage.sum(axis=2)
    r0 = np.ones(cfg.n_resources)
    expected = ((p.uptake * eta * r0).sum(axis=1) - p.mortality) / (p.uptake * r0).sum(axis=1)
    np.testing.assert_allclose(net_retention_proxy(p), expected, rtol=1e-12)


def test_gross_retained_equals_one_minus_leakage_for_any_r0():
    _, p = _params(leak=0.3, temperature=303.0)
    for r0 in (0.5, 1.0, 2.0):
        np.testing.assert_allclose(gross_retained_fraction(p, r0), 0.7, rtol=1e-12)


def test_proxy_is_gross_minus_mortality_term_and_scales_with_r0():
    _, p = _params()
    for r0 in (0.5, 1.0, 2.0):
        expected = gross_retained_fraction(p, r0) - p.mortality / uptake_at_reference(p, r0)
        np.testing.assert_allclose(net_retention_proxy(p, r0), expected, rtol=1e-12)
    # the loss term halves when R0 doubles
    loss = lambda r0: gross_retained_fraction(p, r0) - net_retention_proxy(p, r0)
    np.testing.assert_allclose(loss(2.0), loss(1.0) / 2.0, rtol=1e-12)


def test_community_average_weightings_and_empty():
    values = np.array([0.2, 0.6, 0.9])
    biomass = np.array([1.0, 3.0, 0.0])
    alive = np.array([True, True, False])
    assert np.isclose(community_average(values, biomass, alive), 0.5)
    flux = biomass * np.array([3.0, 1.0, 5.0])
    assert np.isclose(community_average(values, biomass, alive, flux), (0.2 * 3 + 0.6 * 3) / 6)
    assert np.isnan(community_average(values, biomass, np.zeros(3, bool)))


def test_zero_uptake_gives_nan():
    _, p = _params()
    p0 = type(p)(uptake=np.zeros_like(p.uptake), mortality=p.mortality,
                 resource_supply=p.resource_supply, resource_decay=p.resource_decay,
                 leakage=p.leakage, leakage_fraction=p.leakage_fraction)
    assert np.isnan(net_retention_proxy(p0)).all()
