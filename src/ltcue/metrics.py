"""Outcome and convergence metrics for a simulated community."""

from __future__ import annotations

import numpy as np
from digimicpy import MiCRMParameters, micrm_rhs


def survivors(consumers: np.ndarray, threshold: float) -> np.ndarray:
    return consumers > threshold


def reference_cue(
    parameters: MiCRMParameters,
    reference_resource: float = 1.0,
) -> np.ndarray:
    """Per-consumer carbon-use efficiency at a fixed reference resource level.

    For consumer i with uptake u_ij, retained fraction (1 - l_j) and mortality
    m_i (treated as maintenance loss):

        CUE_i = (sum_j u_ij (1 - l_j) R0 - m_i) / (sum_j u_ij R0)

    Evaluating at a fixed R0 (not at the equilibrium resource pool) keeps the
    quantity a property of the traits. At equilibrium a persisting consumer has
    net growth zero by construction, so an equilibrium version would be
    identically zero.
    """

    r0 = np.full(parameters.n_resources, reference_resource)
    uptake_flux = parameters.uptake @ r0
    retained_flux = (parameters.uptake * (1.0 - parameters.leakage_fraction)) @ r0
    with np.errstate(divide="ignore", invalid="ignore"):
        result = (retained_flux - parameters.mortality) / uptake_flux
    return np.where(uptake_flux > 0.0, result, np.nan)


def community_cue(
    per_consumer: np.ndarray,
    consumers: np.ndarray,
    alive: np.ndarray,
    weighting: str = "biomass",
) -> float:
    if not alive.any():
        return float("nan")
    if weighting == "biomass":
        return float(np.average(per_consumer[alive], weights=consumers[alive]))
    if weighting == "equal":
        return float(per_consumer[alive].mean())
    raise ValueError(f"unknown weighting: {weighting}")


def convergence_diagnostics(
    state: np.ndarray,
    parameters: MiCRMParameters,
    alive: np.ndarray,
) -> dict[str, float]:
    """Residuals of the right-hand side at the final state.

    ``resid_abs`` is the largest absolute derivative over all variables;
    ``rate_rel`` is the largest per-capita rate of change among survivors, which
    is scale free and therefore comparable across treatments.
    """

    n = parameters.n_consumers
    rhs = micrm_rhs(0.0, state, parameters)
    consumers = state[:n]
    rate_rel = (
        float(np.max(np.abs(rhs[:n][alive]) / consumers[alive])) if alive.any() else 0.0
    )
    return {"resid_abs": float(np.max(np.abs(rhs))), "rate_rel": rate_rel}
