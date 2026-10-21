"""Outcome and convergence metrics for a simulated community."""

from __future__ import annotations

import numpy as np
from digimicpy import MiCRMParameters, micrm_rhs


def survivors(consumers: np.ndarray, threshold: float) -> np.ndarray:
    return consumers > threshold


def uptake_at_reference(parameters: MiCRMParameters, reference_resource: float) -> np.ndarray:
    """Uptake flux of each consumer at a fixed reference resource level R0."""

    r0 = np.full(parameters.n_resources, reference_resource)
    return parameters.uptake @ r0


def gross_retained_fraction(
    parameters: MiCRMParameters,
    reference_resource: float = 1.0,
) -> np.ndarray:
    """Per-consumer fraction of uptake retained after leakage.

        G_i = sum_j u_ij (1 - l_j) R0 / sum_j u_ij R0

    It does not involve mortality. With a leakage fraction shared by all
    resources it equals 1 - l for every consumer and every R0.
    """

    r0 = np.full(parameters.n_resources, reference_resource)
    uptake_flux = parameters.uptake @ r0
    retained_flux = (parameters.uptake * (1.0 - parameters.leakage_fraction)) @ r0
    with np.errstate(divide="ignore", invalid="ignore"):
        result = retained_flux / uptake_flux
    return np.where(uptake_flux > 0.0, result, np.nan)


def net_retention_proxy(
    parameters: MiCRMParameters,
    reference_resource: float = 1.0,
) -> np.ndarray:
    """Per-consumer net-retention proxy at a normalised reference resource level R0.

        P_i = ( sum_j u_ij (1 - l_j) R0 - m_i ) / ( sum_j u_ij R0 )

    The mortality term m_i is a generic first-order loss (maintenance, death,
    washout or a mixture), so P_i is a trait-based proxy evaluated at a
    normalised reference environment, not a realised carbon-use efficiency. At
    equilibrium a persisting consumer has zero net growth, so an equilibrium
    version of this quantity would be identically zero.
    """

    r0 = np.full(parameters.n_resources, reference_resource)
    uptake_flux = parameters.uptake @ r0
    retained_flux = (parameters.uptake * (1.0 - parameters.leakage_fraction)) @ r0
    with np.errstate(divide="ignore", invalid="ignore"):
        result = (retained_flux - parameters.mortality) / uptake_flux
    return np.where(uptake_flux > 0.0, result, np.nan)


def community_average(
    per_consumer: np.ndarray,
    consumers: np.ndarray,
    alive: np.ndarray,
    weights: np.ndarray | None = None,
) -> float:
    """Weighted mean over surviving consumers.

    ``weights`` defaults to consumer biomass. Pass ``consumers * uptake_flux`` for
    the flux-weighted version, which weights each survivor by its biomass times
    its uptake at the reference resource vector.
    """

    weight = consumers if weights is None else weights
    keep = alive & np.isfinite(per_consumer) & (weight > 0.0)
    if not keep.any():
        return float("nan")
    return float(np.average(per_consumer[keep], weights=weight[keep]))


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
