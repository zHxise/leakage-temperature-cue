"""Community construction.

A community is the set of traits that are fixed across every treatment: uptake
preferences and scaling, the structure of the leakage tensor, baseline
mortality and thermal optima. Temperature and leakage only transform these
fixed traits, so every treatment combination of one replicate is the *same*
community (paired design).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from digimicpy import (
    MiCRMParameters,
    generate_l_tensor,
    modular_uptake,
    temperature_adjusted_parameters,
)

from .config import ExperimentConfig

_STREAMS = ("uptake", "leakage", "mortality", "optima")


@dataclass(frozen=True)
class Community:
    """Treatment-independent traits of one replicate."""

    replicate: int
    uptake: np.ndarray            # (N, M) rates at the reference temperature
    leakage_unit: np.ndarray      # (N, M, M) tensor whose rows sum to 1
    mortality: np.ndarray         # (N,)
    uptake_topt: np.ndarray       # (N,) kelvin


def _generators(cfg: ExperimentConfig, replicate: int) -> dict[str, np.random.Generator]:
    if cfg.rng_mode == "legacy":
        shared = np.random.default_rng(cfg.base_seed + replicate)
        return {name: shared for name in _STREAMS}
    children = np.random.SeedSequence([cfg.base_seed, replicate]).spawn(len(_STREAMS))
    return {name: np.random.default_rng(seq) for name, seq in zip(_STREAMS, children)}


def build_community(cfg: ExperimentConfig, replicate: int) -> Community:
    """Draw the traits of one replicate.

    The leakage tensor is generated once with total leakage 1 and scaled by the
    treatment value afterwards. ``generate_l_tensor`` normalises each row before
    applying the total, so this is identical to generating it at the treatment
    value from the same stream.
    """

    rng = _generators(cfg, replicate)
    n, m = cfg.n_consumers, cfg.n_resources

    uptake = modular_uptake(n, m, cfg.n_modules, cfg.specialization_ratio, rng=rng["uptake"])
    low, high = cfg.uptake_scale_range
    uptake = uptake * rng["uptake"].uniform(low, high, (n, 1))

    leakage_unit = generate_l_tensor(
        n, m, cfg.n_modules, cfg.specialization_ratio, 1.0, rng=rng["leakage"]
    )

    low, high = cfg.mortality_range
    mortality = rng["mortality"].uniform(low, high, n)

    low, high = cfg.thermal.uptake_topt_range
    uptake_topt = rng["optima"].uniform(low, high, n)

    return Community(replicate, uptake, leakage_unit, mortality, uptake_topt)


def parameters_at(
    community: Community,
    cfg: ExperimentConfig,
    leakage: float,
    temperature_k: float,
) -> MiCRMParameters:
    """Model parameters of ``community`` under one treatment combination."""

    m = cfg.n_resources
    reference = MiCRMParameters(
        uptake=community.uptake,
        mortality=community.mortality,
        resource_supply=np.full(m, cfg.resource_supply),
        resource_decay=np.full(m, cfg.resource_decay),
        leakage=community.leakage_unit * leakage,
        leakage_fraction=np.full(m, leakage),
    )
    th = cfg.thermal
    return temperature_adjusted_parameters(
        reference,
        temperature_k,
        cfg.reference_temperature_k,
        uptake_activation_energy=th.uptake_ea,
        uptake_optimum_temperature=community.uptake_topt,
        uptake_deactivation_energy=th.uptake_ed,
        mortality_activation_energy=th.mortality_ea,
        mortality_optimum_temperature=th.mortality_topt,
        mortality_deactivation_energy=th.mortality_ed,
    )
