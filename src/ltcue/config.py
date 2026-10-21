"""Experiment configuration.

Every quantity that defines an experiment lives here so that a run is fully
described by one serialisable object (stored next to the results).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from typing import Literal

import json


@dataclass(frozen=True)
class ThermalTraits:
    """Sharpe-Schoolfield style traits (energies in eV, temperatures in K).

    Uptake rates follow a performance curve with activation energy
    ``uptake_ea``, deactivation energy ``uptake_ed`` and a consumer-specific
    optimum drawn uniformly from ``uptake_topt_range``. Mortality follows the
    same functional form with a shared optimum placed above the simulated range,
    so it increases monotonically with temperature over the experiment.
    """

    uptake_ea: float = 0.65
    uptake_ed: float = 3.0
    uptake_topt_range: tuple[float, float] = (298.0, 315.0)
    mortality_ea: float = 0.30
    mortality_ed: float = 3.0
    mortality_topt: float = 340.0


@dataclass(frozen=True)
class ExperimentConfig:
    # Community structure
    n_consumers: int = 20
    n_resources: int = 10
    n_modules: int = 2
    specialization_ratio: float = 8.0
    uptake_scale_range: tuple[float, float] = (0.5, 1.5)
    mortality_range: tuple[float, float] = (0.15, 0.25)
    resource_supply: float = 1.0
    resource_decay: float = 0.5

    # Treatments
    temperatures_k: tuple[float, ...] = (288.0, 294.0, 300.0, 306.0, 312.0, 318.0)
    leakage_fractions: tuple[float, ...] = (0.0, 0.2, 0.4, 0.6, 0.8)
    reference_temperature_k: float = 303.0
    n_communities: int = 15

    # Randomness
    base_seed: int = 1000
    # "independent": one child stream per trait (uptake, leakage, mortality,
    # optima), so adding a trait never shifts the others.
    # "legacy": a single stream consumed in the order of the pilot script.
    rng_mode: Literal["independent", "legacy"] = "independent"

    # Dynamics
    initial_consumer_biomass: float = 0.01
    initial_resource: float = 1.0
    t_end: float = 60000.0
    method: str = "LSODA"
    rtol: float = 1e-8
    atol: float = 1e-12

    # Outcome definitions
    extinction_threshold: float = 1e-5
    reference_resource: float = 1.0

    thermal: ThermalTraits = field(default_factory=ThermalTraits)

    def with_(self, **changes) -> "ExperimentConfig":
        return replace(self, **changes)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "ExperimentConfig":
        raw = json.loads(text)
        raw["thermal"] = ThermalTraits(
            **{
                k: tuple(v) if isinstance(v, list) else v
                for k, v in raw["thermal"].items()
            }
        )
        for key in (
            "uptake_scale_range",
            "mortality_range",
            "temperatures_k",
            "leakage_fractions",
        ):
            raw[key] = tuple(raw[key])
        return cls(**raw)
