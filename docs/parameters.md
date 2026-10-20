# Parameters

All defaults are in `src/ltcue/config.py`. Energies are in eV, temperatures in
kelvin.

## Thermal traits

Rates follow the Sharpe–Schoolfield-type curve implemented in
`digimicpy.thermal_scaling_factor`, normalised to the rate at the reference
temperature (303 K, 30 °C).

| Process | Activation energy E (eV) | Deactivation energy Ed (eV) | Optimum Tpk (K) | Varies among consumers |
|---|---|---|---|---|
| Uptake | 0.65 | 3.0 | U(298, 315) | optimum only |
| Mortality | 0.30 | 3.0 | 340 (shared) | no |

* Only **uptake** and **mortality** depend on temperature.
* **Leakage fractions, resource supply and resource decay are temperature
  independent.**
* Because the mortality optimum (340 K) lies far above the simulated range,
  mortality increases monotonically with temperature, while uptake peaks near
  each consumer's own optimum and declines above it.
* The values above are **illustrative placeholders**, not fitted estimates.

## Community and environment

| Quantity | Value |
|---|---|
| Consumers / resources | 20 / 10 |
| Modules (uptake and leakage structure) | 2, specialisation ratio 8 |
| Per-consumer uptake scaling at reference temperature | U(0.5, 1.5) |
| Baseline mortality | U(0.15, 0.25) |
| Resource supply (per resource) | 1.0 |
| Resource decay (per resource) | 0.5 |
| Initial consumer biomass / resource | 0.01 / 1.0 |

## Treatments

| Factor | Levels |
|---|---|
| Temperature | 288, 294, 300, 306, 312, 318 K (15–45 °C) |
| Leakage fraction (same for all resources) | 0.0, 0.2, 0.4, 0.6, 0.8 |
| Communities | 15 in the pilot |

## Numerics and definitions

| Setting | Value |
|---|---|
| Solver | LSODA, rtol 1e-8, atol 1e-12 |
| Integration horizon | t = 0 to 60000 (the pilot script used 5000) |
| Extinction threshold | consumer biomass ≤ 1e-5 |
| CUE reference resource level | R0 = 1 for every resource |

