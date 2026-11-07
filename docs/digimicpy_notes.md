# Notes on DigiMicPy documentation and behaviour

Observations from building the temperature x leakage comparison, checked against
DigiMicPy at commit `ac9b9fe` (the version pinned in `requirements.txt`). They
are candidates for documentation updates; none was verified beyond this setting.

## 1. Thermal scaling is normalised at the reference temperature

`thermal_scaling_factor` divides the performance at the target temperature by the
performance at `reference_temperature`, and `temperature_adjusted_parameters`
treats the supplied parameters as the rates at that reference. Consequently the
height of the thermal peak relative to the reference rate depends on how far each
trait's optimum is from the reference temperature. With activation energy 0.65 eV,
deactivation energy 3.0 eV and a reference of 303 K, the factor at the optimum is:

| Optimum (K) | 298 | 303 | 308 | 315 |
|---|---|---|---|---|
| Factor at the optimum | 1.50 | 1.00 | 1.22 | 2.03 |

When optima are drawn independently across consumers, peak height therefore
varies with the optimum even if no such correlation is intended. A sentence in
`thermal.md` would make this explicit.

## 2. The CUE recipe calls mortality "maintenance"

`cue.md` describes `(retained0 - mortality) / uptake0` as maintenance-adjusted
CUE. In the model, mortality is a generic first-order loss of consumer biomass and
cannot be told apart from death or washout. A note that the quantity is a
reference-condition net-retention proxy, and not a realised CUE, would avoid
over-interpretation.

## 3. A uniform reference resource level makes the sensitivity check uninformative

With a uniform `R0` and uniform leakage, the gross retained fraction equals
1 minus the leakage and the net proxy is `gross - K / R0`, where `K` is a weighted
mean that does not depend on `R0`. Signs and orderings across treatments are then
the same for every `R0`; only magnitudes change. A sensitivity check on `R0` is
informative only with non-uniform reference vectors (for example the realised
resource levels of each community). Biomass-weighted and flux-weighted summaries
also coincide closely when leakage is uniform and uptake does not vary much among
survivors.

## 4. Endpoint residual alone can miss slowly declining survivors

`micrm_stability.md` asks that the endpoint residual be small. The absolute
residual can be small while some survivors are still declining at a small
per-capita rate. In 600 runs of 20 communities (20 consumers, 10 resources),
integrating to t = 5000 changed richness in 19 % of runs compared with t = 60000
(3.7 % at t = 20000, 0.8 % at t = 120000). The per-capita rate of change of the
survivors, `|dC_i/dt| / C_i`, together with a check of the integration horizon,
detects this. The package examples integrate to t = 20, which is illustrative;
a note that the horizon has to be checked for each parameter set would help.

## 5. Extinction threshold

The recipe uses `1e-5` as the survivor threshold. Between `1e-8` and `1e-3`,
richness changed in 0.5–1.2 % of runs, by one consumer in each case. Reporting
the threshold, as the recipe already advises, is sufficient; counting consumers
near the threshold is a cheap extra check.

See `examples/paired_temperature_leakage.py` for a compact example that
includes the diagnostics above.
