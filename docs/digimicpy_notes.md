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

## 3. With a uniform reference level the proxy is `gross - K / R0`; signs and rankings can still depend on `R0`

With a uniform `R0` and uniform leakage, the gross retained fraction equals
1 minus the leakage and the net proxy is `gross - K / R0`, where `K` is a weighted
mean of `mortality_i / sum_j uptake_ij` that does not depend on `R0` (for both
the biomass and the flux weighting). This does not make signs or treatment
rankings independent of `R0`:

* the proxy is negative whenever `R0 < K / gross`, so its sign can change with `R0`;
* two treatments with different gross retention and different `K` can swap order,
  at `R0 = (K_a - K_b) / (gross_a - gross_b)` when that value is positive;
* only for treatments with the same leakage (same gross retention) is the order
  fixed, by `K` alone.

In the 1,800 stored runs of this repository (`results/main_reference_consistency.csv`,
biomass-weighted proxy, R0 = 0.5, 1 and 2, comparisons within each community):

| Check | Result |
|---|---|
| Runs with a negative proxy at R0 = 0.5, 1, 2 | 406, 53, 0 of 1,799 with survivors |
| Runs whose sign changes across R0 | 406 of 1,799 (339 of them at leakage 0.8) |
| Treatment pairs at equal leakage whose order changes | 0 of 4,495 |
| Treatment pairs at different leakage whose order changes | 676 of 21,576 (3.1 %) |

For example, the mean proxy without leakage at 45 °C is below the mean at
leakage 0.2 and 21 °C for R0 = 0.5 (0.51 against 0.54) and above it for
R0 = 1 (0.76 against 0.67). With the flux weighting the counts are 400 sign
changes and 569 reversals (2.6 %). A uniform-`R0` sensitivity check is therefore
informative, and qualitative statements should be limited to those that hold at
every `R0` examined. Non-uniform reference vectors (for example the realised
resource levels of each community) would address an additional question and
were not examined. Biomass-weighted and flux-weighted summaries coincide closely
when leakage is uniform and uptake does not vary much among survivors.

## 4. Endpoint residual alone can miss slowly declining survivors

`micrm_stability.md` asks that the endpoint residual be small. The absolute
residual can be small while some survivors are still declining at a small
per-capita rate. In 600 runs of 20 communities (20 consumers, 10 resources),
integrating to t = 5000 changed richness in 19 % of runs compared with t = 60000
(3.7 % at t = 20000, 0.8 % at t = 120000). The per-capita rate of change of the
survivors, `|dC_i/dt| / C_i`, together with a check of the integration horizon,
detects this. The package examples integrate to t = 20, which is illustrative;
a note that the horizon has to be checked for each parameter set would help.
`examples/slow_transient.py` reproduces one such case in about 2 seconds
(output in `examples/slow_transient_output.txt`): at t = 5000 the largest
absolute residual is 3.0e-08, yet one of the 5 consumers above the threshold is
declining at a per-capita rate of 1.9e-03 and is lost by t = 20000.

## 5. Extinction threshold

The recipe uses `1e-5` as the survivor threshold. Between `1e-8` and `1e-3`,
richness changed in 0.5–1.2 % of runs, by one consumer in each case. Reporting
the threshold, as the recipe already advises, is sufficient; counting consumers
near the threshold is a cheap extra check.

See `examples/paired_temperature_leakage.py` for a compact example that
includes the diagnostics above.

A proposed documentation patch for notes 2 to 4 (`docs/content/cue.md` and
`docs/content/micrm_stability.md`) is in the branch
`docs/net-retention-proxy-and-transient-checks` of the fork
<https://github.com/zHxise/DigiMicPy>. It has not been proposed upstream.
