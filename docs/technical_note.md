# Technical note: temperature and leakage in microbial consumer-resource communities

Status: first run. Thermal parameters are illustrative placeholders
(`parameters.md`), so the results describe what the model does under these
assumptions and are not biological conclusions.

## Question and design

How do temperature and leakage jointly affect persistence, biomass and two
retention metrics when the same communities are compared across all
combinations?

60 communities of 20 consumers and 10 resources, each simulated at six
temperatures (15–45 °C) and five leakage fractions (0–0.8): 1,800 simulations,
integrated to t = 60000. Each community is built once; only the leakage scaling
and the temperature adjustment change between treatments. Uptake and mortality
depend on temperature; leakage, resource supply and resource decay do not.
Uncertainty is a 95 % bootstrap over communities.

## Metrics

* **Persistence:** number of consumers with final biomass above 1e-5.
* **Biomass:** total final consumer biomass.
* **Gross retained fraction** G: share of uptake that is retained after leakage,
  `sum_j u_ij (1 - l_j) R0 / sum_j u_ij R0`, averaged over surviving consumers.
  It does not involve mortality and is fixed by leakage: here it equals
  1 − leakage for every consumer, every R0 and both weightings.
* **Net-retention proxy** P: `G_i − m_i / sum_j u_ij R0`. The mortality term is a
  generic first-order loss (maintenance, death, washout or a mixture), so P is a
  trait-based proxy at a normalised reference environment, not a realised
  carbon-use efficiency. It is not an equilibrium quantity: at equilibrium a
  persisting consumer has zero net growth.
* **Weightings:** biomass-weighted, and flux-weighted (each survivor weighted by
  its biomass times its uptake at the reference resource vector, closer to a
  community-level summary of carbon allocation). Neither is presented as a
  realised CUE; that would need explicit accounting of biomass production and
  carbon fate.
* **Reference environment:** R0 = 1 for every resource, treated as a normalised
  reference, not as the realised equilibrium resource. For comparison, supply 1
  and decay 0.5 give a resource level of 2 without consumers. Results are also
  computed at R0 = 0.5 and 2. With a uniform R0 the proxy is exactly
  `P(R0) = G − K/R0`, with K the weighted mean of `m_i / sum_j u_ij`, and the
  weights do not depend on R0, so signs and orderings are the same for every R0
  and only magnitudes change; this is a property of the definition, not an
  empirical robustness test. A non-uniform reference vector (for example the
  realised resource levels) would be a different and informative check, not done here.

## Results (means over 60 communities)

**Persistence** (Fig. 2). Richness has a maximum at 33 °C (about 4.0 of 20 at
zero leakage) and is lower at 15 °C (3.3) and 45 °C (3.0). Leakage lowers
richness only slightly: from 0 to 0.8 leakage the mean loss is 0.15 to 0.35
consumers between 15 and 39 °C and 0.67 at 45 °C. The within-community slope of
richness on leakage is negative at 15, 21, 33, 39 and 45 °C (−0.40, −0.37, −0.33,
−0.28 and −0.73 consumers per unit leakage); its interval includes zero at 27 °C
(−0.15) (Fig. 4). The 45 °C slope is about 1.8 times the 15 °C slope, but the
paired contrast is −0.33 consumers per unit leakage with 95 % interval −0.64 to
−0.01, so the evidence for a stronger effect at high temperature is weak.

**Biomass** (Figs. 2 and 4). Total biomass decreases with temperature at every
leakage level (90 to 25 from 15 to 45 °C without leakage). At steady state the
product of biomass and the mortality scaling factor is nearly constant (49–52
between 15 and 39 °C), so the decline follows the rise in mortality. Leakage
reduces biomass at all temperatures, not linearly (at 15 °C: 90, 88, 84, 75 and 51
for leakage 0 to 0.8), so slopes per unit leakage are summaries only (−46 at
15 °C, −19 at 27 and 45 °C). The relative reduction from 0 to 0.8 leakage is 44 % at
15 °C, 29 % at 27 °C and 65 % at 45 °C.

**Retention metrics** (Figs. 3 and 5).
* Gross retained fraction: 1 − leakage, by construction.
* Net-retention proxy at R0 = 1: its decrease with leakage is exactly −1 per
  consumer (community slope −1.0 at every temperature, −0.98 at 45 °C because
  the set of survivors changes). Without leakage it rises from 0.83 at 15 °C to
  0.88 at 27 °C and falls to 0.76 at 45 °C (flux-weighted: 0.84, 0.88, 0.77); the two weightings differ by at most 0.01
  without leakage.
* Statements that hold at R0 = 0.5, 1 and 2 (`results/main_reference_claims.csv`;
  by construction, see above, so they are consequences of the definition):
  the proxy decreases with leakage at every temperature; without leakage it is
  highest at 27 °C among the simulated temperatures, and 27 °C is higher than
  both 21 and 33 °C at every leakage level (paired intervals above zero); it is
  lower at 15 °C than at 27 °C (−0.10, −0.05, −0.02 for R0 = 0.5, 1, 2, no
  leakage) and lower at 45 °C than at 39 °C (−0.20, −0.10, −0.05), with
  intervals excluding zero. The size of these differences scales with 1/R0.
* Not invariant: the sign and level of the proxy. At 0.8 leakage and 45 °C the
  proxy is negative in 59/59 (100 %), 35/59 (59 %) and 0/59 (0 %) of the
  communities with survivors for R0 = 0.5, 1 and 2; at 15 °C, in 60/60 (100 %),
  11/60 (18 %) and 0/60 (0 %). Statements about whether the proxy is positive or
  negative should not be made.

## Numerical checks (`results/sensitivity_summary.csv`, 20 communities, 600 runs each)

These runs were made before the metrics were renamed; the retention proxy equals
the earlier `cue` column (verified), so the values are unchanged.


* Tighter tolerances and a stiff solver (Radau) reproduce biomass to better than
  3e-9 (relative) and the proxy to better than 2e-11, with no change in
  persistence.
* Different initial biomass and resource levels give identical persistence.
* Alternative extinction thresholds (1e-8 to 1e-3) change persistence in 0.5–1.2 %
  of runs, by one consumer.
* Integrating to 120000 instead of 60000 changes persistence in 0.8 % of runs.
* At the earlier horizon (5000), 19 % of runs differed in persistence, and the sign
  of the mean leakage slope of richness changed at one or more temperatures,
  because some consumers were still declining slowly (20000: 3.7 % of runs).
  The default horizon is therefore 60000.
* No integrator failures in any run.

## Model assumptions and limitations

1. The equations describe a generic open-resource environment: constant resource
   input, first-order resource removal and no explicit dilution of consumers.
   This model is not described as a fermentation. A continuous-culture reading
   would need an explicit dilution rate D (consumers lose D·C_i, resources gain
   D·(R_in − R)), with any genuine chemical decay kept separate; a batch reading
   would have no continuous supply after inoculation. This extension is deferred
   until the intended experimental setting is chosen.
2. Thermal parameters are illustrative. The uptake curves are normalised at
   303 K within the range of optima (298–315 K), so the relative peak height
   depends on the distance between each optimum and 303 K (1.5 at 298 K, 2.0 at
   315 K), an implicit correlation between traits that is not intended. The
   mortality optimum (340 K) makes mortality rise monotonically over the range.
   The decline of biomass follows this rise; the lower proxy at 45 °C combines
   lower uptake of the survivors above their optima (about half of the drop
   between 39 and 45 °C without leakage) and higher mortality (about a third).
3. Leakage is temperature independent and equal for all resources and species.
4. Mortality is a generic loss term; interpretation as maintenance, death or
   washout cannot be identified within the model.
5. The proxy depends on the reference environment R0 through its level and
   magnitude differences (see above).
6. Communities are small (20 consumers) and few persist (2–4 on average), so the
   persistence effects are small relative to variation between communities.
7. Persistence is "persisting consumers over the integration horizon", not a
   stability or coexistence criterion.

## Reproducing the results

See the README. Figures are in `figures/` and tables in `results/`.
