# ltcue: temperature × leakage in microbial consumer-resource communities

Simulations of how temperature and metabolite leakage jointly affect
**persistence** (number of surviving consumers), **total biomass** and two
**reference-condition retention metrics**, built on
[DigiMicPy](https://github.com/DigiMicOrg/DigiMicPy).

The environment is a generic open-resource system (constant input, first-order
resource removal, no consumer dilution), not a fermentation. Thermal parameters
are illustrative placeholders; see [`docs/parameters.md`](docs/parameters.md).
Model assumptions, results and limitations are in
[`docs/technical_note.md`](docs/technical_note.md).

## Design

* Each replicate is one **community**: uptake preferences, leakage structure,
  baseline mortality and thermal optima are drawn once.
* The same community is simulated under **every** combination of leakage
  (0.0–0.8) and temperature (15–45 °C), so effects are estimated within
  communities.
* Random traits come from separate streams (uptake structure and scaling,
  leakage structure, mortality, thermal optima).
* `rng_mode="legacy"` reproduces the earlier single-script pilot
  (15 communities); a regression test compares against its stored results in
  `tests/data/`.

## Outcomes

* **Richness**: consumers with final biomass above 1e-5.
* **Biomass**: sum of final consumer biomass.
* **Gross retained fraction**: per consumer,
  `sum_j u_ij (1 - l_j) R0 / sum_j u_ij R0`. It is fixed by leakage and does not
  involve mortality.
* **Net-retention proxy**: gross retained fraction minus `m_i / sum_j u_ij R0`,
  where `m_i` is a generic first-order loss (maintenance, death, washout or a
  mixture). It is evaluated at a normalised reference environment (R0 = 1 for every
  resource, not the realised equilibrium) and is not a realised carbon-use
  efficiency.
* Both are averaged over survivors with **biomass weights** and with **flux
  weights** (biomass times uptake at the reference resource vector). Sensitivity to
  R0 = 0.5, 1 and 2 is computed from the stored final states.
* Convergence diagnostics are stored for every run: integrator success, largest
  residual of the equations, largest per-capita rate of change of survivors and
  the number of consumers near the extinction threshold.

## Install

Requires Python 3.11 or later.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m pip install -e .
```

`requirements.txt` pins the DigiMicPy commit used.

## Run

```bash
python -m pytest
python scripts/run_experiment.py --n-communities 15 --n-jobs 2 --out results/pilot
```

This writes `results/pilot.csv` (one row per simulation), the final states and
the exact configuration. The analysis and figures are reproduced with
`scripts/make_figures.py` and `scripts/reference_sensitivity.py`; numerical
checks with `scripts/run_sensitivity.py`.

## Layout

```
src/ltcue/    configuration, community construction, simulation, metrics, analysis
scripts/      run_experiment, merge_results, make_figures, reference_sensitivity, run_sensitivity
tests/        unit tests and pilot regression test
docs/         technical note and parameters
results/      summary tables (final states and sensitivity runs are not tracked)
figures/      figures of the technical note
```
