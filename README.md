# ltcue: temperature × leakage in microbial consumer-resource communities

Simulations of how temperature and metabolite leakage jointly affect
**persistence** (number of surviving consumers), **total biomass** and
**community carbon-use efficiency (CUE)**, built on
[DigiMicPy](https://github.com/DigiMicOrg/DigiMicPy).

Thermal parameters are illustrative placeholders; see
[`docs/parameters.md`](docs/parameters.md) for which processes depend on
temperature and with which values.

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
* **CUE**: per consumer, `(sum_j u_ij (1 - l_j) R0 - m_i) / (sum_j u_ij R0)` with
  R0 = 1 (mortality treated as maintenance loss), averaged over survivors with
  biomass weights.
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
the exact configuration.

## Layout

```
src/ltcue/    configuration, community construction, simulation, metrics, analysis
scripts/      run_experiment.py
tests/        unit tests and pilot regression test
docs/         thermal and model parameters
```
