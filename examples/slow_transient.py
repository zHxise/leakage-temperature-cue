"""Slow transient: a small endpoint residual does not show that survivors have settled.

One community of the temperature x leakage comparison (seed 4, no leakage, at
the reference temperature of 303 K, where the thermal scaling factor is 1, so
the illustrative thermal parameters play no role here). The same initial
condition is integrated to four horizons. At t = 5000 the largest absolute
residual of the equations is already small, but one consumer above the survivor
threshold is still declining at a per-capita rate near 2e-3 and is lost at the
longer horizons.

Needs DigiMicPy and NumPy only. Run from the repository root:
    python examples/slow_transient.py
"""

from __future__ import annotations

import json
from importlib.metadata import distribution

import numpy as np
import scipy
from digimicpy import MiCRMParameters, generate_l_tensor, micrm_rhs, modular_uptake, solve_micrm

SEED = 4
N_CONSUMERS, N_RESOURCES = 20, 10
HORIZONS = (5000.0, 20000.0, 60000.0, 120000.0)
SURVIVOR_THRESHOLD = 1e-5
METHOD, RTOL, ATOL = "LSODA", 1e-8, 1e-12


def community(seed: int) -> MiCRMParameters:
    """Same draws, in the same order, as examples/paired_temperature_leakage.py."""

    rng = np.random.default_rng(seed)
    uptake = modular_uptake(N_CONSUMERS, N_RESOURCES, 2, 8.0, rng=rng)
    uptake = uptake * rng.uniform(0.5, 1.5, (N_CONSUMERS, 1))
    leakage_unit = generate_l_tensor(N_CONSUMERS, N_RESOURCES, 2, 8.0, 1.0, rng=rng)
    return MiCRMParameters(
        uptake=uptake,
        mortality=rng.uniform(0.15, 0.25, N_CONSUMERS),
        resource_supply=np.full(N_RESOURCES, 1.0),
        resource_decay=np.full(N_RESOURCES, 0.5),
        leakage=leakage_unit * 0.0,
        leakage_fraction=np.zeros(N_RESOURCES),
    )


def digimicpy_version() -> str:
    dist = distribution("digimicpy")
    try:
        commit = json.loads(dist.read_text("direct_url.json"))["vcs_info"]["commit_id"][:7]
    except Exception:
        commit = "unknown"
    return f"{dist.version} (commit {commit})"


def main() -> None:
    parameters = community(SEED)
    initial = np.concatenate([np.full(N_CONSUMERS, 0.01), np.ones(N_RESOURCES)])

    print(f"DigiMicPy {digimicpy_version()}, NumPy {np.__version__}, SciPy {scipy.__version__}")
    print(f"seed {SEED} (numpy.random.default_rng), {N_CONSUMERS} consumers, {N_RESOURCES} resources")
    print(f"solver {METHOD}, rtol {RTOL:g}, atol {ATOL:g}, survivor threshold {SURVIVOR_THRESHOLD:g}")
    print("initial state: consumers 0.01, resources 1.0\n")
    print(f"{'t_end':>7} {'success':>7} {'survivors':>9} {'max |residual|':>15} "
          f"{'max |dC/dt|/C':>14} {'smallest survivor':>18}")

    for horizon in HORIZONS:
        result = solve_micrm(parameters, initial, (0.0, horizon), method=METHOD, rtol=RTOL, atol=ATOL)
        state = result.y[:, -1]
        consumers = state[:N_CONSUMERS]
        alive = consumers > SURVIVOR_THRESHOLD
        rhs = micrm_rhs(result.t[-1], state, parameters)
        per_capita = np.abs(rhs[:N_CONSUMERS][alive] / consumers[alive])
        print(f"{horizon:7.0f} {str(bool(result.success)):>7} {int(alive.sum()):9d} "
              f"{np.max(np.abs(rhs)):15.1e} {per_capita.max():14.1e} {consumers[alive].min():18.1e}")


if __name__ == "__main__":
    main()
