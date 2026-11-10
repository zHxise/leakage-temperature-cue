import importlib.util
from pathlib import Path

import numpy as np

_PATH = Path(__file__).resolve().parents[1] / "examples" / "paired_temperature_leakage.py"
_SPEC = importlib.util.spec_from_file_location("paired_example", _PATH)
example = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(example)


def test_example_point_is_consistent():
    example.T_END = 3000.0  # short horizon; only structure is checked here
    row = example.run(seed=0, leakage=0.4, temperature=303.0)
    assert row["success"]
    # with uniform leakage the gross retained fraction equals 1 - leakage for any R0
    for r0 in example.REFERENCE_LEVELS:
        assert np.isclose(row[f"gross_{r0}"], 0.6)
    # net proxy = gross - K / R0: the loss term halves when R0 doubles
    loss = {r0: row[f"gross_{r0}"] - row[f"net_bio_{r0}"] for r0 in example.REFERENCE_LEVELS}
    assert np.isclose(loss[0.5], 2 * loss[1.0]) and np.isclose(loss[1.0], 2 * loss[2.0])


def test_slow_transient_example_loses_a_survivor_after_a_small_residual():
    path = _PATH.with_name("slow_transient.py")
    spec = importlib.util.spec_from_file_location("slow_transient", path)
    slow = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(slow)

    parameters = slow.community(slow.SEED)
    initial = np.concatenate([np.full(slow.N_CONSUMERS, 0.01), np.ones(slow.N_RESOURCES)])
    survivors = {}
    for horizon in (5000.0, 20000.0):
        result = example.solve_micrm(parameters, initial, (0.0, horizon), method=slow.METHOD,
                                     rtol=slow.RTOL, atol=slow.ATOL)
        state = result.y[:, -1]
        alive = state[: slow.N_CONSUMERS] > slow.SURVIVOR_THRESHOLD
        rhs = example.micrm_rhs(result.t[-1], state, parameters)
        survivors[horizon] = int(alive.sum())
        if horizon == 5000.0:
            # small absolute residual, yet one survivor is still declining
            assert np.max(np.abs(rhs)) < 1e-6
            assert np.max(np.abs(rhs[: slow.N_CONSUMERS][alive] / state[: slow.N_CONSUMERS][alive])) > 1e-3
    assert survivors == {5000.0: 5, 20000.0: 4}
