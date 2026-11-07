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
