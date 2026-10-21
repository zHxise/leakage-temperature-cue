"""Concatenate chunked runs written by run_experiment.py --reps.

Example:
    python scripts/merge_results.py results/chunks/main_*.csv --out results/main
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chunks", nargs="+", type=Path, help="chunk .csv files")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n-communities", type=int, required=True)
    args = parser.parse_args()

    chunks = sorted(args.chunks, key=lambda p: pd.read_csv(p)["rep"].min())
    table = pd.concat([pd.read_csv(p) for p in chunks], ignore_index=True)
    states = np.vstack([np.load(p.with_suffix(".npz"))["states"] for p in chunks])
    assert len(table) == len(states)
    assert not table.duplicated(["leak", "T", "rep"]).any(), "duplicated simulations"
    assert sorted(table["rep"].unique()) == list(range(args.n_communities)), "non-contiguous replicates"
    assert table["rep"].nunique() == args.n_communities, "missing replicates"

    table.to_csv(args.out.with_suffix(".csv"), index=False)
    np.savez_compressed(args.out.with_suffix(".npz"), states=states)
    config = chunks[0].with_suffix(".config.json").read_text()
    args.out.with_suffix(".config.json").write_text(config)
    print(f"merged {len(chunks)} chunks, {len(table)} simulations")


if __name__ == "__main__":
    main()
