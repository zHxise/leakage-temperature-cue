"""Run the temperature x leakage factorial experiment and store the results.

Example:
    python scripts/run_experiment.py --n-communities 60 --n-jobs 2 --out results/main
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from ltcue import ExperimentConfig, run_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-communities", type=int, default=15)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--rng-mode", choices=["independent", "legacy"], default="independent")
    parser.add_argument("--reps", type=int, nargs=2, metavar=("START", "STOP"), default=None,
                        help="only run community indices START..STOP-1 (runs a subset; outputs are written per prefix)")
    parser.add_argument("--out", type=Path, default=Path("results/main"),
                        help="output prefix (writes <prefix>.csv, .npz and .config.json)")
    args = parser.parse_args()

    cfg = ExperimentConfig(n_communities=args.n_communities, rng_mode=args.rng_mode)
    start = time.time()
    table, states = run_experiment(
        cfg, n_jobs=args.n_jobs,
        replicates=None if args.reps is None else range(*args.reps),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out.with_suffix(".csv"), index=False)
    np.savez_compressed(args.out.with_suffix(".npz"), states=states)
    args.out.with_suffix(".config.json").write_text(cfg.to_json())

    print(f"{len(table)} simulations in {time.time() - start:.0f} s")
    print(f"integrator failures: {int((~table['ok']).sum())}")
    print(f"max |rhs|: {table['resid_abs'].max():.2e}; max per-capita rate: {table['rate_rel'].max():.2e}")


if __name__ == "__main__":
    main()
