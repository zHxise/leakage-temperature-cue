"""Numerical sensitivity of the headline results.

Re-runs the same communities under alternative solver settings, integration
horizons and initial conditions, and recomputes outcomes under alternative
extinction thresholds from the stored baseline final states.

Each variant is run and stored separately so the work can be done in pieces:

    python scripts/run_sensitivity.py run baseline --n-communities 20 --n-jobs 2 --reps 0 10
    python scripts/run_sensitivity.py run baseline --n-communities 20 --n-jobs 2 --reps 10 20
    python scripts/run_sensitivity.py run tight_tol --n-communities 20 --n-jobs 2 --reps 0 10
    ...
    python scripts/run_sensitivity.py summarize --n-communities 20

Variants: baseline, tight_tol, radau, horizon_5000, horizon_20000, horizon_120000,
init_low, init_high. Omitting --reps stores a single <variant>.csv instead of chunks.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ltcue import ExperimentConfig, build_community, parameters_at, run_experiment
from ltcue.analysis import leakage_slopes
from ltcue.metrics import (
    community_average,
    net_retention_proxy,
    survivors,
    uptake_at_reference,
)

KEYS = ["leak", "T", "rep"]


def load_variant(directory: Path, key: str) -> tuple[pd.DataFrame, np.ndarray] | None:
    """Load a stored variant: either ``key.csv`` or chunks ``key__*.csv``, never both."""

    whole = sorted(directory.glob(f"{key}.csv"))
    chunks = sorted(directory.glob(f"{key}__*.csv"))
    if whole and chunks:
        raise RuntimeError(f"both {key}.csv and {key}__*.csv exist; remove one set")
    paths = whole or chunks
    if not paths:
        return None
    table = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    states = np.vstack([np.load(p.with_suffix(".npz"))["states"] for p in paths])
    if table.duplicated(KEYS).any():
        raise RuntimeError(f"duplicated simulations in variant {key}")
    return table, states


def compare(base: pd.DataFrame, other: pd.DataFrame, name: str) -> dict:
    merged = base.merge(other, on=KEYS, suffixes=("", "_alt"))
    d_rich = (merged["richness_alt"] - merged["richness"]).abs()
    rel_bio = ((merged["biomass_alt"] - merged["biomass"]).abs() / merged["biomass"].clip(lower=1e-12))
    d_proxy = (merged["retention_proxy_alt"] - merged["retention_proxy"]).abs()

    summary = {
        "variant": name,
        "n": len(merged),
        "integrator_failures": int((~other["ok"]).sum()),
        "share_richness_changed": float((d_rich > 0).mean()),
        "max_richness_change": int(d_rich.max()),
        "median_rel_biomass_change": float(rel_bio.median()),
        "max_rel_biomass_change": float(rel_bio.max()),
        "max_proxy_change": float(np.nanmax(d_proxy)),
        "max_rate_rel": float(other["rate_rel"].max()),
    }
    for metric in ("richness", "biomass", "retention_proxy"):
        s_base = leakage_slopes(base, metric).groupby("T")["slope"].mean()
        s_alt = leakage_slopes(other, metric).groupby("T")["slope"].mean()
        summary[f"{metric}_slope_sign_kept"] = bool((np.sign(s_base) == np.sign(s_alt)).all())
    return summary


def threshold_variant(cfg, table, states, threshold):
    """Recompute outcomes from stored final states under another extinction threshold."""

    out = table.copy()
    richness, proxy, proxy_flux = [], [], []
    for row, state in zip(table.itertuples(index=False), states):
        parameters = parameters_at(build_community(cfg, row.rep), cfg, row.leak, row.T)
        consumers = state[: cfg.n_consumers]
        alive = survivors(consumers, threshold)
        per_consumer = net_retention_proxy(parameters, cfg.reference_resource)
        flux = consumers * uptake_at_reference(parameters, cfg.reference_resource)
        richness.append(int(alive.sum()))
        proxy.append(community_average(per_consumer, consumers, alive))
        proxy_flux.append(community_average(per_consumer, consumers, alive, flux))
    out["richness"], out["retention_proxy"], out["retention_proxy_flux"] = richness, proxy, proxy_flux
    return out


VARIANTS = {
    "baseline": {},
    "tight_tol": dict(rtol=1e-10, atol=1e-14),
    "radau": dict(method="Radau"),
    "horizon_5000": dict(t_end=5000.0),
    "horizon_20000": dict(t_end=20000.0),
    "horizon_120000": dict(t_end=120000.0),
    "init_low": dict(initial_consumer_biomass=0.001),
    "init_high": dict(initial_consumer_biomass=0.1, initial_resource=3.0),
}
LABELS = {
    "tight_tol": "rtol 1e-10, atol 1e-14",
    "radau": "Radau solver",
    "horizon_5000": "horizon 5000 (pilot)",
    "horizon_20000": "horizon 20000",
    "horizon_120000": "horizon 120000",
    "init_low": "initial biomass 0.001",
    "init_high": "initial biomass 0.1, resource 3",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("variant", choices=sorted(VARIANTS))
    run.add_argument("--n-jobs", type=int, default=1)
    run.add_argument("--reps", type=int, nargs=2, metavar=("START", "STOP"), default=None,
                     help="run only these community indices; saved as <variant>__<START>")
    summ = sub.add_parser("summarize")
    for p in (run, summ):
        p.add_argument("--n-communities", type=int, default=20)
        p.add_argument("--dir", type=Path, default=Path("results/sensitivity"))
    args = parser.parse_args()

    cfg = ExperimentConfig(n_communities=args.n_communities)
    args.dir.mkdir(parents=True, exist_ok=True)

    if args.command == "run":
        reps = None if args.reps is None else range(*args.reps)
        table, states = run_experiment(
            cfg.with_(**VARIANTS[args.variant]), n_jobs=args.n_jobs, replicates=reps
        )
        name = args.variant if reps is None else f"{args.variant}__{args.reps[0]:03d}"
        table.to_csv(args.dir / f"{name}.csv", index=False)
        np.savez_compressed(args.dir / f"{name}.npz", states=states)
        print(f"{name}: {len(table)} runs, failures {int((~table['ok']).sum())}")
        return

    loaded = load_variant(args.dir, "baseline")
    if loaded is None:
        raise SystemExit("run the baseline variant first")
    base, states = loaded
    summaries = []
    for key, label in LABELS.items():
        loaded = load_variant(args.dir, key)
        if loaded is not None:
            table = loaded[0]
            if table["rep"].nunique() == args.n_communities:
                summaries.append(compare(base, table, label))
            else:
                print(f"skipping {label}: only {table['rep'].nunique()} communities")
    for threshold in (1e-8, 1e-6, 1e-4, 1e-3):
        table = threshold_variant(cfg, base, states, threshold)
        summaries.append(compare(base, table, f"extinction threshold {threshold:g}"))

    result = pd.DataFrame(summaries)
    result.to_csv(args.dir.parent / "sensitivity_summary.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
