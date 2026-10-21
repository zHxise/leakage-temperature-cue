"""Run one treatment combination, or the full factorial experiment."""

from __future__ import annotations

from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from digimicpy import solve_micrm

from .communities import Community, build_community, parameters_at
from .config import ExperimentConfig
from .metrics import (
    community_average,
    convergence_diagnostics,
    gross_retained_fraction,
    net_retention_proxy,
    survivors,
    uptake_at_reference,
)


def simulate_point(
    cfg: ExperimentConfig,
    community: Community,
    leakage: float,
    temperature_k: float,
) -> tuple[dict, np.ndarray]:
    """Integrate one community under one treatment; return a result row and final state."""

    n, m = cfg.n_consumers, cfg.n_resources
    parameters = parameters_at(community, cfg, leakage, temperature_k)
    x0 = np.concatenate(
        [np.full(n, cfg.initial_consumer_biomass), np.full(m, cfg.initial_resource)]
    )
    result = solve_micrm(
        parameters,
        x0,
        (0.0, cfg.t_end),
        method=cfg.method,
        rtol=cfg.rtol,
        atol=cfg.atol,
    )
    state = result.y[:, -1]
    consumers = state[:n]
    alive = survivors(consumers, cfg.extinction_threshold)
    r0 = cfg.reference_resource
    proxy = net_retention_proxy(parameters, r0)
    gross = gross_retained_fraction(parameters, r0)
    flux_weight = consumers * uptake_at_reference(parameters, r0)
    threshold = cfg.extinction_threshold

    row = {
        "leak": leakage,
        "T": temperature_k,
        "rep": community.replicate,
        "ok": bool(result.success),
        "richness": int(alive.sum()),
        "biomass": float(consumers.sum()),
        "retention_proxy": community_average(proxy, consumers, alive),
        "retention_proxy_flux": community_average(proxy, consumers, alive, flux_weight),
        "gross_retained": community_average(gross, consumers, alive),
        "gross_retained_flux": community_average(gross, consumers, alive, flux_weight),
        "n_negative_proxy": int((proxy[alive] < 0).sum()),
        "near_threshold": int(
            ((consumers > 0.1 * threshold) & (consumers < 10.0 * threshold)).sum()
        ),
        "nfev": int(result.nfev),
        **convergence_diagnostics(state, parameters, alive),
    }
    return row, state


def _run_replicate(args: tuple[ExperimentConfig, int]) -> tuple[list[dict], list[np.ndarray]]:
    cfg, replicate = args
    community = build_community(cfg, replicate)
    rows, states = [], []
    for leakage in cfg.leakage_fractions:
        for temperature in cfg.temperatures_k:
            row, state = simulate_point(cfg, community, leakage, temperature)
            rows.append(row)
            states.append(state)
    return rows, states


def run_experiment(
    cfg: ExperimentConfig,
    n_jobs: int = 1,
    replicates: Iterable[int] | None = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Run every (leakage, temperature) combination for every replicate.

    Returns a tidy table (one row per simulation) and the matching array of
    final states, ordered identically to the table. ``replicates`` selects a
    subset of community indices (default: all ``cfg.n_communities``); results of
    disjoint subsets can be concatenated.
    """

    indices = range(cfg.n_communities) if replicates is None else replicates
    jobs = [(cfg, rep) for rep in indices]
    if n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            chunks = list(pool.map(_run_replicate, jobs))
    else:
        chunks = [_run_replicate(job) for job in jobs]

    rows = [row for chunk_rows, _ in chunks for row in chunk_rows]
    states = np.vstack([s for _, chunk_states in chunks for s in chunk_states])
    table = pd.DataFrame(rows)
    return table, states
