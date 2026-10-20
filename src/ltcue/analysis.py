"""Summaries and paired comparisons.

Replicates are communities, so uncertainty is quantified by bootstrapping over
communities (the unit that is independently sampled), never over rows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def bootstrap_mean_ci(
    values: np.ndarray,
    n_boot: int = 4000,
    level: float = 0.95,
    seed: int = 0,
) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    if values.size == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(n_boot, values.size), replace=True).mean(axis=1)
    alpha = (1.0 - level) / 2.0
    return float(values.mean()), float(np.quantile(draws, alpha)), float(np.quantile(draws, 1 - alpha))


def summarize(table: pd.DataFrame, metrics: list[str], by: list[str]) -> pd.DataFrame:
    """Mean and bootstrap CI of each metric within each group."""

    rows = []
    for number, (key, group) in enumerate(table.groupby(by)):
        key = key if isinstance(key, tuple) else (key,)
        for metric in metrics:
            # a distinct seed per group keeps resampling indices independent across cells
            mean, lo, hi = bootstrap_mean_ci(group[metric].to_numpy(), seed=number)
            rows.append({**dict(zip(by, key)), "metric": metric, "mean": mean, "lo": lo, "hi": hi,
                         "n": int(group[metric].notna().sum())})
    return pd.DataFrame(rows)


def leakage_slopes(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Per community and temperature: least-squares slope of ``metric`` on leakage.

    Each slope is computed within one community, so it is a paired effect of
    leakage that is free of between-community variation.
    """

    rows = []
    for (rep, temp), group in table.groupby(["rep", "T"]):
        valid = group.dropna(subset=[metric])
        if valid["leak"].nunique() < 2:
            continue
        slope = np.polyfit(valid["leak"], valid[metric], 1)[0]
        rows.append({"rep": rep, "T": temp, "slope": slope})
    return pd.DataFrame(rows)


def slope_summary(slopes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for number, (temp, group) in enumerate(slopes.groupby("T")):
        mean, lo, hi = bootstrap_mean_ci(group["slope"].to_numpy(), seed=number)
        rows.append({"T": temp, "mean": mean, "lo": lo, "hi": hi, "n": len(group)})
    return pd.DataFrame(rows)


def slope_contrast(slopes: pd.DataFrame, temp_a: float, temp_b: float) -> tuple[float, float, float]:
    """Paired difference of leakage slopes (temp_b - temp_a), across communities."""

    wide = slopes.pivot(index="rep", columns="T", values="slope")
    return bootstrap_mean_ci((wide[temp_b] - wide[temp_a]).to_numpy())
