"""Sensitivity of the retention metrics to the reference resource level R0.

Recomputes the net-retention proxy and the gross retained fraction at
R0 in {0.5, 1, 2} from the stored final states of a run (no new simulations),
with biomass weights and with flux weights (biomass x uptake at the reference
resource vector). For a uniform R0 the proxy is exactly G - K/R0, so signs and
orderings do not depend on R0 and only magnitudes (and the share of negative
values) scale; the table makes this explicit.

Example:
    python scripts/reference_sensitivity.py --prefix results/main
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ltcue import ExperimentConfig, build_community, parameters_at
from ltcue.analysis import bootstrap_mean_ci, leakage_slopes, summarize
from ltcue.metrics import (
    community_average,
    gross_retained_fraction,
    net_retention_proxy,
    survivors,
    uptake_at_reference,
)

REFERENCE_LEVELS = (0.5, 1.0, 2.0)


def per_run_table(cfg, table, states) -> pd.DataFrame:
    communities = {rep: build_community(cfg, rep) for rep in table["rep"].unique()}
    rows = []
    for row, state in zip(table.itertuples(index=False), states):
        p = parameters_at(communities[row.rep], cfg, row.leak, row.T)
        consumers = state[: cfg.n_consumers]
        alive = survivors(consumers, cfg.extinction_threshold)
        for r0 in REFERENCE_LEVELS:
            proxy = net_retention_proxy(p, r0)
            gross = gross_retained_fraction(p, r0)
            flux = consumers * uptake_at_reference(p, r0)
            rows.append({
                "R0": r0, "leak": row.leak, "T": row.T, "rep": row.rep,
                "proxy": community_average(proxy, consumers, alive),
                "proxy_flux": community_average(proxy, consumers, alive, flux),
                "gross": community_average(gross, consumers, alive),
                "gross_flux": community_average(gross, consumers, alive, flux),
            })
    return pd.DataFrame(rows)


def paired_difference(frame: pd.DataFrame, value: str, leak: float, t_a: float, t_b: float):
    wide = frame[frame["leak"] == leak].pivot(index="rep", columns="T", values=value)
    return bootstrap_mean_ci((wide[t_b] - wide[t_a]).to_numpy())


def claims_table(frame: pd.DataFrame, cfg) -> pd.DataFrame:
    low, mid, hot_start, hot = 288.0, 300.0, 312.0, 318.0
    rows = []
    for r0, sub in frame.groupby("R0"):
        for value in ("proxy", "proxy_flux"):
            means = sub.groupby(["leak", "T"])[value].mean().unstack()
            slopes = leakage_slopes(sub.rename(columns={value: "y"}), "y")
            slope_ci = [bootstrap_mean_ci(g["slope"].to_numpy()) for _, g in slopes.groupby("T")]
            drop = paired_difference(sub, value, 0.0, hot_start, hot)
            cold = paired_difference(sub, value, 0.0, mid, low)
            negative = sub[(sub["leak"] == max(cfg.leakage_fractions))]
            neg_low = negative[negative["T"] == low][value]
            neg_hot = negative[negative["T"] == hot][value]
            # 27 C against both neighbours, at every leakage level (paired by community)
            peak_ok = all(
                paired_difference(sub, value, lk, t_nb, mid)[1] > 0
                for lk in cfg.leakage_fractions for t_nb in (294.0, 306.0)
            )
            rows.append({
                "R0": r0, "metric": value,
                "T_of_max_no_leakage_K": float(means.loc[0.0].idxmax()),
                "leakage_slope_all_T_negative": bool(all(hi < 0 for _, _, hi in slope_ci)),
                "change_39_to_45C_no_leakage": drop[0], "lo_39_45": drop[1], "hi_39_45": drop[2],
                "change_27_to_15C_no_leakage": cold[0], "lo_27_15": cold[1], "hi_27_15": cold[2],
                "mid_higher_than_neighbours_all_leakage": bool(peak_ok),
                "n_negative_at_max_leakage_15C": int((neg_low < 0).sum()),
                "n_with_survivors_at_max_leakage_15C": int(neg_low.notna().sum()),
                "share_negative_at_max_leakage_15C": float((neg_low < 0).sum() / neg_low.notna().sum()),
                "n_negative_at_max_leakage_45C": int((neg_hot < 0).sum()),
                "n_with_survivors_at_max_leakage_45C": int(neg_hot.notna().sum()),
                "share_negative_at_max_leakage_45C": float((neg_hot < 0).sum() / neg_hot.notna().sum()),
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=Path("results/main"))
    args = parser.parse_args()

    table = pd.read_csv(args.prefix.with_suffix(".csv"))
    states = np.load(args.prefix.with_suffix(".npz"))["states"]
    cfg = ExperimentConfig.from_json(args.prefix.with_suffix(".config.json").read_text())

    frame = per_run_table(cfg, table, states)
    out = args.prefix.parent / args.prefix.name
    summary = summarize(frame, ["proxy", "proxy_flux", "gross", "gross_flux"], ["R0", "leak", "T"])
    summary.to_csv(f"{out}_reference_summary.csv", index=False)
    claims = claims_table(frame, cfg)
    claims.to_csv(f"{out}_reference_claims.csv", index=False)

    for value in ("proxy", "proxy_flux"):
        for r0 in REFERENCE_LEVELS:
            sub = summary[(summary["metric"] == value) & (summary["R0"] == r0)]
            wide = sub.pivot(index="leak", columns="T", values="mean")
            print(f"\n{value}, R0 = {r0}\n{wide.round(3).to_string()}")
    print("\n" + claims.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
