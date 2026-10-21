"""Figures and summary tables from a stored experiment.

Example:
    python scripts/make_figures.py --prefix results/main --figdir figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from digimicpy import thermal_scaling_factor

from ltcue import ExperimentConfig, build_community
from ltcue.analysis import leakage_slopes, slope_contrast, slope_summary, summarize

INK = "#1f2933"
MUTED = "#6b7280"
GRID = "#e5e7eb"
MORTALITY = "#c2410c"
UPTAKE = "#1d4ed8"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "font.size": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False,
})


def leak_colors(values):
    ramp = plt.get_cmap("Blues")
    return {v: ramp(0.3 + 0.65 * i / max(len(values) - 1, 1)) for i, v in enumerate(values)}


def celsius(t):
    return np.asarray(t) - 273.15


def thermal_figure(cfg, path):
    community = build_community(cfg, 0)
    th = cfg.thermal
    temps = np.linspace(cfg.temperatures_k[0] - 3, cfg.temperatures_k[-1] + 3, 200)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for i, topt in enumerate(community.uptake_topt):
        factor = [thermal_scaling_factor(t, th.uptake_ea, topt, th.uptake_ed, cfg.reference_temperature_k)
                  for t in temps]
        ax.plot(celsius(temps), factor, color=UPTAKE, alpha=0.35, lw=1,
                label="Uptake (one line per consumer)" if i == 0 else None)
    mort = [thermal_scaling_factor(t, th.mortality_ea, th.mortality_topt, th.mortality_ed,
                                   cfg.reference_temperature_k) for t in temps]
    ax.plot(celsius(temps), mort, color=MORTALITY, lw=2, label="Mortality (shared)")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(f"Rate relative to {cfg.reference_temperature_k - 273.15:.0f} °C")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def lines_by_leakage(summary, metrics, titles, ylabels, path, cfg):
    colors = leak_colors(cfg.leakage_fractions)
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.4 * len(metrics), 3.4), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, metric, title, ylabel in zip(axes, metrics, titles, ylabels):
        for leak in cfg.leakage_fractions:
            sub = summary[(summary["metric"] == metric) & (summary["leak"] == leak)].sort_values("T")
            x = celsius(sub["T"])
            ax.fill_between(x, sub["lo"], sub["hi"], color=colors[leak], alpha=0.18, lw=0)
            ax.plot(x, sub["mean"], color=colors[leak], lw=2, marker="o", ms=4, label=f"{leak:.1f}")
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(ylabel)
    axes[-1].legend(title="Leakage", loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def reference_figure(reference, cfg, path):
    """Net-retention proxy against temperature for three reference resource levels."""

    levels = sorted(reference["R0"].unique())
    shades = {r0: plt.get_cmap("Blues")(0.35 + 0.6 * i / max(len(levels) - 1, 1))
              for i, r0 in enumerate(levels)}
    leaks = (min(cfg.leakage_fractions), max(cfg.leakage_fractions))
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.4), sharey=True)
    for ax, leak in zip(axes, leaks):
        for r0 in levels:
            sub = reference[(reference["metric"] == "proxy") & (reference["R0"] == r0)
                            & (reference["leak"] == leak)].sort_values("T")
            ax.fill_between(celsius(sub["T"]), sub["lo"], sub["hi"], color=shades[r0], alpha=0.18, lw=0)
            ax.plot(celsius(sub["T"]), sub["mean"], color=shades[r0], lw=2, marker="o", ms=4,
                    label=f"{r0:g}")
        ax.axhline(0, color=MUTED, lw=1)
        ax.set_title(f"Leakage {leak:g}", loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Temperature (°C)")
    axes[0].set_ylabel("Net-retention proxy (biomass-weighted)")
    axes[-1].legend(title="Reference R0", loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def slopes_figure(slope_tables, titles, ylabels, path):
    fig, axes = plt.subplots(1, len(slope_tables), figsize=(4.4 * len(slope_tables), 3.4), sharex=True)
    for ax, table, title, ylabel in zip(np.atleast_1d(axes), slope_tables, titles, ylabels):
        x = celsius(table["T"])
        ax.axhline(0, color=MUTED, lw=1)
        ax.fill_between(x, table["lo"], table["hi"], color=UPTAKE, alpha=0.2, lw=0)
        ax.plot(x, table["mean"], color=UPTAKE, lw=2, marker="o", ms=4)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=Path("results/main"))
    parser.add_argument("--figdir", type=Path, default=Path("figures"))
    args = parser.parse_args()

    table = pd.read_csv(args.prefix.with_suffix(".csv"))
    cfg = ExperimentConfig.from_json(args.prefix.with_suffix(".config.json").read_text())
    args.figdir.mkdir(parents=True, exist_ok=True)

    metrics = ["richness", "biomass", "retention_proxy", "retention_proxy_flux",
               "gross_retained", "gross_retained_flux"]
    summary = summarize(table, metrics, ["leak", "T"])
    summary.to_csv(args.prefix.parent / f"{args.prefix.name}_summary.csv", index=False)

    thermal_figure(cfg, args.figdir / "fig1_thermal_traits.png")
    lines_by_leakage(summary, ["richness", "biomass"], ["Persistence", "Total consumer biomass"],
                     ["Surviving consumers (of %d)" % cfg.n_consumers, "Biomass (arbitrary units)"],
                     args.figdir / "fig2_persistence_biomass.png", cfg)
    lines_by_leakage(
        summary, ["retention_proxy", "retention_proxy_flux"],
        ["Net-retention proxy, biomass-weighted", "Net-retention proxy, flux-weighted"],
        [f"Proxy at R0 = {cfg.reference_resource:g}"] * 2,
        args.figdir / "fig3_net_retention_proxy.png", cfg,
    )
    reference = args.prefix.parent / f"{args.prefix.name}_reference_summary.csv"
    if reference.exists():
        reference_figure(pd.read_csv(reference), cfg, args.figdir / "fig5_reference_level.png")

    slope_tables, rows, contrasts = [], [], []
    for metric in ("richness", "biomass", "retention_proxy"):
        slopes = leakage_slopes(table, metric)
        s = slope_summary(slopes)
        s.insert(0, "metric", metric)
        rows.append(s)
        if metric != "retention_proxy":
            slope_tables.append(s)
        low, high = min(cfg.temperatures_k), max(cfg.temperatures_k)
        mean, lo, hi = slope_contrast(slopes, low, high)
        contrasts.append({"metric": metric, "T_low": low, "T_high": high,
                          "contrast": mean, "lo": lo, "hi": hi})
        print(f"{metric}: slope(T={high:.0f}) - slope(T={low:.0f}) = {mean:.3f} [{lo:.3f}, {hi:.3f}]")
    pd.concat(rows).to_csv(args.prefix.parent / f"{args.prefix.name}_leakage_slopes.csv", index=False)
    pd.DataFrame(contrasts).to_csv(
        args.prefix.parent / f"{args.prefix.name}_slope_contrasts.csv", index=False
    )
    slopes_figure(slope_tables, ["Effect of leakage on persistence", "Effect of leakage on biomass"],
                  ["Change in richness per unit leakage", "Change in biomass per unit leakage"],
                  args.figdir / "fig4_leakage_effect_by_temperature.png")


if __name__ == "__main__":
    main()
