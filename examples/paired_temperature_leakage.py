"""Paired temperature x leakage comparison with convergence checks.

Each community is drawn once. Temperature and leakage only transform that
community, so every treatment combination of one community shares its uptake
structure, leakage structure, baseline mortality and thermal optima.

For each run the script reports richness, total biomass, the gross retained
fraction and a reference-condition net-retention proxy (biomass- and
flux-weighted) at several reference resource levels, together with two
convergence diagnostics: the largest absolute residual of the equations and the
largest per-capita rate of change of the survivors. Thermal parameters are
illustrative, not fitted to data. The proxy is not a realised community CUE.

Run from the repository root:
    python examples/paired_temperature_leakage.py
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from digimicpy import (
    MiCRMParameters,
    generate_l_tensor,
    micrm_rhs,
    modular_uptake,
    solve_micrm,
    temperature_adjusted_parameters,
)

N_CONSUMERS, N_RESOURCES = 20, 10
N_COMMUNITIES = 5
TEMPERATURES = (288.0, 303.0, 318.0)       # kelvin
LEAKAGE = (0.0, 0.4, 0.8)
REFERENCE_TEMPERATURE = 303.0
REFERENCE_LEVELS = (0.5, 1.0, 2.0)         # uniform reference resource levels
EXTINCTION = 1e-5
T_END = 60000.0                            # long: some survivors decay slowly


def draw_community(seed: int) -> dict:
    """Traits that stay fixed across all treatments of one community."""

    rng = np.random.default_rng(seed)
    uptake = modular_uptake(N_CONSUMERS, N_RESOURCES, 2, 8.0, rng=rng)
    uptake = uptake * rng.uniform(0.5, 1.5, (N_CONSUMERS, 1))
    # Total leakage 1: rows are scaled by the treatment value afterwards.
    leakage_unit = generate_l_tensor(N_CONSUMERS, N_RESOURCES, 2, 8.0, 1.0, rng=rng)
    return {
        "uptake": uptake,
        "leakage_unit": leakage_unit,
        "mortality": rng.uniform(0.15, 0.25, N_CONSUMERS),
        "uptake_topt": rng.uniform(298.0, 315.0, N_CONSUMERS),
    }


def parameters_at(community: dict, leakage: float, temperature: float) -> MiCRMParameters:
    reference = MiCRMParameters(
        uptake=community["uptake"],
        mortality=community["mortality"],
        resource_supply=np.full(N_RESOURCES, 1.0),
        resource_decay=np.full(N_RESOURCES, 0.5),
        leakage=community["leakage_unit"] * leakage,
        leakage_fraction=np.full(N_RESOURCES, leakage),
    )
    return temperature_adjusted_parameters(
        reference,
        temperature,
        REFERENCE_TEMPERATURE,
        uptake_activation_energy=0.65,
        uptake_optimum_temperature=community["uptake_topt"],
        uptake_deactivation_energy=3.0,
        mortality_activation_energy=0.30,
        mortality_optimum_temperature=340.0,
        mortality_deactivation_energy=3.0,
    )


def retention(parameters: MiCRMParameters, consumers, alive, r0: float):
    """Gross retained fraction and net-retention proxy per consumer, and weights."""

    uptake_ref = parameters.uptake.sum(axis=1) * r0
    retained_ref = (parameters.uptake * (1.0 - parameters.leakage.sum(axis=2))).sum(axis=1) * r0
    with np.errstate(divide="ignore", invalid="ignore"):
        gross = retained_ref / uptake_ref
        net = (retained_ref - parameters.mortality) / uptake_ref
    ok = alive & (uptake_ref > 0)
    if not ok.any():
        return np.nan, np.nan, np.nan, np.nan
    biomass_w, flux_w = consumers[ok], consumers[ok] * uptake_ref[ok]
    return (
        np.average(gross[ok], weights=biomass_w), np.average(gross[ok], weights=flux_w),
        np.average(net[ok], weights=biomass_w), np.average(net[ok], weights=flux_w),
    )


def run(seed: int, leakage: float, temperature: float) -> dict:
    parameters = parameters_at(draw_community(seed), leakage, temperature)
    x0 = np.concatenate([np.full(N_CONSUMERS, 0.01), np.ones(N_RESOURCES)])
    result = solve_micrm(parameters, x0, (0.0, T_END), method="LSODA", rtol=1e-8, atol=1e-12)
    state = result.y[:, -1]
    consumers = state[:N_CONSUMERS]
    alive = consumers > EXTINCTION
    rhs = micrm_rhs(0.0, state, parameters)
    out = {
        "success": bool(result.success),
        "richness": int(alive.sum()),
        "biomass": float(consumers.sum()),
        "resid_abs": float(np.max(np.abs(rhs))),
        # per-capita rate: an absolute residual can be tiny while survivors still decline
        "rate_rel": float(np.max(np.abs(rhs[:N_CONSUMERS][alive] / consumers[alive]))) if alive.any() else 0.0,
    }
    for r0 in REFERENCE_LEVELS:
        g_b, g_f, n_b, n_f = retention(parameters, consumers, alive, r0)
        out[f"gross_{r0}"], out[f"net_bio_{r0}"], out[f"net_flux_{r0}"] = g_b, n_b, n_f
    return out


def main() -> None:
    print(f"{N_COMMUNITIES} communities x {len(TEMPERATURES)} temperatures x {len(LEAKAGE)} leakage values")
    results = {(s, l, t): run(s, l, t) for s in range(N_COMMUNITIES) for l in LEAKAGE for t in TEMPERATURES}

    assert all(r["success"] for r in results.values()), "integrator failure"
    print(f"largest per-capita rate of change of survivors: {max(r['rate_rel'] for r in results.values()):.2e}")
    print(f"largest absolute residual: {max(r['resid_abs'] for r in results.values()):.2e}\n")

    print("mean over communities (same communities in every cell)")
    print(f"{'T (K)':>6} {'leak':>5} {'richness':>9} {'biomass':>8} {'gross':>6}  "
          + "  ".join(f"net(R0={r:g})" for r in REFERENCE_LEVELS))
    for t in TEMPERATURES:
        for l in LEAKAGE:
            cell = [results[(s, l, t)] for s in range(N_COMMUNITIES)]
            mean = lambda key: np.nanmean([c[key] for c in cell])
            print(f"{t:6.0f} {l:5.1f} {mean('richness'):9.2f} {mean('biomass'):8.2f} {mean('gross_1.0'):6.2f}  "
                  + "  ".join(f"{mean(f'net_bio_{r}'):11.3f}" for r in REFERENCE_LEVELS))

    # With a uniform R0 the proxy equals gross - K/R0 (K does not depend on R0).
    # Its sign changes at R0 = K/gross, and treatments with different gross and K
    # can swap order as R0 changes; only the order of treatments with the same
    # leakage is fixed. Count what happens in these runs, within each community.
    net = {key: np.array([r[f"net_bio_{r0}"] for r0 in REFERENCE_LEVELS]) for key, r in results.items()}
    sign_changes = sum(len(set(np.sign(v))) > 1 for v in net.values())
    pairs = {"same": 0, "different": 0}
    reversed_ = {"same": 0, "different": 0}
    for s in range(N_COMMUNITIES):
        cells = [(l, t) for l in LEAKAGE for t in TEMPERATURES]
        for a, b in combinations(cells, 2):
            kind = "same" if a[0] == b[0] else "different"
            pairs[kind] += 1
            reversed_[kind] += len(set(np.sign(net[(s, *a)] - net[(s, *b)]))) > 1
    print(f"\nnet proxy (biomass-weighted) across R0 = {REFERENCE_LEVELS}, within community:")
    print(f"  runs whose sign changes with R0: {sign_changes} of {len(net)}")
    print(f"  treatment pairs with the same leakage whose order changes: "
          f"{reversed_['same']} of {pairs['same']}")
    print(f"  treatment pairs with different leakage whose order changes: "
          f"{reversed_['different']} of {pairs['different']}")
    print("Rankings at equal leakage cannot depend on R0; signs and rankings across leakage "
          "levels can.\nReport the counts for the runs at hand (in the full grid of this "
          "repository some rankings\nacross leakage levels do reverse; see docs/technical_note.md).")


if __name__ == "__main__":
    main()
