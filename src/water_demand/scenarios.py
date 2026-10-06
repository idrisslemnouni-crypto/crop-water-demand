"""Apply an explicitly hypothetical root zone to previously calculated real weather."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from water_demand.balance import RootZone, daily_balance


def run(root: Path):
    path = root / "reports/daily-water-demand.csv"
    frame = pd.read_csv(path, parse_dates=["date"])
    config = json.loads((root / "configs/root-zone-scenario.json").read_text(encoding="utf-8"))
    selected = frame.loc[frame.date.between(config["start_date"], config["end_date"])].copy()
    calendar = pd.date_range(config["start_date"], config["end_date"], freq="D")
    zone = RootZone(
        **{
            key: config[key]
            for key in ("theta_fc", "theta_wp", "root_depth_m", "depletion_fraction")
        }
    )
    rows = []
    for site, group in selected.groupby("site"):
        group = group.sort_values("date")
        if not pd.DatetimeIndex(group.date).equals(calendar):
            raise ValueError("A complete unique scenario calendar is required for every site")
        depletion = config["initial_depletion_mm"]
        for row in group.itertuples(index=False):
            balance = daily_balance(
                zone,
                depletion,
                float(row.et0_pm_mm),
                config["kc"],
                float(row.PRECTOTCORR),
                config["runoff_mm"],
                config["net_irrigation_mm"],
                config["capillary_rise_mm"],
            )
            rows.append(
                dict(
                    site=site,
                    date=row.date,
                    et0_mm=row.et0_pm_mm,
                    precipitation_mm=row.PRECTOTCORR,
                    **balance,
                )
            )
            depletion = balance["depletion_end_mm"]
    if not rows:
        raise ValueError("No scenario observations")
    output = pd.DataFrame(rows)
    summaries = {}
    for site, group in output.groupby("site"):
        summaries[site] = {
            "days": len(group),
            "precipitation_mm": float(group.precipitation_mm.sum()),
            "potential_etc_mm": float(group.potential_etc_mm.sum()),
            "adjusted_etc_scenario_mm": float(group.adjusted_etc_scenario_mm.sum()),
            "deep_percolation_mm": float(group.deep_percolation_mm.sum()),
            "days_with_start_ks_below_one": int(group.ks_start.lt(1).sum()),
            "end_depletion_mm": float(group.depletion_end_mm.iloc[-1]),
            "max_absolute_balance_residual_mm": float(group.balance_residual_mm.abs().max()),
        }
    if (
        not np.isfinite(output.select_dtypes("number")).all().all()
        or output.balance_residual_mm.abs().max() > 1e-8
    ):
        raise ValueError("Invalid flux or unclosed scenario water budget")
    result = dict(
        config=config,
        source_daily_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        scope="Hypothetical rainfed fixed-root-zone response to real gridded meteorology; no field water-stress or irrigation validation",
        sites=summaries,
    )
    reports = root / "reports"
    output.to_csv(reports / "root-zone-scenarios.csv", index=False)
    (reports / "root-zone-scenarios.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    fig = Figure(figsize=(11, 6))
    axes = fig.subplots(2, 1, sharex=True)
    for site, group in output.groupby("site"):
        axes[0].plot(group.date, group.depletion_end_mm, label=site)
        axes[1].plot(group.date, group.ks_start, label=site)
    axes[0].axhline(zone.raw_mm, color="black", linestyle="--", label="Assumed RAW")
    axes[0].set(
        ylabel="Root-zone depletion (mm)",
        ylim=(0, zone.taw_mm + 10),
        title="Hypothetical fixed soil/crop scenario · real 2024 gridded weather",
    )
    axes[1].set(ylabel="Scenario Ks at day start", xlabel="Date", ylim=(-0.03, 1.03))
    axes[0].legend(ncol=4)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(reports / "figures/root-zone-scenarios.png", dpi=150)
    return result


if __name__ == "__main__":
    print(json.dumps(run(Path.cwd()), indent=2))
