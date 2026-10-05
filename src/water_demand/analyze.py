"""Physical reference demand and explicitly hypothetical crop-coefficient scenarios."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from water_demand.physics import hargreaves, penman_monteith
from water_demand.source import UNITS, read_sources


def run(root):
    config = json.loads((root / "configs/default.json").read_text())
    source = read_sources(root)
    frame = source[source.date.dt.year == config["year"]].copy()
    required = list(UNITS)
    if frame[required].isna().any().any() or frame.groupby("site").size().ne(366).any():
        raise ValueError("Complete leap-year 2024 inputs required; no silent imputation")
    result = penman_monteith(
        frame.T2M_MAX,
        frame.T2M_MIN,
        frame.RH2M,
        frame.WS2M,
        frame.ALLSKY_SFC_SW_DWN,
        frame.PS,
        frame.latitude,
        frame.date.dt.dayofyear,
        frame.elevation_m,
    )
    for key, value in result.items():
        frame[key] = value
    frame["et0_hargreaves_mm"] = hargreaves(frame.T2M_MAX, frame.T2M_MIN, frame.ra_mj_m2_day)
    frame["climatic_balance_mm"] = frame.PRECTOTCORR - frame.et0_pm_mm
    for kc in config["kc_scenarios"]:
        frame[f"etc_scenario_kc_{kc:g}_mm"] = kc * frame.et0_pm_mm
    scenario = [c for c in frame if c.startswith("etc_scenario")]
    monthly = (
        frame.assign(month=frame.date.dt.month)
        .groupby(["site", "month"])[
            ["et0_pm_mm", "et0_hargreaves_mm", "PRECTOTCORR", "climatic_balance_mm", *scenario]
        ]
        .sum()
        .reset_index()
    )
    comparisons = []
    for site, group in frame.groupby("site"):
        error = group.et0_hargreaves_mm - group.et0_pm_mm
        comparisons.append(
            {
                "site": site,
                "days": len(group),
                "pm_total_mm": float(group.et0_pm_mm.sum()),
                "hargreaves_total_mm": float(group.et0_hargreaves_mm.sum()),
                "mean_method_difference_mm_day": float(error.mean()),
                "method_mae_mm_day": float(abs(error).mean()),
                "method_rmse_mm_day": float(np.sqrt((error**2).mean())),
            }
        )
    comparison = pd.DataFrame(comparisons)
    reports = root / "reports"
    (reports / "figures").mkdir(parents=True, exist_ok=True)
    frame.to_csv(reports / "daily-water-demand.csv", index=False)
    monthly.to_csv(reports / "monthly-water-demand.csv", index=False)
    comparison.to_csv(reports / "method-comparison.csv", index=False)
    evidence = {
        "year": 2024,
        "rows": len(frame),
        "sites": sorted(frame.site.unique()),
        "cloud_ratio_bounded_days": int(frame.cloud_ratio_bounded.sum()),
        "emissivity_bounded_days": int(frame.emissivity_bounded.sum()),
        "negative_raw_et0_clipped_days": int((frame.raw_et0_pm_mm < 0).sum()),
        "config": config,
        "annual_totals": frame.groupby("site")[
            ["et0_pm_mm", "PRECTOTCORR", "climatic_balance_mm", *scenario]
        ]
        .sum()
        .to_dict(orient="index"),
        "scope": "Reference ET0 from gridded weather and hypothetical constant Kc scenarios; no observed crop ET or irrigation need",
    }
    (reports / "summary.json").write_text(json.dumps(evidence, indent=2, allow_nan=False))
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    for ax, (site, group) in zip(axes, monthly.groupby("site"), strict=True):
        ax.plot(group.month, group.et0_pm_mm, label="FAO-56 PM", color="#17684e", marker="o")
        ax.plot(
            group.month,
            group.et0_hargreaves_mm,
            label="Hargreaves",
            color="#cf872c",
            linestyle="--",
        )
        ax.bar(group.month, group.PRECTOTCORR, color="#7ea5c7", alpha=0.45, label="Precipitation")
        ax.set(title=site, xlabel="2024 month", xticks=range(1, 13))
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("Monthly reference ET0 or precipitation (mm)")
    axes[-1].legend(fontsize=9)
    fig.suptitle(
        "Gridded weather → reference evapotranspiration\nMethod comparison is disagreement, not measured accuracy or an irrigation prescription"
    )
    fig.tight_layout()
    fig.savefig(reports / "figures/monthly-reference-demand.png", dpi=150)
    plt.close(fig)
    return evidence


if __name__ == "__main__":
    print(json.dumps(run(Path.cwd()), indent=2))
