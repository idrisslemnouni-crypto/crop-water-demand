"""Fixed-root-zone daily scenarios in millimetres; no calibrated irrigation advice."""

import math
from dataclasses import dataclass


def finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Finite numeric {name} required")
    return float(value)


@dataclass(frozen=True)
class RootZone:
    theta_fc: float
    theta_wp: float
    root_depth_m: float
    depletion_fraction: float

    def __post_init__(self):
        for name in ("theta_fc", "theta_wp", "root_depth_m", "depletion_fraction"):
            finite_number(getattr(self, name), name)
        if not 0 <= self.theta_wp < self.theta_fc <= 1 or self.root_depth_m <= 0:
            raise ValueError(
                "Require 0 <= wilting point < field capacity <= 1 and positive root depth"
            )
        if not 0 < self.depletion_fraction < 1:
            raise ValueError("Depletion fraction must be strictly between zero and one")
        if not math.isfinite(self.taw_mm) or self.taw_mm <= 0:
            raise ValueError("Finite positive root-zone capacity required")

    @property
    def taw_mm(self):
        return 1000 * (self.theta_fc - self.theta_wp) * self.root_depth_m

    @property
    def raw_mm(self):
        return self.depletion_fraction * self.taw_mm


def daily_balance(
    zone: RootZone,
    depletion_mm,
    et0_mm,
    kc,
    precipitation_mm=0.0,
    runoff_mm=0.0,
    net_irrigation_mm=0.0,
    capillary_rise_mm=0.0,
):
    """FAO-56 Eq.82–86/88 with beginning-of-day stress and same-day drainage.

    Water supplied is net infiltration; losses/efficiency outside the root zone are not inferred.
    Single Kc stress approximation and fixed roots/p apply. Explicitly cap ET by available water.
    """
    names = (
        "depletion",
        "et0",
        "kc",
        "precipitation",
        "runoff",
        "net_irrigation",
        "capillary_rise",
    )
    d, et0, coefficient, rain, runoff, irrigation, capillary = [
        finite_number(value, name)
        for value, name in zip(
            (
                depletion_mm,
                et0_mm,
                kc,
                precipitation_mm,
                runoff_mm,
                net_irrigation_mm,
                capillary_rise_mm,
            ),
            names,
            strict=True,
        )
    ]
    if not 0 <= d <= zone.taw_mm or min(et0, coefficient, rain, runoff, irrigation, capillary) < 0:
        raise ValueError("Depletion must be within TAW; physical fluxes and Kc must be nonnegative")
    if runoff > rain:
        raise ValueError("Runoff cannot exceed precipitation")
    ks = 1.0 if d <= zone.raw_mm else (zone.taw_mm - d) / (zone.taw_mm - zone.raw_mm)
    potential = et0 * coefficient
    inflow = rain - runoff + irrigation + capillary
    if not math.isfinite(potential) or not math.isfinite(inflow):
        raise ValueError("Flux calculation overflow")
    actual = min(ks * potential, zone.taw_mm - d + inflow)
    drainage = max(inflow - actual - d, 0.0)
    end = d - inflow + actual + drainage
    # Suppress floating-point excursions only, never discard real water losses.
    end = min(max(end, 0.0), zone.taw_mm)
    return {
        "depletion_start_mm": d,
        "taw_mm": zone.taw_mm,
        "raw_mm": zone.raw_mm,
        "ks_start": ks,
        "potential_etc_mm": potential,
        "adjusted_etc_scenario_mm": actual,
        "deep_percolation_mm": drainage,
        "depletion_end_mm": end,
        "storage_end_mm": zone.taw_mm - end,
        "net_irrigation_mm": irrigation,
        "unmet_potential_etc_mm": potential - actual,
        "balance_residual_mm": end - d + inflow - actual - drainage,
    }
