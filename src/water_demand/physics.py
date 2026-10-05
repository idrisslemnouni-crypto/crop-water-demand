"""Daily FAO-56 reference calculations, with explicit approximations and units."""

import numpy as np


def saturation_pressure(t_c):
    return 0.6108 * np.exp(17.27 * np.asarray(t_c) / (np.asarray(t_c) + 237.3))


def extraterrestrial_radiation(latitude_deg, day_of_year):
    latitude = np.asarray(latitude_deg, dtype=float)
    j = np.asarray(day_of_year, dtype=float)
    if (
        (abs(latitude) > 60).any()
        or ((j < 1) | (j > 366)).any()
        or not np.isfinite(latitude).all()
        or not np.isfinite(j).all()
    ):
        raise ValueError("Supported nonpolar latitude and day-of-year required")
    phi = np.deg2rad(latitude)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365)
    declination = 0.409 * np.sin(2 * np.pi * j / 365 - 1.39)
    sunset = np.arccos(np.clip(-np.tan(phi) * np.tan(declination), -1, 1))
    return (
        24
        * 60
        / np.pi
        * 0.0820
        * dr
        * (
            sunset * np.sin(phi) * np.sin(declination)
            + np.cos(phi) * np.cos(declination) * np.sin(sunset)
        )
    )


def penman_monteith(
    tmax_c,
    tmin_c,
    rh_percent,
    wind_2m_ms,
    solar_mj_m2_day,
    pressure_kpa,
    latitude_deg,
    day_of_year,
    elevation_m,
    actual_vapour_pressure_kpa=None,
):
    tmax, tmin, rh, wind, solar, pressure, elevation = [
        np.asarray(v, dtype=float)
        for v in [
            tmax_c,
            tmin_c,
            rh_percent,
            wind_2m_ms,
            solar_mj_m2_day,
            pressure_kpa,
            elevation_m,
        ]
    ]
    for value in [tmax, tmin, rh, wind, solar, pressure, elevation]:
        if not np.isfinite(value).all():
            raise ValueError("Finite physical meteorological inputs required")
    if (
        (tmin > tmax).any()
        or (tmin < -60).any()
        or (tmax > 60).any()
        or ((rh < 0) | (rh > 100)).any()
        or (wind < 0).any()
        or (solar < 0).any()
        or (pressure <= 0).any()
        or ((elevation < -500) | (elevation > 5000)).any()
    ):
        raise ValueError("Unsupported physical input range or units")
    t = (tmax + tmin) / 2
    es = (saturation_pressure(tmax) + saturation_pressure(tmin)) / 2
    ea = (
        es * rh / 100
        if actual_vapour_pressure_kpa is None
        else np.asarray(actual_vapour_pressure_kpa, dtype=float)
    )
    if not np.isfinite(ea).all() or (ea < 0).any() or (ea > es).any():
        raise ValueError("Invalid actual vapour pressure")
    slope = 4098 * saturation_pressure(t) / (t + 237.3) ** 2
    gamma = 0.000665 * pressure
    ra = extraterrestrial_radiation(latitude_deg, day_of_year)
    rso = (0.75 + 2e-5 * elevation) * ra
    ratio = solar / rso
    bounded = np.clip(ratio, 0.3, 1)
    emissivity = np.maximum(0.34 - 0.14 * np.sqrt(ea), 0)
    rnl = (
        4.903e-9
        * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4)
        / 2
        * emissivity
        * (1.35 * bounded - 0.35)
    )
    rn = 0.77 * solar - rnl
    raw = (0.408 * slope * rn + gamma * 900 / (t + 273) * wind * (es - ea)) / (
        slope + gamma * (1 + 0.34 * wind)
    )
    return {
        "et0_pm_mm": np.maximum(raw, 0),
        "raw_et0_pm_mm": raw,
        "ra_mj_m2_day": ra,
        "net_radiation_mj_m2_day": rn,
        "cloud_ratio": ratio,
        "cloud_ratio_bounded": (ratio < 0.3) | (ratio > 1),
        "emissivity_bounded": (0.34 - 0.14 * np.sqrt(ea) < 0),
        "ea_estimated_kpa": ea,
    }


def hargreaves(tmax_c, tmin_c, ra_mj_m2_day):
    tmax, tmin, ra = [np.asarray(v, dtype=float) for v in [tmax_c, tmin_c, ra_mj_m2_day]]
    if (
        not all(np.isfinite(v).all() for v in [tmax, tmin, ra])
        or (tmin > tmax).any()
        or (ra < 0).any()
    ):
        raise ValueError("Finite ordered temperatures and radiation required")
    # FAO Eq.52 uses radiation expressed as equivalent evaporation, mm/day.
    return np.maximum(0.0023 * ((tmax + tmin) / 2 + 17.8) * np.sqrt(tmax - tmin) * (0.408 * ra), 0)
