"""Independent published examples and physical conservation contracts."""

import numpy as np
import pytest

from water_demand.balance import RootZone, daily_balance


def test_fao_example_36_available_water():
    # Maize/silt combination in the independent rounded reference table.
    zone = RootZone(0.32, 0.15, 1.2, 0.55)
    assert zone.taw_mm == pytest.approx(204)
    assert zone.raw_mm == pytest.approx(112, abs=0.3)


def test_fao_example_37_ten_day_stress_sequence():
    zone = RootZone(0.32, 0.12, 0.8, 0.4)
    expected = [61.0, 67.0, 72.8, 78.3, 83.4, 88.2, 92.6, 96.9, 100.8, 104.5]
    depletion = 55.0
    for rounded in expected:
        result = daily_balance(zone, depletion, 5, 1.2)
        assert result["depletion_end_mm"] == pytest.approx(rounded, abs=0.08)
        assert abs(result["balance_residual_mm"]) < 1e-10
        depletion = result["depletion_end_mm"]


def test_rain_drainage_and_wilting_balance():
    zone = RootZone(0.3, 0.1, 1, 0.5)
    wet = daily_balance(zone, 10, 5, 1, precipitation_mm=40, runoff_mm=5)
    assert wet["depletion_end_mm"] == 0
    assert wet["deep_percolation_mm"] == 20
    dry = daily_balance(zone, zone.taw_mm, 10, 1)
    assert dry["ks_start"] == dry["adjusted_etc_scenario_mm"] == 0
    assert dry["depletion_end_mm"] == zone.taw_mm
    assert wet["balance_residual_mm"] == dry["balance_residual_mm"] == 0


def test_large_demand_cannot_extract_more_than_stored_water():
    zone = RootZone(0.3, 0.1, 1, 0.5)
    result = daily_balance(zone, 120, 1000, 1)
    assert result["adjusted_etc_scenario_mm"] == pytest.approx(80)
    assert result["depletion_end_mm"] == pytest.approx(zone.taw_mm)
    assert result["balance_residual_mm"] == pytest.approx(0)


def test_derived_flux_overflow_rejected():
    with pytest.raises(ValueError, match="overflow"):
        daily_balance(RootZone(0.3, 0.1, 1, 0.5), 0, 1e308, 1e308)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"theta_fc": np.nan},
        {"theta_wp": 0.4},
        {"root_depth_m": -1},
        {"depletion_fraction": 1},
        {"theta_fc": True},
    ],
)
def test_invalid_soil_contract(kwargs):
    params = dict(theta_fc=0.3, theta_wp=0.1, root_depth_m=1, depletion_fraction=0.5) | kwargs
    with pytest.raises(ValueError):
        RootZone(**params)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"depletion_mm": -1},
        {"kc": np.nan},
        {"precipitation_mm": -1},
        {"runoff_mm": 1},
        {"et0_mm": True},
    ],
)
def test_invalid_flux_contract(kwargs):
    with pytest.raises(ValueError):
        daily_balance(RootZone(0.3, 0.1, 1, 0.5), **(dict(depletion_mm=0, et0_mm=5, kc=1) | kwargs))
