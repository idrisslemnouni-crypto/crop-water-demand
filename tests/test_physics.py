import numpy as np
import pytest

from water_demand.physics import (
    extraterrestrial_radiation,
    hargreaves,
    penman_monteith,
    saturation_pressure,
)


def test_fao_published_daily_example_18():
    # Published Uccle 6 July case, FAO-56 Ch.4 example 18; rounded inputs.
    ea = (saturation_pressure(12.3) * 0.84 + saturation_pressure(21.5) * 0.63) / 2
    result = penman_monteith(21.5, 12.3, 73.5, 2.078, 22.07, 100.1, 50.8, 187, 100, ea)
    assert result["ra_mj_m2_day"] == pytest.approx(41.09, abs=0.08)
    assert result["et0_pm_mm"] == pytest.approx(3.88, abs=0.03)


def test_fao_hargreaves_example_20_and_energy_conversion():
    # Lyon July example 20, FAO-56 Ch.4: reference result 5.0 mm/day.
    assert hargreaves(26.6, 14.8, 40.55) == pytest.approx(5.0, abs=0.1)


def test_wrong_units_and_nonfinite_inputs_rejected():
    with pytest.raises(ValueError):
        penman_monteith(25, 15, 150, 2, 20, 101, 33, 180, 100)
    with pytest.raises(ValueError):
        penman_monteith(25, 15, 50, 2, np.nan, 101, 33, 180, 100)
    with pytest.raises(ValueError):
        hargreaves(15, 25, 40)
    with pytest.raises(ValueError):
        extraterrestrial_radiation(85, 180)


def test_explicit_radiation_bounds_and_nonnegative_output():
    result = penman_monteith(25, 15, 50, 2, 100, 101, 33, 180, 100)
    assert result["cloud_ratio_bounded"] and result["et0_pm_mm"] >= 0
