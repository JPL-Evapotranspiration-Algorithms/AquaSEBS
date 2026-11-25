"""
Unit tests verifying offline_mode enforcement and zero external network calls
when all required forcings are provided.
"""

from datetime import datetime
import importlib
import numpy as np
import pytest
from rasters import RasterGrid
from AquaSEBS import AquaSEBS, MissingOfflineParameter, water_heat_flux

aquasebs_mod = importlib.import_module("AquaSEBS.AquaSEBS")
whf_mod = importlib.import_module("AquaSEBS.water_heat_flux")


def test_offline_mode_raises_when_Ta_C_missing():
    """Verify MissingOfflineParameter is raised when offline_mode=True and Ta_C is omitted."""
    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    time_UTC = datetime(2023, 7, 15, 18, 0)

    with pytest.raises(MissingOfflineParameter, match="Ta_C"):
        AquaSEBS(
            WST_C=np.full(geometry.shape, 22.0),
            Rn_Wm2=np.full(geometry.shape, 500.0),
            water=np.ones(geometry.shape, dtype=bool),
            geometry=geometry,
            time_UTC=time_UTC,
            offline_mode=True,
        )


def test_offline_mode_raises_when_SWin_Wm2_missing():
    """Verify MissingOfflineParameter is raised when offline_mode=True and SWin_Wm2 is omitted."""
    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    time_UTC = datetime(2023, 7, 15, 18, 0)

    with pytest.raises(MissingOfflineParameter, match="SWin_Wm2"):
        AquaSEBS(
            WST_C=np.full(geometry.shape, 22.0),
            Ta_C=np.full(geometry.shape, 25.0),
            albedo=np.full(geometry.shape, 0.06),
            emissivity=np.full(geometry.shape, 0.985),
            water=np.ones(geometry.shape, dtype=bool),
            geometry=geometry,
            time_UTC=time_UTC,
            offline_mode=True,
        )


def test_offline_mode_raises_when_water_mask_missing_and_masking_enabled():
    """Verify MissingOfflineParameter is raised when mask_non_water_pixels=True and water is omitted offline."""
    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )

    with pytest.raises(MissingOfflineParameter, match="water"):
        AquaSEBS(
            WST_C=np.full(geometry.shape, 22.0),
            Ta_C=np.full(geometry.shape, 25.0),
            Td_C=np.full(geometry.shape, 14.0),
            windspeed_mps=np.full(geometry.shape, 3.0),
            SWnet=np.full(geometry.shape, 550.0),
            Rn_Wm2=np.full(geometry.shape, 500.0),
            geometry=geometry,
            mask_non_water_pixels=True,
            offline_mode=True,
        )


def test_water_heat_flux_offline_mode_raises_missing_offline_parameter():
    """Verify water_heat_flux() raises MissingOfflineParameter when required forcings are omitted offline."""
    with pytest.raises(MissingOfflineParameter, match="windspeed_mps"):
        water_heat_flux(
            WST_C=np.array([20.0]),
            Td_C=np.array([10.0]),
            SWnet=np.array([500.0]),
            mask_non_water_pixels=False,
            offline_mode=True,
        )


def test_offline_mode_does_not_instantiate_GEOS5FP(monkeypatch):
    """Verify GEOS5FP() is never instantiated when offline_mode=True or all inputs are supplied."""
    def fail_init(*args, **kwargs):
        raise AssertionError("GEOS5FP() should not be instantiated when offline or inputs are provided")

    monkeypatch.setattr(aquasebs_mod, "GEOS5FP", fail_init)
    monkeypatch.setattr(whf_mod, "GEOS5FP", fail_init)

    results = AquaSEBS(
        WST_C=np.array([[22.0]]),
        Ta_C=np.array([[25.0]]),
        RH=np.array([[0.50]]),
        windspeed_mps=np.array([[3.0]]),
        SWin_Wm2=np.array([[600.0]]),
        albedo=np.array([[0.06]]),
        Rn_Wm2=np.array([[520.0]]),
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    assert np.all(np.isfinite(results["LE_Wm2"]))
    assert np.all(np.isfinite(results["W_Wm2"]))


def test_AquaSEBS_forwards_RH_SWin_and_albedo_to_water_heat_flux():
    """Verify AquaSEBS() forwards RH, SWin_Wm2, and albedo to water_heat_flux() without network calls."""
    WST_C = np.array([[22.0]])
    Ta_C = np.array([[25.0]])
    RH = np.array([[0.50]])
    windspeed_mps = np.array([[3.0]])
    SWin_Wm2 = np.array([[600.0]])
    albedo = np.array([[0.06]])
    Rn_Wm2 = np.array([[520.0]])

    results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        RH=RH,
        windspeed_mps=windspeed_mps,
        SWin_Wm2=SWin_Wm2,
        albedo=albedo,
        Rn_Wm2=Rn_Wm2,
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    expected_whf = water_heat_flux(
        WST_C=WST_C,
        Ta_C=Ta_C,
        RH=RH,
        windspeed_mps=windspeed_mps,
        SWin_Wm2=SWin_Wm2,
        albedo=albedo,
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    assert np.allclose(results["W_Wm2"], expected_whf["W_Wm2"])
    assert np.allclose(results["SWnet"], SWin_Wm2 * (1.0 - albedo))


def test_offline_execution_makes_zero_network_calls(monkeypatch):
    """Verify full offline execution never triggers GEOS5FP or NASADEM.swb lookups."""

    def fail_network(*args, **kwargs):
        raise AssertionError("External network lookup was invoked during offline execution")

    monkeypatch.setattr(aquasebs_mod.NASADEM, "swb", fail_network)
    monkeypatch.setattr(whf_mod.NASADEM, "swb", fail_network)

    class FailingGEOS5FP:
        def Ta_C(self, *args, **kwargs):
            fail_network()

        def RH(self, *args, **kwargs):
            fail_network()

        def wind_speed_mps(self, *args, **kwargs):
            fail_network()

        def SWin(self, *args, **kwargs):
            fail_network()

        def ALBEDO(self, *args, **kwargs):
            fail_network()

    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    shape = geometry.shape

    results = AquaSEBS(
        WST_C=np.full(shape, 21.0),
        Ta_C=np.full(shape, 24.0),
        Td_C=np.full(shape, 12.0),
        windspeed_mps=np.full(shape, 3.0),
        SWnet=np.full(shape, 520.0),
        Rn_Wm2=np.full(shape, 480.0),
        water=np.ones(shape, dtype=bool),
        geometry=geometry,
        GEOS5FP_connection=FailingGEOS5FP(),
        mask_non_water_pixels=True,
        offline_mode=True,
    )

    assert np.all(np.isfinite(results["LE_Wm2"]))
    assert np.all(np.isfinite(results["W_Wm2"]))


def test_water_heat_flux_uses_geos5fp_wind_speed_mps():
    """Verify online windspeed lookup invokes GEOS5FP.wind_speed_mps (Issue #15)."""

    class MockGEOS5FP:
        def __init__(self):
            self.called = False

        def wind_speed_mps(self, time_UTC, geometry, resampling):
            self.called = True
            return np.full(geometry.shape, 4.2)

    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    mock_geos = MockGEOS5FP()

    out = water_heat_flux(
        WST_C=np.full(geometry.shape, 20.0),
        Td_C=np.full(geometry.shape, 10.0),
        SWnet=np.full(geometry.shape, 500.0),
        geometry=geometry,
        time_UTC=datetime(2023, 7, 15, 18, 0),
        GEOS5FP_connection=mock_geos,
        mask_non_water_pixels=False,
        offline_mode=False,
    )

    assert mock_geos.called
    assert np.allclose(out["windspeed_mps"], 4.2)
