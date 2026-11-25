"""
Unit tests for spatial water masking, Raster geometry inference, and daylight
evapotranspiration upscaling (`upscale_to_daylight=True`).
"""

from datetime import datetime
import numpy as np
from rasters import Raster, RasterGrid
from AquaSEBS import AquaSEBS


def test_water_masking_excludes_land_pixels():
    """Verify non-water pixels become NaN when mask_non_water_pixels=True."""
    water_mask = np.array([[True, False], [False, True]])
    WST_C = np.array([[22.0, 35.0], [34.0, 21.0]])
    Ta_C = np.full((2, 2), 25.0)
    Td_C = np.full((2, 2), 14.0)
    windspeed_mps = np.full((2, 2), 3.0)
    SWnet = np.full((2, 2), 550.0)
    Rn_Wm2 = np.full((2, 2), 500.0)

    masked_results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        Td_C=Td_C,
        windspeed_mps=windspeed_mps,
        SWnet=SWnet,
        Rn_Wm2=Rn_Wm2,
        water=water_mask,
        mask_non_water_pixels=True,
        offline_mode=True,
    )

    for key in ("W_Wm2", "LE_Wm2", "Tn", "η", "β", "Te"):
        arr = np.asarray(masked_results[key])
        assert np.isfinite(arr[0, 0]) and np.isfinite(arr[1, 1])
        assert np.isnan(arr[0, 1]) and np.isnan(arr[1, 0])

    unmasked_results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        Td_C=Td_C,
        windspeed_mps=windspeed_mps,
        SWnet=SWnet,
        Rn_Wm2=Rn_Wm2,
        water=water_mask,
        mask_non_water_pixels=False,
        offline_mode=True,
    )
    assert np.all(np.isfinite(unmasked_results["LE_Wm2"]))


def test_raster_geometry_inference_and_daylight_upscaling():
    """Verify automatic RasterGeometry inference and daylight ET upscaling."""
    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    shape = geometry.shape
    water_raster = Raster(np.array([[1, 0], [1, 1]], dtype=bool), geometry=geometry)
    wst_raster = Raster(np.full(shape, 22.0), geometry=geometry)

    # Solar noon (~19:30 UTC at -120° longitude)
    time_UTC = datetime(2023, 7, 15, 19, 30)

    results = AquaSEBS(
        WST_C=wst_raster,
        Ta_C=Raster(np.full(shape, 25.0), geometry=geometry),
        Td_C=Raster(np.full(shape, 14.0), geometry=geometry),
        windspeed_mps=Raster(np.full(shape, 3.0), geometry=geometry),
        SWnet=Raster(np.full(shape, 600.0), geometry=geometry),
        Rn_Wm2=Raster(np.full(shape, 550.0), geometry=geometry),
        water=water_raster,
        time_UTC=time_UTC,
        upscale_to_daylight=True,
        mask_non_water_pixels=True,
        offline_mode=True,
    )

    assert "LE_Wm2" in results
    assert any("daylight" in k.lower() or "et_" in k.lower() for k in results.keys())

    # Verify water pixels are positive and masked land pixel [0, 1] remains NaN
    le_arr = np.asarray(results["LE_Wm2"])
    assert le_arr[0, 0] > 0.0
    assert np.isnan(le_arr[0, 1])
