"""
Unit tests for the main AquaSEBS open-water surface energy balance model
(`AquaSEBS.AquaSEBS`).
"""

from datetime import datetime
import numpy as np
import pytest
from priestley_taylor import epsilon_from_Ta_C, GAMMA_PA, PT_ALPHA
from rasters import RasterGrid
from AquaSEBS import AquaSEBS


def test_AquaSEBS_offline_energy_balance_matches_priestley_taylor():
    """Verify end-to-end LE_Wm2 matches alpha * epsilon * (Rn - W) offline."""
    WST_C = np.array([[22.0, 24.5]])
    Ta_C = np.array([[25.0, 26.0]])
    Td_C = np.array([[14.0, 15.0]])
    windspeed_mps = np.array([[3.2, 4.1]])
    SWnet = np.array([[560.0, 580.0]])
    Rn_Wm2 = np.array([[510.0, 530.0]])

    results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        Td_C=Td_C,
        windspeed_mps=windspeed_mps,
        SWnet=SWnet,
        Rn_Wm2=Rn_Wm2,
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    expected_epsilon = epsilon_from_Ta_C(Ta_C=Ta_C, gamma_Pa=GAMMA_PA)
    expected_LE = PT_ALPHA * expected_epsilon * (Rn_Wm2 - results["W_Wm2"])

    assert "LE_Wm2" in results
    assert "W_Wm2" in results
    assert "G_Wm2" in results
    assert np.allclose(results["epsilon"], expected_epsilon, rtol=1e-12)
    assert np.allclose(results["LE_Wm2"], expected_LE, rtol=1e-12)


def test_AquaSEBS_custom_alpha_and_gamma():
    """Verify custom Priestley-Taylor alpha and psychrometric constant gamma scale LE_Wm2."""
    WST_C = np.array([[20.0]])
    Ta_C = np.array([[22.0]])
    Rn_Wm2 = np.array([[500.0]])
    W_Wm2 = np.array([[200.0]])
    custom_alpha = 1.15
    custom_gamma = 64.0

    results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        Rn_Wm2=Rn_Wm2,
        W_Wm2=W_Wm2,
        α=custom_alpha,
        γ_Pa=custom_gamma,
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    expected_epsilon = epsilon_from_Ta_C(Ta_C=Ta_C, gamma_Pa=custom_gamma)
    expected_LE = custom_alpha * expected_epsilon * (Rn_Wm2 - W_Wm2)

    assert np.allclose(results["LE_Wm2"], expected_LE, rtol=1e-12)


def test_AquaSEBS_precomputed_water_heat_flux_passthrough():
    """Verify providing W_Wm2 directly bypasses water_heat_flux() calculation."""
    WST_C = np.array([[21.0]])
    Ta_C = np.array([[23.0]])
    Rn_Wm2 = np.array([[450.0]])
    W_Wm2 = np.array([[150.0]])

    results = AquaSEBS(
        WST_C=WST_C,
        Ta_C=Ta_C,
        Rn_Wm2=Rn_Wm2,
        W_Wm2=W_Wm2,
        mask_non_water_pixels=False,
        offline_mode=True,
    )

    assert np.allclose(results["G_Wm2"], W_Wm2)
    assert "Tn" not in results
    assert "LE_Wm2" in results


def test_AquaSEBS_verma_net_radiation_integration():
    """Verify Rn_Wm2 is computed via verma_net_radiation when radiation inputs are supplied."""
    geometry = RasterGrid.from_bbox(
        xmin=-120.0, ymin=35.0, xmax=-119.98, ymax=35.02, cell_size=0.01
    )
    shape = geometry.shape
    time_UTC = datetime(2023, 7, 15, 18, 0)

    results = AquaSEBS(
        WST_C=np.full(shape, 22.0),
        Ta_C=np.full(shape, 25.0),
        RH=np.full(shape, 0.45),
        Td_C=np.full(shape, 12.0),
        windspeed_mps=np.full(shape, 3.0),
        SWin_Wm2=np.full(shape, 750.0),
        SWnet=np.full(shape, 705.0),
        albedo=np.full(shape, 0.06),
        emissivity=np.full(shape, 0.985),
        water=np.ones(shape, dtype=bool),
        geometry=geometry,
        time_UTC=time_UTC,
        mask_non_water_pixels=True,
        offline_mode=True,
    )

    assert "Rn_Wm2" in results
    assert results["Rn_Wm2"].shape == shape
    assert np.all(np.isfinite(results["Rn_Wm2"]))
    assert np.all(np.isfinite(results["LE_Wm2"]))


def test_AquaSEBS_missing_net_radiation_raises_value_error():
    """Verify ValueError when neither Rn_Wm2 nor sufficient radiation inputs are provided."""
    with pytest.raises(ValueError, match="Rn_Wm2"):
        AquaSEBS(
            WST_C=np.array([[20.0]]),
            Ta_C=np.array([[22.0]]),
            W_Wm2=np.array([[100.0]]),
            mask_non_water_pixels=False,
            offline_mode=True,
        )
