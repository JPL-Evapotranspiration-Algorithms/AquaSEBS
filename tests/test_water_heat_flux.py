"""
Unit tests for the AquaSEBS Equilibrium Temperature Model (ETM) water heat flux
parameterization (`AquaSEBS.water_heat_flux`).

References:
    - Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).
      Surface Energy Balance of Fresh and Saline Waters: AquaSEBS.
      Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
    - Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W.,
      Pearson, C., & Huntington, J. L. (2023). Remotely sensed terrestrial open
      water evaporation. Scientific Reports, 13, 8174.
      https://doi.org/10.1038/s41598-023-34921-2
"""

import numpy as np
import pytest
from AquaSEBS.water_heat_flux import water_heat_flux


def test_water_heat_flux_etm_reference_values():
    """Verify exact numerical agreement across all six ETM equations (Eqs. 8-13)."""
    WST_C = np.array([20.0])
    Td_C = np.array([10.0])
    windspeed_mps = np.array([3.0])
    SWnet = np.array([500.0])

    # Hand-calculated ETM equations (Abdelrady et al., 2016, Eqs. 8-13):
    # Tn = 0.5 * (20 - 10) = 5.0
    # eta = 0.35 + 0.015 * 20 + 0.0012 * (5.0 ** 2) = 0.35 + 0.30 + 0.03 = 0.68
    # S = 3.3 * 3.0 = 9.9
    # beta = 4.5 + 0.05 * 20 + (0.68 + 0.47) * 9.9 = 4.5 + 1.0 + 11.385 = 16.885
    # Te = 10.0 + (500.0 / 16.885) = 39.61208172934557
    # W_Wm2 = 16.885 * (39.61208172934557 - 20.0) = 331.15
    expected_Tn = 5.0
    expected_eta = 0.68
    expected_S = 9.9
    expected_beta = 16.885
    expected_Te = 10.0 + (500.0 / 16.885)
    expected_W = 331.15

    result = water_heat_flux(
        WST_C=WST_C,
        Td_C=Td_C,
        windspeed_mps=windspeed_mps,
        SWnet=SWnet,
        mask_non_water_pixels=False,
    )

    assert np.isclose(result["Tn"][0], expected_Tn, rtol=1e-12)
    assert np.isclose(result["η"][0], expected_eta, rtol=1e-12)
    assert np.isclose(result["S"][0], expected_S, rtol=1e-12)
    assert np.isclose(result["β"][0], expected_beta, rtol=1e-12)
    assert np.isclose(result["Te"][0], expected_Te, rtol=1e-12)
    assert np.isclose(result["W_Wm2"][0], expected_W, rtol=1e-12)


def test_water_heat_flux_sign_and_equilibrium_conditions():
    """
    Verify physical sign behavior:
    - Daytime warming (Te > WST_C) produces positive W_Wm2 (heat storage).
    - Nocturnal cooling (Te < WST_C, e.g., SWnet = 0) produces negative W_Wm2 (heat release).
    - Thermal equilibrium identity W = beta * (Te - WST_C).
    """
    # Daytime warming
    day = water_heat_flux(
        WST_C=np.array([18.0]),
        Td_C=np.array([12.0]),
        windspeed_mps=np.array([2.5]),
        SWnet=np.array([600.0]),
        mask_non_water_pixels=False,
    )
    assert day["Te"][0] > 18.0
    assert day["W_Wm2"][0] > 0.0

    # Nighttime / low-radiation cooling (preserves negative flux without clipping)
    night = water_heat_flux(
        WST_C=np.array([22.0]),
        Td_C=np.array([8.0]),
        windspeed_mps=np.array([4.0]),
        SWnet=np.array([0.0]),
        mask_non_water_pixels=False,
    )
    assert np.isclose(night["Te"][0], 8.0)
    assert night["W_Wm2"][0] < 0.0

    # Algebraic identity check at updated surface temperature
    eq_temp = float(day["Te"][0])
    eq_res = water_heat_flux(
        WST_C=np.array([eq_temp]),
        Td_C=np.array([12.0]),
        windspeed_mps=np.array([2.5]),
        SWnet=np.array([600.0]),
        mask_non_water_pixels=False,
    )
    recomputed_W = eq_res["β"][0] * (eq_res["Te"][0] - eq_res["WST_C"][0])
    assert np.isclose(eq_res["W_Wm2"][0], recomputed_W, atol=1e-12)


def test_water_heat_flux_2d_array_and_return_keys():
    """Verify 2D NumPy array propagation and complete dictionary key schema."""
    shape = (3, 4)
    WST_C = np.full(shape, 21.5)
    Td_C = np.full(shape, 11.0)
    windspeed_mps = np.full(shape, 3.5)
    SWnet = np.full(shape, 480.0)

    out = water_heat_flux(
        WST_C=WST_C,
        Td_C=Td_C,
        windspeed_mps=windspeed_mps,
        SWnet=SWnet,
        mask_non_water_pixels=False,
    )

    expected_keys = {
        "WST_C",
        "Ta_C",
        "RH",
        "Td_C",
        "windspeed_mps",
        "SWnet",
        "Tn",
        "η",
        "S",
        "β",
        "Te",
        "W_Wm2",
    }
    assert expected_keys.issubset(out.keys())
    for key in ("Tn", "η", "S", "β", "Te", "W_Wm2"):
        assert out[key].shape == shape
        assert np.all(np.isfinite(out[key]))


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (
            dict(WST_C=None, Td_C=np.array([10.0]), windspeed_mps=np.array([3.0]), SWnet=np.array([500.0])),
            "WST_C",
        ),
        (
            dict(WST_C=np.array([20.0]), Td_C=None, windspeed_mps=np.array([3.0]), SWnet=np.array([500.0])),
            "Td_C",
        ),
        (
            dict(WST_C=np.array([20.0]), Td_C=np.array([10.0]), windspeed_mps=None, SWnet=np.array([500.0])),
            "windspeed_mps",
        ),
        (
            dict(WST_C=np.array([20.0]), Td_C=np.array([10.0]), windspeed_mps=np.array([3.0]), SWnet=None),
            "SWnet",
        ),
    ],
)
def test_water_heat_flux_missing_inputs_raise_value_error(kwargs, match):
    """Verify ValueError is raised when required inputs cannot be resolved."""
    with pytest.raises(ValueError, match=match):
        water_heat_flux(mask_non_water_pixels=False, **kwargs)
