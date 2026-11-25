"""
Equilibrium Temperature Model (ETM) for Open-Water Heat Flux in AquaSEBS.

This module implements the remote-sensing water heat storage and exchange flux
(`W_Wm2`, denoted `G0w` in the literature) parameterization used in the Aquatic
Surface Energy Balance System (AquaSEBS). Because satellite thermal infrared
sensors observe only the water skin temperature rather than vertical depth
temperature profiles, AquaSEBS estimates instantaneous air-water heat flux
via the Equilibrium Temperature Model (Edinger et al., 1968; Abdelrady et al.,
2016; Fisher et al., 2023).

Primary References:
    - Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).
      Surface Energy Balance of Fresh and Saline Waters: AquaSEBS.
      Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
    - Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W.,
      Pearson, C., & Huntington, J. L. (2023). Remotely sensed terrestrial open
      water evaporation. Scientific Reports, 13, 8174.
      https://doi.org/10.1038/s41598-023-34921-2
    - Edinger, J. E., Duttweiler, D. W., & Geyer, J. C. (1968). The response of
      water temperatures to meteorological conditions. Water Resources Research,
      4(5), 1137-1143. https://doi.org/10.1029/WR004i005p01137
    - Abualnaja, Y. (2009). Estimation of the net surface heat flux in the
      Arabian Gulf based on the equilibrium temperature. Journal of King
      Abdulaziz University: Marine Sciences, 20, 21-29.
    - Ahmad, F., & Sultan, S. A. R. (1994). Equilibrium temperature as a
      parameter for estimating the net heat-flux at the air-sea interface in
      the central Red Sea. Oceanologica Acta, 17(3), 341-343.
"""

from typing import Union, Dict
import logging
import numpy as np
from datetime import datetime
import rasters as rt
from rasters import Raster, RasterGeometry
from NASADEM import NASADEM
from GEOS5FP import GEOS5FP
from check_distribution import check_distribution

from .constants import (
    RESAMPLING_METHOD,
    MASK_NON_WATER_PIXELS,
    DEWPOINT_METHOD,
    HIGH_WIND_THRESHOLD_MPS,
)
from .exceptions import MissingOfflineParameter

logger = logging.getLogger(__name__)


def _nan_like(template: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """Create an all-NaN array or Raster matching the shape, dtype, and geometry of `template`."""
    arr = np.asarray(template)
    dtype = arr.dtype if np.issubdtype(arr.dtype, np.floating) else np.float32
    nan_arr = np.full(arr.shape, np.nan, dtype=dtype)
    if isinstance(template, Raster):
        return Raster(nan_arr, geometry=template.geometry)
    return nan_arr


def water_heat_flux(
        WST_C: Union[Raster, np.ndarray],
        albedo: Union[Raster, np.ndarray] = None,
        Ta_C: Union[Raster, np.ndarray] = None,
        RH: Union[Raster, np.ndarray] = None,
        Td_C: Union[Raster, np.ndarray] = None,
        windspeed_mps: Union[Raster, np.ndarray] = None,
        SWnet: Union[Raster, np.ndarray] = None,
        SWin_Wm2: Union[Raster, np.ndarray] = None,
        geometry: RasterGeometry = None,
        time_UTC: datetime = None,
        water: Union[Raster, np.ndarray] = None,
        GEOS5FP_connection: GEOS5FP = None,
        resampling: str = RESAMPLING_METHOD,
        mask_non_water_pixels: bool = MASK_NON_WATER_PIXELS,
        dewpoint_method: str = DEWPOINT_METHOD,
        offline_mode: bool = False) -> Dict[str, Union[Raster, np.ndarray, float]]:
    """
    Calculate instantaneous water heat flux (`W_Wm2`) using the AquaSEBS
    Equilibrium Temperature Model (ETM).

    Scientific Methodology:
        In situ water heat storage flux is defined vertically as the depth-integrated
        rate of thermal energy change:

            G0w = rho_w * c_w * integral_0^Z (dT / dt) dZ

        (Abdelrady et al., 2016, Eq. 7). Because spaceborne thermal infrared
        instruments only measure skin water surface temperature (`WST_C`, T0),
        AquaSEBS employs the Equilibrium Temperature Model (Edinger et al., 1968;
        Abdelrady et al., 2016, Section 2.2.1, Eqs. 8-13; Fisher et al., 2023,
        Methods, Eqs. 1-6).

        The equilibrium temperature (`Te`) is defined as the hypothetical water
        surface temperature at which net heat exchange across the air-water
        interface equals zero:

        1. **Half Temperature Gradient (`Tn`, °C)** (Abdelrady et al., 2016, Eq. 8):
               Tn = 0.5 * (WST_C - Td_C)
        2. **Evaporation Efficiency (`η`, dimensionless)** (Abdelrady et al., 2016, Eq. 9):
               η = 0.35 + 0.015 * WST_C + 0.0012 * Tn^2
        3. **Wind Function (`S`, m s^-1)** (Abdelrady et al., 2016, Eq. 10):
               S = 3.3 * windspeed_mps
        4. **Thermal Exchange Coefficient (`β`, W m^-2 °C^-1)** (Abdelrady et al., 2016, Eq. 11):
               β = 4.5 + 0.05 * WST_C + (η + 0.47) * S
        5. **Equilibrium Temperature (`Te`, °C)** (Abdelrady et al., 2016, Eq. 12):
               Te = Td_C + (SWnet / β)
        6. **Water Heat Flux (`W_Wm2`, W m^-2)** (Abdelrady et al., 2016, Eq. 13):
               W_Wm2 = β * (Te - WST_C)

        Positive `W_Wm2` indicates net radiative/sensible heating stored into the
        water column (typical during daytime solar loading), while negative `W_Wm2`
        indicates heat released from water storage back to the atmosphere (typical
        at night or during cold-air outbreaks).

    Scientific Caveats & Sensitivity Benchmarks:
        - **Dominant Error Source**: Abdelrady et al. (2016, Section 4.5-4.6)
          demonstrated that total AquaSEBS flux uncertainty propagates primarily
          from water heat flux estimation (rRMSE = 4.3-8.5% against ECMWF across
          five marine basins). Sensitivity ranks from highest to lowest as:
          (1) water surface temperature (`WST_C`), (2) dew-point temperature
          (`Td_C`), and (3) wind speed (`windspeed_mps`). Thermal sensors with
          accuracy < 0.3 K are recommended (Abdelrady et al., 2016, Section 5).
        - **Dew-Point Approximation**: When `Td_C` is omitted, it is estimated
          from `Ta_C` and `RH` using either `dewpoint_method="linear"` (default,
          `Td_C = Ta_C - ((100 - RH * 100) / 5.0)`, preserving ECOSTRESS
          Collection 3 baseline parity) or `dewpoint_method="magnus_tetens"`
          (WMO Magnus-Tetens saturation vapor pressure formula, recommended
          under arid conditions where RH < 50%).
        - **High-Wind Regime Shift**: Instantaneous open-water evaporation shifts
          from radiatively controlled to aerodynamically controlled above
          `HIGH_WIND_THRESHOLD_MPS = 7.5` m s^-1 (Fisher et al., 2023).
          A boolean `"high_wind_flag"` mask (`windspeed_mps > 7.5`) is returned
          as diagnostic metadata without modifying flux values.
        - **Bathymetry & Mixing Depth**: ETM does not explicitly parameterize
          water column depth or vertical mixing layer thickness, which governs
          diurnal heat storage hysteresis in shallow lakes versus deep oceans
          (Abdelrady et al., 2016, Section 4.4-4.5; Fisher et al., 2023).

    Parameters:
        WST_C (Union[Raster, np.ndarray]): Water surface skin temperature in degrees
            Celsius (°C). Required.
        albedo (Union[Raster, np.ndarray], optional): Broadband shortwave surface
            albedo (dimensionless, 0-1). Used to compute `SWnet = SWin_Wm2 * (1 - albedo)`
            when `SWnet` is omitted; retrieved from GEOS-5 FP if needed and
            `offline_mode=False`.
        Ta_C (Union[Raster, np.ndarray], optional): Near-surface air temperature in
            degrees Celsius (°C). Retrieved from GEOS-5 FP if omitted when `Td_C`
            must be estimated and `offline_mode=False`.
        RH (Union[Raster, np.ndarray], optional): Near-surface relative humidity as
            a fraction (0-1). Retrieved from GEOS-5 FP if omitted when `Td_C` must
            be estimated and `offline_mode=False`.
        Td_C (Union[Raster, np.ndarray], optional): Near-surface air dew-point
            temperature in degrees Celsius (°C). Estimated from `Ta_C` and `RH` if
            omitted.
        windspeed_mps (Union[Raster, np.ndarray], optional): Wind speed at reference
            height in meters per second (m s^-1). Retrieved from GEOS-5 FP if
            omitted and `offline_mode=False`.
        SWnet (Union[Raster, np.ndarray], optional): Net shortwave solar radiation
            absorbed by the water surface in Watts per square meter (W m^-2).
            Computed from `SWin_Wm2` (or GEOS-5 FP incoming shortwave) and `albedo`
            if omitted.
        SWin_Wm2 (Union[Raster, np.ndarray], optional): Downwelling shortwave solar
            radiation in Watts per square meter (W m^-2). Used with `albedo` to
            compute `SWnet` when `SWnet` is omitted; retrieved from GEOS-5 FP if
            omitted and `offline_mode=False`.
        geometry (RasterGeometry, optional): Spatial grid definition (`rasters.RasterGeometry`).
            Inferred automatically from `WST_C.geometry` if `WST_C` is a `Raster`.
        time_UTC (datetime, optional): UTC timestamp of the observation, required
            when retrieving missing meteorological/radiative inputs from GEOS-5 FP.
        water (Union[Raster, np.ndarray], optional): Boolean or binary open-water
            mask (True/1 = water, False/0 = land). Retrieved from NASADEM Surface
            Water Body (`NASADEM.swb`) when `mask_non_water_pixels=True`,
            `geometry` is provided, and `offline_mode=False`.
        GEOS5FP_connection (GEOS5FP, optional): Existing `GEOS5FP` assimilation
            client instance. Created lazily only when a network retrieval is needed.
        resampling (str, optional): Spatial resampling method when reprojecting
            GEOS-5 FP grids onto `geometry`. Defaults to `RESAMPLING_METHOD` ("cubic").
        mask_non_water_pixels (bool, optional): If True, masks non-water pixels in
            `WST_C` to NaN using the `water` mask so land pixels are excluded from
            ETM outputs. If no finite water pixels remain after masking, returns
            geometry-matched NaN outputs without raising `BlankOutputError`.
            Defaults to `MASK_NON_WATER_PIXELS` (True).
        dewpoint_method (str, optional): Method for estimating `Td_C` from `Ta_C`
            and `RH` when `Td_C` is omitted: `"linear"` (default, Collection 3
            baseline) or `"magnus_tetens"`. Defaults to `DEWPOINT_METHOD` ("linear").
        offline_mode (bool, optional): If True, disables network retrievals from
            GEOS-5 FP and NASADEM, raising `MissingOfflineParameter` if required
            inputs are missing. Defaults to False.

    Returns:
        Dict[str, Union[Raster, np.ndarray, float]]: Dictionary containing masked
        inputs, intermediate ETM parameters, diagnostics, and water heat flux:
            - `WST_C`: Water surface temperature (°C) (masked to water if enabled).
            - `Ta_C`: Air temperature (°C) (or None if `Td_C` was given directly).
            - `RH`: Relative humidity fraction (0-1) (or None if `Td_C` was given).
            - `Td_C`: Dew-point temperature (°C).
            - `windspeed_mps`: Wind speed (m s^-1).
            - `high_wind_flag`: Boolean mask where `windspeed_mps > 7.5` m s^-1.
            - `SWnet`: Net shortwave radiation (W m^-2).
            - `Tn`: Half temperature gradient `0.5 * (WST_C - Td_C)` (°C).
            - `η`: Dimensionless evaporation efficiency.
            - `S`: Scaled wind function `3.3 * windspeed_mps` (m s^-1).
            - `β`: Thermal exchange coefficient (W m^-2 °C^-1).
            - `Te`: Equilibrium water surface temperature (°C).
            - `W_Wm2`: Instantaneous water heat flux (W m^-2).

    Raises:
        ValueError: If `WST_C`, `Td_C`, `windspeed_mps`, or `SWnet` cannot be
            resolved when `offline_mode=False`, or if `dewpoint_method` is invalid.
        MissingOfflineParameter: If `offline_mode=True` and required inputs
            (`water` when masking, `Ta_C`, `RH`, `windspeed_mps`, `SWin_Wm2`,
            or `albedo`) are omitted.

    Examples:
        >>> import numpy as np
        >>> from AquaSEBS import water_heat_flux
        >>> out = water_heat_flux(
        ...     WST_C=np.array([20.0]),
        ...     Td_C=np.array([10.0]),
        ...     windspeed_mps=np.array([3.0]),
        ...     SWnet=np.array([500.0]),
        ...     mask_non_water_pixels=False
        ... )
        >>> round(float(out["W_Wm2"][0]), 2)
        331.15
    """
    def _get_geos5fp() -> GEOS5FP:
        nonlocal GEOS5FP_connection
        if GEOS5FP_connection is None:
            GEOS5FP_connection = GEOS5FP()
        return GEOS5FP_connection

    if dewpoint_method not in ("linear", "magnus_tetens"):
        raise ValueError(
            f"unsupported dewpoint_method: {dewpoint_method!r} (expected 'linear' or 'magnus_tetens')"
        )

    # Infer spatial grid geometry from the thermal water surface temperature Raster
    # when an explicit RasterGeometry is not supplied
    if geometry is None and isinstance(WST_C, Raster):
        geometry = WST_C.geometry

    # Water surface skin temperature (T0 / WST_C) is the primary remote-sensing state
    # variable and dominant sensitivity driver of the Equilibrium Temperature Model
    # (Abdelrady et al., 2016, Section 4.6)
    if WST_C is None:
        raise ValueError("water surface temperature (WST_C) not given")

    check_distribution(WST_C, "WST_C")

    # Retrieve static open-water extent from the NASADEM Surface Water Body (SWB)
    # dataset to prevent terrestrial land thermal signatures from contaminating
    # open-water heat flux estimates (Fisher et al., 2023, Methods)
    if mask_non_water_pixels:
        if water is None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but water mask is not provided")
            if geometry is not None:
                water = NASADEM.swb(geometry=geometry)
                check_distribution(water, "water")
            else:
                raise ValueError("water mask not given and geometry not available for NASADEM.swb lookup")

        # Apply open-water mask by setting non-water pixels to NaN; NaNs propagate
        # cleanly through all downstream ETM equations
        WST_C = rt.where(water, WST_C, np.nan)

    # Gracefully handle land-only tiles or scenes where all open-water pixels are
    # cloud-obscured (zero finite water pixels remain after masking) so callers
    # such as JET3 do not abort with BlankOutputError on valid terrestrial scenes
    if not np.any(np.isfinite(np.asarray(WST_C))):
        logger.warning("no valid unmasked water pixels present; returning NaN water heat flux outputs")
        nan_output = _nan_like(WST_C)
        high_wind_flag = (
            windspeed_mps > HIGH_WIND_THRESHOLD_MPS
            if windspeed_mps is not None
            else np.zeros(np.asarray(WST_C).shape, dtype=bool)
        )
        return {
            "WST_C": WST_C,
            "Ta_C": Ta_C,
            "RH": RH,
            "Td_C": Td_C if Td_C is not None else _nan_like(WST_C),
            "windspeed_mps": windspeed_mps if windspeed_mps is not None else _nan_like(WST_C),
            "high_wind_flag": high_wind_flag,
            "SWnet": SWnet if SWnet is not None else _nan_like(WST_C),
            "Tn": _nan_like(WST_C),
            "η": _nan_like(WST_C),
            "S": _nan_like(WST_C),
            "β": _nan_like(WST_C),
            "Te": _nan_like(WST_C),
            "W_Wm2": nan_output,
        }

    # Estimate near-surface air dew-point temperature (Td_C) from air temperature
    # and relative humidity if Td_C is not provided directly; Td_C sets the lower
    # thermodynamic baseline for equilibrium water temperature (Abdelrady et al.,
    # 2016, Eq. 8 & Eq. 12)
    if Td_C is None:
        if Ta_C is None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but Ta_C (or Td_C) is not provided")
            if geometry is not None and time_UTC is not None:
                Ta_C = _get_geos5fp().Ta_C(
                    time_UTC=time_UTC,
                    geometry=geometry,
                    resampling=resampling
                )

        if Ta_C is not None:
            check_distribution(Ta_C, "Ta_C")

        if RH is None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but RH (or Td_C) is not provided")
            if geometry is not None and time_UTC is not None:
                RH = _get_geos5fp().RH(
                    time_UTC=time_UTC,
                    geometry=geometry,
                    resampling=resampling
                )

        if RH is not None:
            check_distribution(RH, "RH")

        if Ta_C is not None and RH is not None:
            if dewpoint_method == "linear":
                # Approximate dew-point temperature (°C) using the standard linear psychrometric
                # rule-of-thumb: Td ≈ Ta - ((100 - RH%) / 5), where RH% = RH * 100
                Td_C = Ta_C - ((100 - RH * 100) / 5.0)
            elif dewpoint_method == "magnus_tetens":
                # WMO / Alduchov & Eskridge (1996) Magnus-Tetens saturation vapor pressure formula
                a = 17.625
                b = 243.04
                rh_clipped = np.clip(RH, 1e-6, 1.0)
                gamma_m = np.log(rh_clipped) + (a * Ta_C) / (b + Ta_C)
                Td_C = (b * gamma_m) / (a - gamma_m)

    if Td_C is None:
        raise ValueError("dew-point temperature (Td_C) not given")

    check_distribution(Td_C, "Td_C")

    # Retrieve reference-height wind speed (u, m s^-1) from GEOS-5 FP if omitted;
    # wind speed controls turbulent mechanical mixing and convective vapor removal
    # across the air-water boundary layer (Abdelrady et al., 2016, Eq. 10;
    # Fisher et al., 2023, Table 2)
    if windspeed_mps is None:
        if offline_mode:
            raise MissingOfflineParameter("offline mode is enabled but windspeed_mps is not provided")
        if geometry is not None and time_UTC is not None:
            geos5fp = _get_geos5fp()
            wind_fn = getattr(geos5fp, "wind_speed_mps", None) or getattr(geos5fp, "wind_speed")
            windspeed_mps = wind_fn(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )

    if windspeed_mps is None:
        raise ValueError("wind-speed (windspeed_mps) not given")

    check_distribution(windspeed_mps, "windspeed_mps")
    high_wind_flag = windspeed_mps > HIGH_WIND_THRESHOLD_MPS

    # Derive net shortwave solar radiation absorbed by the water surface:
    #   SWnet = SWin_Wm2 * (1 - albedo)
    # where SWin_Wm2 (W m^-2) is downwelling shortwave irradiance and albedo (0-1) is
    # broadband water surface reflectance (Abdelrady et al., 2016, Section 2.2.1;
    # Fisher et al., 2023, Methods)
    if SWnet is None:
        if SWin_Wm2 is None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but SWin_Wm2 (or SWnet) is not provided")
            if geometry is not None and time_UTC is not None:
                SWin_Wm2 = _get_geos5fp().SWin(
                    time_UTC=time_UTC,
                    geometry=geometry,
                    resampling=resampling
                )

        if SWin_Wm2 is not None:
            check_distribution(SWin_Wm2, "SWin_Wm2")

        if albedo is None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but albedo (or SWnet) is not provided")
            if geometry is not None and time_UTC is not None:
                albedo = _get_geos5fp().ALBEDO(
                    time_UTC=time_UTC,
                    geometry=geometry,
                    resampling=resampling
                )

        if albedo is not None:
            check_distribution(albedo, "albedo")

        if SWin_Wm2 is not None and albedo is not None:
            SWnet = SWin_Wm2 * (1 - albedo)

    if SWnet is None:
        raise ValueError("net shortwave radiation (SWnet) not given")

    check_distribution(SWnet, "SWnet")

    # Step 1 — Half temperature gradient Tn (°C) (Abdelrady et al., 2016, Eq. 8;
    # Fisher et al., 2023, Eq. 1): characterizes the vertical air-water vapor and
    # thermal gradient driving latent and sensible heat exchange
    Tn = 0.5 * (WST_C - Td_C)
    check_distribution(Tn, "Tn")

    # Step 2 — Dimensionless evaporation efficiency η (Abdelrady et al., 2016, Eq. 9;
    # Fisher et al., 2023, Eq. 2): parameterizes the non-linear increase in saturation
    # vapor pressure slope and evaporative conductance with warmer water temperature
    # (0.015 * WST_C) and strong air-water thermal gradients (0.0012 * Tn^2) above
    # the neutral baseline of 0.35 (Edinger et al., 1968)
    η = 0.35 + 0.015 * WST_C + 0.0012 * (Tn ** 2)
    check_distribution(η, "η")

    # Step 3 — Empirical wind function S(W) (m s^-1) (Abdelrady et al., 2016, Eq. 10;
    # Fisher et al., 2023, Eq. 3): scales reference-height wind speed by 3.3 to
    # represent turbulent aerodynamic transfer efficiency over open water
    S = 3.3 * windspeed_mps
    check_distribution(S, "S")

    # Step 4 — Thermal exchange coefficient β (W m^-2 °C^-1) (Abdelrady et al., 2016,
    # Eq. 11; Fisher et al., 2023, Eq. 4): total linearized surface heat-transfer
    # conductance combining longwave radiative linearization (4.5 + 0.05 * WST_C)
    # with wind-driven turbulent latent and sensible exchange ((η + 0.47) * S)
    β = 4.5 + 0.05 * WST_C + (η + 0.47) * S
    check_distribution(β, "β")

    # Step 5 — Equilibrium temperature Te (°C) (Abdelrady et al., 2016, Eq. 12;
    # Fisher et al., 2023, Eq. 5): hypothetical steady-state water skin temperature
    # at which solar radiative heating (SWnet) balances turbulent and radiative
    # heat losses to the atmosphere (Edinger et al., 1968)
    Te = Td_C + (SWnet / β)
    check_distribution(Te, "Te")

    # Step 6 — Instantaneous water heat flux W_Wm2 / G0w (W m^-2) (Abdelrady et al.,
    # 2016, Eq. 13; Fisher et al., 2023, Eq. 6): net heat flux into (> 0) or out of
    # (< 0) the water column proportional to thermal departure from equilibrium (Te - WST_C)
    W_Wm2 = β * (Te - WST_C)
    check_distribution(W_Wm2, "W_Wm2")

    return {
        "WST_C": WST_C,
        "Ta_C": Ta_C,
        "RH": RH,
        "Td_C": Td_C,
        "windspeed_mps": windspeed_mps,
        "high_wind_flag": high_wind_flag,
        "SWnet": SWnet,
        "Tn": Tn,
        "η": η,
        "S": S,
        "β": β,
        "Te": Te,
        "W_Wm2": W_Wm2
    }
