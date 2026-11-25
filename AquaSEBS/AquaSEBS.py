"""
Aquatic Surface Energy Balance System (AquaSEBS) for Open-Water Evaporation.

This module implements the remote-sensing open-water surface energy balance and
evaporation algorithm combining the Equilibrium Temperature Model (ETM) for water
heat flux (Abdelrady et al., 2016) with the energy-limited Priestley-Taylor
evaporation formulation used operationally in NASA's ECOSTRESS and OpenET
systems (Fisher et al., 2023).

Primary References:
    - Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).
      Surface Energy Balance of Fresh and Saline Waters: AquaSEBS.
      Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
    - Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W.,
      Pearson, C., & Huntington, J. L. (2023). Remotely sensed terrestrial open
      water evaporation. Scientific Reports, 13, 8174.
      https://doi.org/10.1038/s41598-023-34921-2
    - Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface
      heat flux and evaporation using large-scale parameters. Monthly Weather
      Review, 100(2), 81-92.
    - Edinger, J. E., Duttweiler, D. W., & Geyer, J. C. (1968). The response of
      water temperatures to meteorological conditions. Water Resources Research,
      4(5), 1137-1143.
    - Verma, M., Fisher, J. B., Mallick, K., Ryu, Y., Kobayashi, H., Guillaume, A.,
      Moore, G., Ramakrishnan, L., Hendrix, V., Wolf, S., Sikka, M., Kiely, G.,
      Wohlfahrt, G., Gielen, B., Roupsard, O., Toscano, P., Arain, M. A., &
      Cescatti, A. (2016). Global surface net-radiation at 5 km from MODIS Terra.
      Remote Sensing, 8(9), 739. https://doi.org/10.3390/rs8090739
    - Bisht, G., Venturini, V., Islam, S., & Jiang, L. (2005). Estimation of the
      net radiation using MODIS (Moderate Resolution Imaging Spectroradiometer).
      Remote Sensing of Environment, 97(1), 52-67.
"""

from typing import Union, Dict
import numpy as np
from datetime import datetime
import logging
from pytictoc import TicToc
import rasters as rt
from rasters import Raster, RasterGeometry
from NASADEM import NASADEM
from GEOS5FP import GEOS5FP
from check_distribution import check_distribution
from priestley_taylor import epsilon_from_Ta_C, GAMMA_PA, PT_ALPHA
from verma_net_radiation import verma_net_radiation
from priestley_taylor import priestley_taylor
from daylight_evapotranspiration import daylight_ET_from_instantaneous_LE

from .constants import *
from .exceptions import *

from .water_heat_flux import water_heat_flux

logger = logging.getLogger(__name__)


def AquaSEBS(
        WST_C: Union[Raster, np.ndarray],
        emissivity: Union[Raster, np.ndarray] = None,
        albedo: Union[Raster, np.ndarray] = None,
        Ta_C: Union[Raster, np.ndarray] = None,
        RH: Union[Raster, np.ndarray] = None,
        Td_C: Union[Raster, np.ndarray] = None,
        windspeed_mps: Union[Raster, np.ndarray] = None,
        SWnet: Union[Raster, np.ndarray] = None,
        Rn_Wm2: Union[Raster, np.ndarray] = None,
        W_Wm2: Union[Raster, np.ndarray] = None,
        SWin_Wm2: Union[Raster, np.ndarray] = None,
        geometry: RasterGeometry = None,
        time_UTC: datetime = None,
        water: Union[Raster, np.ndarray] = None,
        GEOS5FP_connection: GEOS5FP = None,
        α: Union[Raster, np.ndarray, float] = PT_ALPHA,
        γ_Pa: Union[Raster, np.ndarray, float] = GAMMA_PA,
        resampling: str = RESAMPLING_METHOD,
        upscale_to_daylight: bool = UPSCALE_TO_DAYLIGHT,
        mask_non_water_pixels: bool = MASK_NON_WATER_PIXELS,
        offline_mode: bool = False) -> Dict[str, Union[Raster, np.ndarray, float]]:
    """
    Estimate instantaneous and daylight open-water latent heat flux (evaporation)
    and surface energy balance components using the AquaSEBS model.

    Scientific Methodology:
        Over open water bodies, surface energy partitioning follows:

            Rn = LE + H + W

        where Rn is net radiation (W m^-2), LE is latent heat flux (W m^-2),
        H is sensible heat flux (W m^-2), and W is water heat storage/exchange
        flux (W m^-2) (Abdelrady et al., 2016, Eq. 1). Unlike terrestrial land
        surfaces where daily ground heat flux is often assumed negligible, water
        bodies absorb solar radiation at depth and store/release substantial heat
        through vertical mixing and thermal inertia (Abdelrady et al., 2016).

        This implementation couples two validated parameterizations:
        1. **Water Heat Flux (W)** via the Equilibrium Temperature Model (ETM) of
           Edinger et al. (1968) as adapted for remote sensing by Abdelrady et al.
           (2016, Eqs. 8-13) and Fisher et al. (2023, Eqs. 1-6). ETM computes the
           hypothetical equilibrium water temperature (Te) at which net air-water
           heat exchange vanishes and derives W = beta * (Te - WST_C).
        2. **Latent Heat Flux (LE)** via the energy-limited Priestley-Taylor
           (1972) formulation adopted in ECOSTRESS Collection 2 and OpenET
           (Fisher et al., 2023, Eq. 7):

               LE = alpha * (Delta / (Delta + gamma)) * (Rn - W)

           with default Priestley-Taylor coefficient alpha = 1.26 over open water
           where surface resistance is zero.
        3. **Daylight Upscaling** (optional) via sinusoidal diurnal net radiation
           integration and constant evaporative fraction preservation between
           satellite overpass and daily daylight totals (Bisht et al., 2005;
           Fisher, 2018; Fisher et al., 2023).

    Scientific Caveats & Validation Benchmarks:
        - **High-Wind Regime Shift**: Across 19 global validation sites, Fisher
          et al. (2023) found that instantaneous AquaSEBS accuracy degrades during
          high-wind events (u > 7.5 m s^-1 daily mean; r^2 drops from 0.71 to 0.47,
          RMSE rises from 38% to 62% of mean) as evaporation shifts from radiatively
          controlled to atmospherically/aerodynamically controlled. Temporal
          integration to daily totals buffers this sensitivity (daily RMSE =
          1.2-1.5 mm day^-1).
        - **Salinity**: This function computes freshwater potential evaporation.
          For saline or hypersaline waters, evaporation is reduced by Raoult's law
          vapor-pressure depression (up to ~27% in the Great Salt Lake; Abdelrady
          et al., 2016, Eqs. 19-20).
        - **Thermal Sensitivity**: Water heat flux W is the primary source of
          propagated uncertainty and is most sensitive to water surface temperature
          WST_C (Abdelrady et al., 2016, Section 4.6).

    Parameters:
        WST_C (Union[Raster, np.ndarray]): Water surface skin temperature in degrees
            Celsius (°C), typically derived from thermal infrared remote sensing
            (e.g., ECOSTRESS, Landsat, MODIS). Required.
        emissivity (Union[Raster, np.ndarray], optional): Broadband thermal surface
            emissivity (dimensionless, 0-1). Required if `Rn_Wm2` is not provided.
        albedo (Union[Raster, np.ndarray], optional): Broadband shortwave surface
            albedo (dimensionless, 0-1; typical open water values 0.05-0.15).
            Required if `Rn_Wm2` is not provided.
        Ta_C (Union[Raster, np.ndarray], optional): Near-surface air temperature in
            degrees Celsius (°C). Automatically retrieved from GEOS-5 FP if omitted
            and `geometry` + `time_UTC` are provided.
        RH (Union[Raster, np.ndarray], optional): Near-surface relative humidity as
            a fraction (0-1). Automatically retrieved from GEOS-5 FP if omitted.
        Td_C (Union[Raster, np.ndarray], optional): Near-surface dew-point
            temperature in degrees Celsius (°C). Estimated from `Ta_C` and `RH` if
            omitted.
        windspeed_mps (Union[Raster, np.ndarray], optional): Wind speed at reference
            height in meters per second (m s^-1). Automatically retrieved from
            GEOS-5 FP if omitted.
        SWnet (Union[Raster, np.ndarray], optional): Net shortwave radiation absorbed
            by the water surface in Watts per square meter (W m^-2). Computed from
            incoming shortwave and albedo if omitted.
        Rn_Wm2 (Union[Raster, np.ndarray], optional): Instantaneous broadband net
            radiation at the water surface in Watts per square meter (W m^-2).
            Computed via `verma_net_radiation` if omitted.
        W_Wm2 (Union[Raster, np.ndarray], optional): Instantaneous water heat flux
            in Watts per square meter (W m^-2). Positive values indicate heat storage
            into the water column; negative values indicate heat release to the
            atmosphere. Computed via `water_heat_flux` (ETM) if omitted.
        SWin_Wm2 (Union[Raster, np.ndarray], optional): Incoming solar (shortwave)
            radiation in Watts per square meter (W m^-2). Retrieved from GEOS-5 FP
            if omitted when computing `Rn_Wm2`.
        geometry (RasterGeometry, optional): Spatial grid definition (`rasters.RasterGeometry`).
            Inferred automatically from `WST_C.geometry` when `WST_C` is a `Raster`.
        time_UTC (datetime, optional): UTC timestamp of the observation/overpass,
            used for meteorological reanalysis lookup and solar zenith/daylight
            integration.
        water (Union[Raster, np.ndarray], optional): Boolean or binary open-water
            mask (True/1 = water, False/0 = land). Automatically retrieved from
            NASADEM Surface Water Body (`NASADEM.swb`) when `geometry` is provided.
        GEOS5FP_connection (GEOS5FP, optional): Existing `GEOS5FP` data client
            instance for meteorological forcing queries. Created if None.
        α (Union[Raster, np.ndarray, float], optional): Priestley-Taylor advection-
            aridity coefficient (dimensionless). Defaults to `PT_ALPHA` (1.26)
            for open water surfaces (Priestley & Taylor, 1972; Fisher et al., 2023).
        γ_Pa (Union[Raster, np.ndarray, float], optional): Psychrometric constant in
            Pascals per Kelvin (Pa K^-1). Defaults to `GAMMA_PA` (~66 Pa K^-1).
        resampling (str, optional): Spatial resampling method when reprojecting
            coarse meteorological reanalysis grids onto `geometry`. Defaults to
            `RESAMPLING_METHOD` ("cubic").
        upscale_to_daylight (bool, optional): If True and `time_UTC` is provided,
            integrates instantaneous `LE_Wm2` to daylight/daily evapotranspiration
            using sinusoidal net radiation scaling (Bisht et al., 2005; Fisher et al.,
            2023). Defaults to `UPSCALE_TO_DAYLIGHT` (False).
        mask_non_water_pixels (bool, optional): If True, masks non-water pixels to
            NaN using the `water` mask inside `water_heat_flux`. Defaults to
            `MASK_NON_WATER_PIXELS` (True).
        offline_mode (bool, optional): If True, disables network retrievals from
            GEOS-5 FP and raises `MissingOfflineParameter` if required inputs are
            missing. Defaults to False.

    Returns:
        Dict[str, Union[Raster, np.ndarray, float]]: Dictionary containing instantaneous
        energy balance fluxes, intermediate ETM/Priestley-Taylor variables, and
        (when `upscale_to_daylight=True`) daylight integrated ET products:
            - `LE_Wm2`: Instantaneous latent heat flux (W m^-2).
            - `LE_potential_Wm2`: Potential latent heat flux from Priestley-Taylor (W m^-2).
            - `W_Wm2`: Instantaneous water heat flux (W m^-2).
            - `Rn_Wm2`: Instantaneous net radiation (W m^-2).
            - `epsilon`: Psychrometric energy partitioning factor Delta / (Delta + gamma).
            - `Tn`, `η`, `S`, `β`, `Te`: Intermediate ETM parameters (when `W_Wm2`
              is computed internally).

    Raises:
        ValueError: If required inputs (`WST_C` or `Rn_Wm2` when radiation inputs are
            insufficient) are missing.
        MissingOfflineParameter: If `offline_mode=True` and required meteorological
            or radiative forcings (`Ta_C`, `SWin_Wm2`, etc.) are not supplied.

    Examples:
        >>> import numpy as np
        >>> from datetime import datetime
        >>> from AquaSEBS import AquaSEBS
        >>> results = AquaSEBS(
        ...     WST_C=np.array([[22.0, 23.5]]),
        ...     Ta_C=np.array([[25.0, 25.0]]),
        ...     Td_C=np.array([[14.0, 14.0]]),
        ...     windspeed_mps=np.array([[3.5, 4.0]]),
        ...     SWnet=np.array([[550.0, 560.0]]),
        ...     Rn_Wm2=np.array([[500.0, 510.0]]),
        ...     mask_non_water_pixels=False,
        ...     offline_mode=True
        ... )
        >>> results["LE_Wm2"]
    """
    # Initialize dictionary to accumulate all radiation, water heat flux, latent heat,
    # and optional daylight-upscaled output products
    results = {}

    # Infer spatial grid geometry directly from the thermal surface temperature Raster
    # when an explicit RasterGeometry is not provided
    if geometry is None and isinstance(WST_C, Raster):
        geometry = WST_C.geometry

    # Initialize connection to NASA GEOS-5 FP atmospheric assimilation system
    # for meteorological and radiative forcing retrieval when inputs are omitted
    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    # Water surface skin temperature (T0 / WST_C) is the primary remote-sensing thermal
    # state variable driving both longwave emission and air-water heat exchange
    # (Abdelrady et al., 2016, Section 2.2.1; Fisher et al., 2023)
    if WST_C is None:
        raise ValueError("water surface temperature (WST_C) not given")
    
    check_distribution(WST_C, "WST_C")

    # Retrieve the static open-water mask from the NASADEM Surface Water Body (SWB)
    # product so terrestrial land pixels can be excluded from open-water energy balance
    # calculations (Fisher et al., 2023, Methods)
    if water is None and geometry is not None:
        water = NASADEM.swb(geometry=geometry)
        check_distribution(water, "water")

    # Retrieve near-surface air temperature (Ta_C) from GEOS-5 FP reanalysis if not
    # supplied; Ta_C governs saturation vapor pressure slope (Delta) and dew-point estimation
    if Ta_C is None and geometry is not None and time_UTC is not None:
        if offline_mode:
            raise MissingOfflineParameter("offline mode is enabled but Ta_C is not provided")
        
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    check_distribution(Ta_C, "Ta_C")

    # Compute broadband net radiation Rn = (1 - albedo) * SWin + LWin - LWout
    # (Abdelrady et al., 2016, Eq. 1) using the Verma et al. (2016) parameterization
    # when Rn_Wm2 is not supplied directly
    if Rn_Wm2 is None and albedo is not None and WST_C is not None and emissivity is not None and geometry is not None and time_UTC is not None:
        # Retrieve incoming downwelling shortwave radiation (SWin) from GEOS-5 FP if omitted
        if SWin_Wm2 is None and geometry is not None and time_UTC is not None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but SWin_Wm2 is not provided")
        
            SWin_Wm2 = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )
        
        check_distribution(SWin_Wm2, "SWin_Wm2")

        # Calculate net shortwave, net longwave, and total net radiation (Rn_Wm2)
        # using surface skin temperature (WST_C), surface emissivity, and albedo
        # (Verma et al., 2016)
        Rn_results = verma_net_radiation(
            SWin_Wm2=SWin_Wm2,
            albedo=albedo,
            ST_C=WST_C,
            emissivity=emissivity,
            Ta_C=Ta_C,
            RH=RH,
            geometry=geometry,
            time_UTC=time_UTC,
            resampling=resampling,
            GEOS5FP_connection=GEOS5FP_connection,
            offline_mode=offline_mode
        )

        Rn_Wm2 = Rn_results["Rn_Wm2"]
        results.update(Rn_results)

    if Rn_Wm2 is None:
        raise ValueError("net radiation (Rn_Wm2) not given")

    check_distribution(Rn_Wm2, "Rn_Wm2")

    # Compute instantaneous water heat storage/exchange flux (W_Wm2, denoted G0w in
    # Abdelrady et al., 2016, Section 2.2.1, Eqs. 8-13 and Fisher et al., 2023, Eqs. 1-6)
    # via the Equilibrium Temperature Model (ETM; Edinger et al., 1968). Unlike land
    # surfaces where G0 is small, water bodies absorb radiation at depth and exhibit
    # substantial diurnal heat storage and release; no artificial clipping bounds are
    # imposed on W_Wm2 so negative nighttime/cooling heat release is preserved.
    if W_Wm2 is None:
        water_heat_flux_results = water_heat_flux(
            WST_C=WST_C,
            Ta_C=Ta_C,
            Td_C=Td_C,
            windspeed_mps=windspeed_mps,
            SWnet=SWnet,
            geometry=geometry,
            time_UTC=time_UTC,
            water=water,
            GEOS5FP_connection=GEOS5FP_connection,
            resampling=resampling,
            mask_non_water_pixels=mask_non_water_pixels
        )
        
        W_Wm2 = water_heat_flux_results["W_Wm2"]
        results.update(water_heat_flux_results)
    
    check_distribution(W_Wm2, "W_Wm2")

    # Calculate the dimensionless psychrometric energy-partitioning ratio:
    #   epsilon = Delta / (Delta + gamma)
    # where Delta (Pa K^-1) is the slope of the saturation vapor pressure curve at air
    # temperature Ta_C and gamma (Pa K^-1) is the psychrometric constant
    # (Priestley & Taylor, 1972; Fisher et al., 2023, Eq. 7)
    epsilon = epsilon_from_Ta_C(
        Ta_C=Ta_C,
        gamma_Pa=γ_Pa
    )

    check_distribution(epsilon, "epsilon")

    # Compute instantaneous open-water latent heat flux (evaporation energy, W m^-2)
    # using the Priestley-Taylor (1972) equation with available energy (Rn_Wm2 - W_Wm2)
    # and Priestley-Taylor coefficient alpha = 1.26 (Fisher et al., 2023, Eq. 7):
    #   LE = alpha * (Delta / (Delta + gamma)) * (Rn - G0w)
    # Over open water surfaces without moisture limitation, evaporation is primarily
    # radiation/energy-controlled under calm-to-moderate wind speeds (u <= 7.5 m s^-1).
    priestley_taylor_results = priestley_taylor(
        Rn_Wm2=Rn_Wm2,
        G_Wm2=W_Wm2,
        Ta_C=Ta_C,
        epsilon=epsilon,
        PT_alpha=α
    )

    LE_Wm2 = priestley_taylor_results["LE_potential_Wm2"]
    check_distribution(LE_Wm2, "LE_Wm2")
    results.update(priestley_taylor_results)
    results["LE_Wm2"] = LE_Wm2
    
    # Optionally upscale instantaneous overpass latent heat flux (W m^-2) to daylight
    # integrated evapotranspiration totals by modeling diurnal net radiation sinusoidally
    # from date and solar geometry while preserving the instantaneous evaporative fraction
    # EF = LE / (Rn - W) (Bisht et al., 2005; Abdelrady et al., 2016, Section 4.3;
    # Fisher, 2018; Fisher et al., 2023, Methods)
    if upscale_to_daylight and time_UTC is not None:
        logger.info("started daylight ET upscaling")
        t_et = TicToc()
        t_et.tic()

        # Integrate instantaneous latent heat flux across daylight hours using
        # available energy (Rn_instantaneous_Wm2 - G_instantaneous_Wm2)
        daylight_results = daylight_ET_from_instantaneous_LE(
            LE_instantaneous_Wm2=LE_Wm2,
            Rn_instantaneous_Wm2=Rn_Wm2,
            G_instantaneous_Wm2=W_Wm2,
            time_UTC=time_UTC,
            geometry=geometry
        )
        # Add all returned daylight results to output
        results.update(daylight_results)

        elapsed_et = t_et.tocvalue()
        logger.info(f"completed daylight ET upscaling (elapsed: {elapsed_et:.2f} seconds)")

    return results