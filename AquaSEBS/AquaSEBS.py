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

## TODO use NASADEM surface water body extent to mask out land when processing on rasters

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
        # If geometry is not provided, try to infer from surface temperature raster
    results = {}

    if geometry is None and isinstance(WST_C, Raster):
        geometry = WST_C.geometry

    # Create GEOS5FP connection if not provided
    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    if WST_C is None:
        raise ValueError("water surface temperature (WST_C) not given")
    
    check_distribution(WST_C, "WST_C")

    if water is None and geometry is not None:
        water = NASADEM.swb(geometry=geometry)
        check_distribution(water, "water")

    # Retrieve air temperature if not provided, using GEOS5FP and geometry/time
    if Ta_C is None and geometry is not None and time_UTC is not None:
        if offline_mode:
            raise MissingOfflineParameter("offline mode is enabled but Ta_C is not provided")
        
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    check_distribution(Ta_C, "Ta_C")

    # Compute net radiation if not provided, using albedo, ST_C, and emissivity
    if Rn_Wm2 is None and albedo is not None and WST_C is not None and emissivity is not None and geometry is not None and time_UTC is not None:
        # Retrieve incoming shortwave if not provided
        if SWin_Wm2 is None and geometry is not None and time_UTC is not None:
            if offline_mode:
                raise MissingOfflineParameter("offline mode is enabled but SWin_Wm2 is not provided")
        
            SWin_Wm2 = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )
        
        check_distribution(SWin_Wm2, "SWin_Wm2")

        # Calculate net radiation using Verma et al. method
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

    if W_Wm2 is None:
        # Calculate water heat flux using validated AquaSEBS methodology
        # No artificial constraints applied - trust the physics-based equations
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

    epsilon = epsilon_from_Ta_C(
        Ta_C=Ta_C,
        gamma_Pa=γ_Pa
    )

    check_distribution(epsilon, "epsilon")

    # LE_Wm2 = PT_ALPHA * epsilon * (Rn_Wm2 - W_Wm2)

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
    
    if upscale_to_daylight and time_UTC is not None:
        logger.info("started daylight ET upscaling")
        t_et = TicToc()
        t_et.tic()

        # Use new upscaling function from daylight_evapotranspiration
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
