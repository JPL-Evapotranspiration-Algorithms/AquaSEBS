[![CI](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS/actions/workflows/ci.yml/badge.svg)](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS/actions/workflows/ci.yml)
[![PyPI version](https://badge.fury.io/py/AquaSEBS.svg)](https://badge.fury.io/py/AquaSEBS)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# AquaSEBS Water Surface Evaporation Python Package

**AquaSEBS** (Aquatic Surface Energy Balance System) is a Python implementation of the Abdelrady et al. (2016) and Fisher et al. (2023) methodology for estimating water surface evaporation and energy balance components over open water bodies using satellite remote sensing data and meteorological inputs. Developed and implemented by **Gregory H. Halverson** at NASA's Jet Propulsion Laboratory for the **ECOSTRESS** mission—and globally validated across 19 open-water sites by **Fisher et al. (2023)**—the package combines Abdelrady et al.'s Equilibrium Temperature Model (ETM) for water heat flux with the energy-limited Priestley–Taylor approach for open-water evaporation.

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 321H<br>

## Citation

If you use `AquaSEBS` in your research, please cite the software package and the foundational scientific publications:

```bibtex
@software{halverson2026aquasebs,
  title={AquaSEBS: Aquatic Surface Energy Balance System Water Surface Evaporation Model},
  author={Halverson, Gregory H.},
  version={1.4.0},
  year={2026},
  publisher={GitHub},
  url={https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS}
}

@article{abdelrady2016aquasebs,
  title={Surface Energy Balance of Fresh and Saline Waters: AquaSEBS},
  author={Abdelrady, Ahmed and Timmermans, Joris and Vekerdy, Zolt{\'a}n and Salama, Mhd. Suhyb},
  journal={Remote Sensing},
  volume={8},
  number={7},
  pages={583},
  year={2016},
  publisher={MDPI},
  doi={10.3390/rs8070583}
}

@article{fisher2023remotely,
  title={Remotely sensed terrestrial open water evaporation},
  author={Fisher, Joshua B and Dohlen, Matthew B and Halverson, Gregory H and Collison, Jacob W and Pearson, Christopher and Huntington, Justin L},
  journal={Scientific Reports},
  volume={13},
  number={1},
  pages={8174},
  year={2023},
  publisher={Nature Publishing Group},
  doi={10.1038/s41598-023-34921-2}
}
```

## Table of Contents

1. [User Guide](#1-user-guide)
2. [Scientific Background and Theoretical Framework](#2-scientific-background-and-theoretical-framework)
3. [Computational Methodology](#3-computational-methodology)
4. [Validation Benchmark Summary](#4-validation-benchmark-summary)
5. [Scientific Caveats and Methodological Limitations](#5-scientific-caveats-and-methodological-limitations)
6. [Repository Release and Commit History](#6-repository-release-and-commit-history)
7. [Contributing, License, and Links](#7-contributing-license-and-links)

## 1. User Guide

### 1.1 Overview

Developed at NASA JPL for **ECOSTRESS Collection 2+** (`L3_ET_PT-JPL` / `L3T_JET`) and deployed in the **OpenET** platform, `AquaSEBS` implements the operational hybrid open-water surface energy balance formulation validated across 19 global open-water sites in **Fisher et al. (2023)**:

- **Water Heat Flux ($W$):** Computed via the **Equilibrium Temperature Model (ETM)** of Edinger et al. (1968) as adapted for remote sensing by Abdelrady et al. (2016, Eqs. 8–13).
- **Latent Heat Flux ($LE$):** Computed via the energy-limited **Priestley–Taylor (1972)** equation ($\alpha = 1.26$) applied to net available energy $(R_n - W)$.
- **Optional Salinity Scaling ($\sigma$):** Supports Turk's (1970) exponential salinity reduction factor via `salinity_gL` (default `None`, preserving freshwater potential evaporation parity with ECOSTRESS Collection 3 `BuildID 0700`).
- **Automated Data & Masking Integration:** Automatically retrieves missing meteorological and shortwave radiative forcings from NASA's **GEOS-5 FP** assimilation system (`wind_speed_mps`, `Ta_C`, `RH`, `SWin`, `ALBEDO`), masks non-water pixels using the **NASADEM Surface Water Body (`NASADEM.swb`)** product, and gracefully returns shape/geometry-matched `NaN` arrays when a scene contains zero valid open-water pixels.

### 1.2 Installation

#### Install from PyPI

```bash
pip install AquaSEBS
```

#### Install from Source (Development Environment)

```bash
git clone https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS.git
cd AquaSEBS

# Create a Conda/Mamba environment (Python 3.11+ supported; 3.14 default)
mamba create -y -n AquaSEBS -c conda-forge python=3.14
mamba activate AquaSEBS

# Install in editable mode with test and build dependencies
pip install -e .[dev]
```

#### Core Dependencies

- `numpy` — Array math and elementwise operations
- `rasters` — Geospatial `Raster` and `RasterGeometry` handling
- `GEOS5FP` — Automated NASA GEOS-5 FP meteorological and radiation retrieval
- `NASADEM` — Automated static surface water body (`swb`) masking
- `priestley-taylor` ($\ge 1.1.2$) — Psychrometric slope ($\Delta / (\Delta + \gamma)$) and Priestley–Taylor flux solver
- `verma-net-radiation` — Broadband net radiation ($R_n$) calculation
- `daylight-evapotranspiration` — Sinusoidal diurnal integration from instantaneous $LE$ to daylight totals
- `check-distribution` — Input/output numerical range and distribution logging
- `pytictoc` — Execution timing for diurnal upscaling

### 1.3 Usage Workflows

#### Standalone Offline Execution (NumPy Arrays / Scalars)

When you already have all meteorological and radiative forcing inputs (from flux towers, weather stations, or pre-loaded arrays), you can run `AquaSEBS()` completely offline without network access:

```python
import numpy as np
from AquaSEBS import AquaSEBS

# Example inputs for two open-water pixels
WST_C = np.array([[22.0, 24.5]])          # Water surface temperature [°C]
Ta_C = np.array([[25.0, 26.0]])           # Air temperature [°C]
Td_C = np.array([[14.0, 15.0]])           # Dew-point temperature [°C]
windspeed_mps = np.array([[3.2, 4.1]])    # Wind speed [m s^-1]
SWnet = np.array([[560.0, 580.0]])        # Net shortwave radiation [W m^-2]
Rn_Wm2 = np.array([[510.0, 530.0]])       # Net broadband radiation [W m^-2]

results = AquaSEBS(
    WST_C=WST_C,
    Ta_C=Ta_C,
    Td_C=Td_C,
    windspeed_mps=windspeed_mps,
    SWnet=SWnet,
    Rn_Wm2=Rn_Wm2,
    mask_non_water_pixels=False,
    offline_mode=True
)

print("Water heat flux W [W m^-2]:", results["W_Wm2"])
print("Latent heat flux LE [W m^-2]:", results["LE_Wm2"])
print("Equilibrium temp Te [°C]:", results["Te"])
```

#### Automated Satellite Raster Workflow (with GEOS-5 FP & NASADEM)

When passing a georeferenced `rasters.Raster` for `WST_C` (such as ECOSTRESS L2 LSTE or Landsat ARD Surface Temperature) and an observation timestamp `time_UTC`, `AquaSEBS()` automatically infers `geometry`, downloads the `NASADEM.swb` open-water mask, retrieves missing forcings from `GEOS5FP`, and computes `Rn_Wm2` via `verma_net_radiation()`:

```python
from datetime import datetime
import rasters as rt
from AquaSEBS import AquaSEBS

# Load satellite thermal, albedo, and emissivity rasters
WST_C = rt.Raster.open("ECOSTRESS_WST_C.tif")
albedo = rt.Raster.open("albedo.tif")
emissivity = rt.Raster.open("emissivity.tif")
time_UTC = datetime(2023, 7, 15, 18, 30)

results = AquaSEBS(
    WST_C=WST_C,
    albedo=albedo,
    emissivity=emissivity,
    time_UTC=time_UTC,
    mask_non_water_pixels=True,
    upscale_to_daylight=True
)

# Save output rasters
results["LE_Wm2"].to_geotiff("AquaSEBS_LE_instantaneous_Wm2.tif")
results["W_Wm2"].to_geotiff("AquaSEBS_W_instantaneous_Wm2.tif")
```

#### Computing Water Heat Flux Directly (`water_heat_flux`)

```python
import numpy as np
from AquaSEBS import water_heat_flux

whf = water_heat_flux(
    WST_C=np.array([20.0]),
    Td_C=np.array([10.0]),
    windspeed_mps=np.array([3.0]),
    SWnet=np.array([500.0]),
    mask_non_water_pixels=False
)

print("Half temp gradient Tn [°C]:", whf["Tn"])
print("Evaporation efficiency η:", whf["η"])
print("Thermal exchange coeff β [W m^-2 °C^-1]:", whf["β"])
print("Equilibrium temp Te [°C]:", whf["Te"])
print("Water heat flux W [W m^-2]:", whf["W_Wm2"])
```

#### Daylight / Daily Evapotranspiration Upscaling

Set `upscale_to_daylight=True` alongside `time_UTC` (and `geometry`) to integrate instantaneous overpass `LE_Wm2` across daylight hours using sinusoidal diurnal net radiation scaling:

```python
results = AquaSEBS(
    WST_C=WST_C,
    Rn_Wm2=Rn_Wm2,
    Ta_C=Ta_C,
    Td_C=Td_C,
    windspeed_mps=windspeed_mps,
    SWnet=SWnet,
    geometry=geometry,
    time_UTC=datetime(2023, 7, 15, 17, 30),
    upscale_to_daylight=True,
    mask_non_water_pixels=False
)
```

### 1.4 API Reference

#### `AquaSEBS()`

```python
AquaSEBS(
    WST_C,
    emissivity=None,
    albedo=None,
    Ta_C=None,
    RH=None,
    Td_C=None,
    windspeed_mps=None,
    SWnet=None,
    Rn_Wm2=None,
    W_Wm2=None,
    SWin_Wm2=None,
    salinity_gL=None,
    geometry=None,
    time_UTC=None,
    water=None,
    GEOS5FP_connection=None,
    α=PT_ALPHA,
    γ_Pa=GAMMA_PA,
    resampling="cubic",
    upscale_to_daylight=False,
    mask_non_water_pixels=True,
    dewpoint_method="linear",
    offline_mode=False
) -> Dict[str, Union[Raster, np.ndarray, float]]
```

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `WST_C` | `Raster \| ndarray` | **Required** | Water surface skin temperature ($^\circ\text{C}$). |
| `emissivity` | `Raster \| ndarray \| None` | `None` | Surface thermal emissivity ($0\text{–}1$). Required if `Rn_Wm2` is omitted. |
| `albedo` | `Raster \| ndarray \| None` | `None` | Broadband surface albedo ($0\text{–}1$). Required if `Rn_Wm2` is omitted; also forwarded to `water_heat_flux` when `SWnet` is omitted. |
| `Ta_C` | `Raster \| ndarray \| None` | `None` | Near-surface air temperature ($^\circ\text{C}$). Auto-fetched from `GEOS5FP` if omitted (when `offline_mode=False`). |
| `RH` | `Raster \| ndarray \| None` | `None` | Near-surface relative humidity fraction ($0\text{–}1$). Auto-fetched if omitted (when `offline_mode=False`) and forwarded to `verma_net_radiation` and `water_heat_flux`. |
| `Td_C` | `Raster \| ndarray \| None` | `None` | Near-surface dew-point temperature ($^\circ\text{C}$). Estimated from `Ta_C` and `RH` if omitted. |
| `windspeed_mps` | `Raster \| ndarray \| None` | `None` | Wind speed ($\text{m}\cdot\text{s}^{-1}$). Auto-fetched from `GEOS5FP` (`wind_speed_mps`) if omitted (when `offline_mode=False`). |
| `SWnet` | `Raster \| ndarray \| None` | `None` | Net shortwave radiation ($\text{W}\cdot\text{m}^{-2}$). Computed from `SWin_Wm2` and `albedo` if omitted. |
| `Rn_Wm2` | `Raster \| ndarray \| None` | `None` | Broadband net radiation ($\text{W}\cdot\text{m}^{-2}$). Computed via `verma_net_radiation` if omitted. |
| `W_Wm2` | `Raster \| ndarray \| None` | `None` | Water heat flux ($\text{W}\cdot\text{m}^{-2}$). Computed via `water_heat_flux` if omitted. |
| `SWin_Wm2` | `Raster \| ndarray \| None` | `None` | Incoming shortwave radiation ($\text{W}\cdot\text{m}^{-2}$). Auto-fetched if omitted (when `offline_mode=False`) and forwarded to `water_heat_flux`. |
| `salinity_gL` | `Raster \| ndarray \| float \| None` | `None` | Optional water salinity ($\text{g}\cdot\text{L}^{-1}$). When `None`, freshwater potential evaporation is computed (Collection 3 baseline); when provided, scales `LE_Wm2` by Turk's (1970) factor $\sigma$. |
| `geometry` | `RasterGeometry \| None` | `None` | Spatial grid definition; inferred from `WST_C.geometry` if `WST_C` is a `Raster`. |
| `time_UTC` | `datetime \| None` | `None` | Observation UTC timestamp for reanalysis queries and daylight integration. |
| `water` | `Raster \| ndarray \| None` | `None` | Binary/boolean water mask (`True`/`1` = water). Auto-fetched via `NASADEM.swb` if `mask_non_water_pixels=True`, `geometry` is given, and `offline_mode=False`. |
| `GEOS5FP_connection` | `GEOS5FP \| None` | `None` | Reusable `GEOS5FP` connection object (instantiated lazily only if a network retrieval is needed). |
| `α` | `Raster \| ndarray \| float` | `1.26` | Priestley–Taylor coefficient (dimensionless). |
| `γ_Pa` | `Raster \| ndarray \| float` | `GAMMA_PA` | Psychrometric constant ($\text{Pa}\cdot\text{K}^{-1}$). |
| `resampling` | `str` | `"cubic"` | Spatial resampling method for reprojecting meteorological grids. |
| `upscale_to_daylight` | `bool` | `False` | Enable diurnal integration to daylight/daily ET totals. |
| `mask_non_water_pixels` | `bool` | `True` | Mask non-water pixels in `WST_C` to `NaN` using `water` (returning all-`NaN` outputs cleanly when zero valid water pixels exist). |
| `dewpoint_method` | `str` | `"linear"` | Dew-point estimation method when `Td_C` is omitted: `"linear"` (Collection 3 baseline) or `"magnus_tetens"`. |
| `offline_mode` | `bool` | `False` | Raise `MissingOfflineParameter` rather than querying `GEOS5FP` or `NASADEM` when inputs are absent. |

#### Return Dictionary Keys

- **Primary Fluxes:** `LE_Wm2`, `LE_potential_Wm2`, `W_Wm2`, `G_Wm2`, `Rn_Wm2` ($\text{W}\cdot\text{m}^{-2}$)
- **ETM Intermediates (when `W_Wm2` is calculated):** `Tn` ($^\circ\text{C}$), `η` (dimensionless), `S` ($\text{m}\cdot\text{s}^{-1}$), `β` ($\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$), `Te` ($^\circ\text{C}$)
- **Diagnostics & Optional Scaling Factors:** `high_wind_flag` (boolean mask where `windspeed_mps > 7.5`), `salinity_reduction_factor` ($\sigma$, included when `salinity_gL` is provided)
- **Psychrometric & Radiation Intermediates:** `epsilon`, plus `verma_net_radiation` outputs (`SWnet`, `LWin`, `LWout`, etc.) when `Rn_Wm2` is calculated internally
- **Daylight Upscaling Outputs (when `upscale_to_daylight=True`):** `daylight_evapotranspiration` outputs (`LE_daylight_Wm2`, `ET_daylight_mm`, etc.)

### 1.5 Development & Release CLI

Common development, testing, and release tasks are available via `makefile`:

```bash
make install                   # Install in editable mode with [dev] dependencies
make test                      # Run pytest suite
make set-version VERSION=X.Y.Z # Update version across pyproject.toml, CITATION.cff, and README.md
make build                     # Build wheel and sdist
make clean                     # Remove build artifacts and caches
```

## 2. Scientific Background and Theoretical Framework

### 2.1 Motivation & Operational Context

Evaporation from open water bodies—spanning oceans, natural lakes, reservoirs, and inland seas—plays a governing role in global and regional hydrological cycles. Oceans account for roughly $86\%$ of global evaporation feeding the atmospheric water cycle, while terrestrial reservoirs and lakes are critical freshwater resources increasingly stressed by climate variability and human abstraction.

Measuring open-water evaporation in situ is logistically difficult and spatially unrepresentative:

- **Instrumentation challenges:** Eddy covariance towers require stable platforms and strict fetch/footprint constraints; Bowen ratio systems lose sub-daily temporal fidelity; and bulk aerodynamic/mass-transfer formulations depend heavily on site-specific empirical wind coefficients.
- **Spatial heterogeneity:** Evaporation rates vary markedly across a single water body due to underlying bathymetry, subsurface heat storage and circulation, shoreline geometry, and horizontal advection of warm, dry air from surrounding land.

Satellite thermal infrared remote sensing resolves spatial gradients in water surface temperature ($T_0$) at fine scales ($30\text{ m}$ Landsat, $70\text{ m}$ ECOSTRESS, $1\text{ km}$ MODIS). **AquaSEBS** adapts terrestrial remote-sensing energy balance models to open water surfaces and serves as the operational open-water evaporation algorithm for **NASA's ECOSTRESS mission** (`L3_ET_PT-JPL` Collection 2+ / `L3T_JET`) and the **OpenET** platform.

### 2.2 Evolution of the AquaSEBS Framework

1. **Original Formulation — Abdelrady et al. (2016):** Adapted Su's (2002) one-source **Surface Energy Balance System (SEBS)** for fresh and saline water bodies by replacing constant-fraction ground heat flux approximations ($G_0 = 0.5 R_n$) with Edinger et al.'s (1968) **Equilibrium Temperature Model (ETM)** for instantaneous water heat flux ($G_{0w}$), calibrating momentum and heat roughness lengths over Lake IJsselmeer ($z_{0m} = 0.0002\text{ m}$, $z_{0h} = 0.0001\text{ m}$, $kB^{-1} \approx 0.3$), defining water-specific wet/dry evaporative fraction bounds, and incorporating Turk's (1970) salinity reduction factor ($\sigma$).
2. **JPL ECOSTRESS Adaptation & Global Validation — Halverson & Fisher et al. (2023):** Adapted by Gregory H. Halverson at NASA JPL for the ECOSTRESS mission and validated globally by Joshua B. Fisher et al. (2023), coupling Abdelrady et al.'s (2016) exact 6-equation Equilibrium Temperature Model (ETM) for instantaneous water heat flux ($W$, denoted $G_{0w}$) with the **Priestley–Taylor (1972)** energy-limited formulation ($\alpha = 1.26$) for instantaneous latent heat flux ($LE$), paired with sinusoidal diurnal radiation scaling (Bisht et al., 2005; Fisher, 2018) for daily totals.

### 2.3 Mathematical Formulation

#### Surface Energy Balance

$$
R_n = LE + H + W
$$

where $R_n$ is net broadband radiation ($\text{W}\cdot\text{m}^{-2}$), $LE$ (or $\lambda E$) is latent heat flux consumed by evaporation ($\text{W}\cdot\text{m}^{-2}$), $H$ is sensible heat flux exchanged with the atmosphere ($\text{W}\cdot\text{m}^{-2}$), and $W$ (or $G_{0w}$) is water heat flux stored in ($W > 0$) or released from ($W < 0$) the water column ($\text{W}\cdot\text{m}^{-2}$).

#### Net Radiation ($R_n$)

When not supplied directly, $R_n$ is computed via the Verma et al. (2016) surface radiation balance:

$$
R_n = \text{SW}_{\text{in}}(1 - \alpha_{\text{surf}}) + \text{LW}_{\text{in}} - \text{LW}_{\text{out}}
$$

#### Water Heat Flux via Equilibrium Temperature Model ($W$)

In situ water heat flux is defined vertically by the depth-integrated heat storage rate $W = \rho_w c_w \int_{0}^{Z} \frac{\Delta T}{\Delta t}\, dZ$. Because thermal infrared satellites observe only skin temperature ($T_0$, upper $\sim 10\text{--}100\,\mu\text{m}$), AquaSEBS solves for the hypothetical **equilibrium temperature** ($T_e$) at which net air–water heat exchange equals zero (Edinger et al., 1968; Abdelrady et al., 2016, Eqs. 8–13; Fisher et al., 2023, Eqs. 1–6):

1. **Half Temperature Gradient ($T_n$, $^\circ\text{C}$):**
   $$T_n = 0.5\,(T_0 - T_d)$$
2. **Evaporation Efficiency Factor ($\eta$, dimensionless):**
   $$\eta = 0.35 + 0.015\,T_0 + 0.0012\,T_n^2$$
3. **Scaled Wind Speed Function ($S(W)$, $\text{m}\cdot\text{s}^{-1}$):**
   $$S(W) = 3.3\,u$$
4. **Thermal Exchange Coefficient ($\beta$, $\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$):**
   $$\beta = 4.5 + 0.05\,T_0 + (\eta + 0.47)\,S(W)$$
5. **Equilibrium Temperature ($T_e$, $^\circ\text{C}$):**
   $$T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}$$
6. **Instantaneous Water Heat Flux ($W$, $\text{W}\cdot\text{m}^{-2}$):**
   $$W = \beta\,(T_e - T_0)$$

#### Latent Heat Flux via Priestley–Taylor ($LE$)

Over open water surfaces where moisture supply is unlimited and surface resistance is zero (Priestley & Taylor, 1972; Fisher et al., 2023, Eq. 7):

$$
LE = \alpha \left(\frac{\Delta}{\Delta + \gamma}\right)(R_n - W)
$$

where $\alpha = 1.26$ is the Priestley–Taylor coefficient, $\Delta$ is the slope of the saturation vapor pressure–temperature curve at air temperature $T_a$, $\gamma$ is the psychrometric constant ($\approx 66\text{ Pa}\cdot^\circ\text{C}^{-1}$ or $0.066\text{ kPa}\cdot^\circ\text{C}^{-1}$), and $\epsilon = \frac{\Delta}{\Delta + \gamma}$ is the dimensionless psychrometric energy-partitioning factor.

#### Salinity Reduction Factor ($\sigma$)

Dissolved salts lower saturation vapor pressure according to Raoult's law. When `salinity_gL` ($S$, $\text{g}\cdot\text{L}^{-1}$) is provided, `AquaSEBS()` applies Turk's (1970) empirical brine evaporation reduction factor (Abdelrady et al., 2016, Eqs. 19–20):

$$
\sigma = \text{clip}\!\left(1.025 - 0.0246 \exp(0.00879 \cdot S),\; 0.0,\; 1.0\right), \qquad LE_s = \sigma \cdot LE_{\text{fresh}}
$$

yielding $< 1\%$ reduction in fresh/brackish water ($0\text{–}20\text{ g}\cdot\text{L}^{-1}$), $\sim 1\%$ in open ocean ($34.7\text{ g}\cdot\text{L}^{-1}$), $\sim 22\%\text{–}27\%$ in hypersaline lakes such as the Great Salt Lake ($240\text{ g}\cdot\text{L}^{-1}$), and $\sim 32\%$ at $300\text{ g}\cdot\text{L}^{-1}$.

## 3. Computational Methodology

### 3.1 Step-by-Step Execution Pipeline

1. **Spatial Geometry & Water Pixel Masking (`AquaSEBS.py`, `water_heat_flux.py`):**
   - Infers `geometry` from `WST_C.geometry` when `WST_C` is a `rasters.Raster` and `geometry` is omitted.
   - When `mask_non_water_pixels=True` (default), retrieves the static binary surface water body extent from `NASADEM.swb(geometry=geometry)` when `water is None`, `geometry` is provided, and `offline_mode=False` (raising `MissingOfflineParameter` if `offline_mode=True`), and sets non-water pixels in `WST_C` to `NaN` via `rt.where(water, WST_C, np.nan)`.
   - If zero finite water pixels remain after masking (e.g., land-only MGRS tiles or cloud-masked water bodies), returns shape- and `RasterGeometry`-matched `NaN` arrays cleanly without raising `BlankOutputError`.
2. **Meteorological & Radiative Forcing Acquisition (`GEOS5FP` & `verma_net_radiation`):**
   - Lazily instantiates `GEOS5FP` only when an omitted forcing must be retrieved and `offline_mode=False`.
   - Retrieves missing `Ta_C`, `RH`, `windspeed_mps` (via `wind_speed_mps`), `SWin`, and `ALBEDO` from NASA **GEOS-5 FP** (default `"cubic"` resampling) when `geometry` and `time_UTC` are provided.
   - Computes `Rn_Wm2` via `verma_net_radiation()` when `Rn_Wm2` is omitted and `albedo`, `emissivity`, `geometry`, and `time_UTC` are supplied.
   - Forwards `RH`, `SWin_Wm2`, `albedo`, `dewpoint_method`, and `offline_mode` from `AquaSEBS()` to `water_heat_flux()`. Estimates dew-point temperature `Td_C` inside `water_heat_flux()` when omitted using either `dewpoint_method="linear"` ($T_d = T_a - \frac{100 - 100 \cdot \text{RH}}{5.0}$, default) or `dewpoint_method="magnus_tetens"`, and computes net shortwave radiation as $\text{SW}_{\text{net}} = \text{SW}_{\text{in}} \cdot (1 - \alpha_{\text{surf}})$.
3. **Water Heat Flux & High-Wind Diagnostics (`water_heat_flux.py`):**
   - Evaluates ETM Equations 8–13 ($T_n, \eta, S, \beta, T_e, W_{\text{Wm2}}$) without artificial clipping bounds so negative nighttime/cooling heat release is preserved, and flags `high_wind_flag = windspeed_mps > 7.5`.
4. **Instantaneous Latent Heat Flux & Salinity Scaling (`AquaSEBS.py`):**
   - Computes $\epsilon = \Delta / (\Delta + \gamma)$ via `priestley_taylor.epsilon_from_Ta_C()` and evaluates potential open-water latent heat flux via `priestley_taylor()`, applying Turk's (1970) factor $\sigma$ if `salinity_gL` is provided. Sensible heat flux ($H$) remains the implicit residual $H = R_n - W - \text{LE}$.
5. **Optional Diurnal Upscaling (`AquaSEBS.py`):**
   - When `upscale_to_daylight=True` and `time_UTC` is provided, invokes `daylight_ET_from_instantaneous_LE()` using instantaneous `LE_Wm2`, `Rn_Wm2`, and `W_Wm2` (as `G_instantaneous_Wm2`) to compute daylight-integrated ET totals.

### 3.2 Literature vs. Codebase Implementation Matrix

| Model Component | Abdelrady et al. (2016) | Fisher et al. (2023) | Implemented Code (`v1.4.0`) |
| :--- | :--- | :--- | :--- |
| **Water Masking** | Manual wind-direction fetch / static mask | MODIS/Landsat land-water mask | Automated `NASADEM.swb` (`rt.where`, all-`NaN` safe) |
| **Net Radiation ($R_n$)** | 4-component radiation input | FLiES/BESS + surface energy balance | `verma_net_radiation` (`GEOS5FP` + thermal) |
| **Water Heat Flux ($W$)** | ETM (Eqs. 8–13) | ETM (Eqs. 1–6) | **Exact match** (`water_heat_flux.py`) |
| **Dew Point ($T_d$)** | Measured / ECMWF reanalysis | MODIS `MOD07_L2` / NCEP reanalysis | `"linear"` default (Col. 3 parity) + `"magnus_tetens"` opt-in |
| **Sensible Heat ($H$)** | Iterative MOS ($z_{0m}=0.0002$, $z_{0h}=0.0001$) | Not solved explicitly | Not solved explicitly |
| **Latent Heat ($\text{LE}$)** | SEBS evaporative fraction ($EF_r$) | Priestley–Taylor ($\alpha = 1.26$) | **Exact match** (`priestley_taylor`) |
| **Salinity Correction ($\sigma$)** | Turk (1970) exponential factor | Noted as future sensitivity | Opt-in `salinity_gL` (`salinity_reduction_factor`) |
| **High-Wind Regime** | — | $u > 7.5\text{ m}\cdot\text{s}^{-1}$ regime shift | Diagnostic `"high_wind_flag"` mask (`u > 7.5`) |
| **Daily Upscaling** | Daily $EF$ scaling | Sinusoidal diurnal $R_n$ + $EF$ | Optional `daylight_ET_from_instantaneous_LE` |

## 4. Validation Benchmark Summary

### 4.1 Validation Datasets (Fisher et al., 2023)

Fisher et al. (2023) evaluated the operational AquaSEBS formulation across **19 global open-water sites** spanning **1986–2019** and **seven Köppen–Geiger climate zones** (`Cfa` $n=6$, `Dfb` $n=5$, `BWh` $n=3$, `BSk` $n=2$, `BSh` $n=1$, `BWk` $n=1$, `Dfa` $n=1$), compiled from the **Great Lakes Evaporation Network (GLEN)**, **USBR Open Water Evaporation Network (OWEN)**, **Zhao & Gao (2019)**, and **Collison (2019)**:

- **In Situ Methods:** Eddy covariance ($n = 11$), Bowen ratio energy balance / BREB ($n = 5$), bulk mass transfer ($n = 2$), and floating evaporation pan ($n = 1$) across half-hourly ($n = 9$), hourly ($n = 1$), daily ($n = 1$), and monthly ($n = 8$) intervals.
- **Sites:** American Falls (ID), Calm (FL), Cochiti (NM), Five-O (FL), Granite Island (MI), Kasumigaura (Japan), Kinneret (Israel), Lahontan (NV), Long Point (ON), Mead (NV/AZ), Mohave (NV/AZ), Nasser (Egypt), Ross Barnett (MS), Spectacle Reef (MI), Stannard Rock (MI), Starr (FL), Taihu (China), White Bear (MN), and White Shoal (MI).
- **Satellite & Reanalysis Inputs:** Landsat 5/7/8 ARD ($30\text{ m}$, $5 \times 5$ pixel windows) + NCEP-NCAR Reanalysis I; Terra/Aqua MODIS (`MOD11A1` $1\text{ km}$ $3 \times 3$ windows, `MCD43A3`, `MOD07_L2`, `MOD04_L2`/`MOD06_L2` FLiES/BESS) + NCEP-NCAR wind speed across **686 cloud-free instantaneous scenes** and **52 cloud-free daily/monthly scenes**.

### 4.2 Quantitative Performance Metrics

#### Fisher et al. (2023) Global Satellite Validation

| Evaluation Regime | Sensor | $r^2$ | RMSE | Bias | Linear Fit |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Instantaneous (High-wind filtered, $u \le 7.5\text{ m}\cdot\text{s}^{-1}$)** | MODIS ($1\text{ km}$) | $0.71$ | $53.7\text{ W}\cdot\text{m}^{-2}$ ($38\%$ of mean) | $-19.14\text{ W}\cdot\text{m}^{-2}$ ($13\%$ of mean) | $y = 0.8x + 3.0$ |
| **Instantaneous (Unfiltered, all wind speeds)** | MODIS ($1\text{ km}$) | $0.47$ | $84.4\text{–}95.4\text{ W}\cdot\text{m}^{-2}$ ($62\%$ of mean) | $-49.5\text{ W}\cdot\text{m}^{-2}$ ($36\%$ of mean) | $y = 0.6x + 13.0$ |
| **Daily Evaporation** | MODIS ($1\text{ km}$) | $0.47$ | $1.5\text{ mm}\cdot\text{day}^{-1}$ ($41\%$ of mean) | $-0.03\text{ to }+0.19\text{ mm}\cdot\text{day}^{-1}$ ($1\%$ of mean) | $y = 0.96x + 0.4$ |
| **Daily Evaporation** | Landsat ($30\text{ m}$) | $0.56$ | $1.2\text{ mm}\cdot\text{day}^{-1}$ ($38\%$ of mean) | $-0.82\text{ mm}\cdot\text{day}^{-1}$ ($26\%$ of mean) | $y = 0.78x + 0.1$ |

#### Abdelrady et al. (2016) Reanalysis & Field Benchmarks

- **Water Heat Flux ($G_{0w}$) vs. ECMWF ERA-Interim:** $\text{rRMSE} = 4.3\%\text{–}8.5\%$ ($\text{RMSE} = 26.5\text{–}68.5\text{ W}\cdot\text{m}^{-2}$) across five seas/oceans ($n \approx 1005$).
- **Freshwater Lakes:** Lake Tana in situ $\text{LE}$ $R^2 = 0.98$, $\text{RMSE} = 35.6\text{ W}\cdot\text{m}^{-2}$ ($\text{rRMSE} = 10.3\%$); Lake Victoria daily $E$ $\text{RMSE} = 1.5\text{ mm}\cdot\text{day}^{-1}$ ($\text{rRMSE} = 7.4\%$).
- **Hypersaline Water:** Great Salt Lake ($S = 240\text{ g}\cdot\text{L}^{-1}$) salinity correction reduced evaporation by $26.5\%$, improving 3-hourly $\text{RMSE}$ from $0.63$ to $0.38\text{ mm}\cdot(3\text{ h})^{-1}$.

### 4.3 Key Findings & Machine Learning Benchmark

- **High-Wind Regime Shift ($u > 7.5\text{ m}\cdot\text{s}^{-1}$):** Instantaneous evaporation shifts from radiatively controlled to aerodynamically controlled during strong winds, causing the Priestley–Taylor formulation to underestimate high-wind peaks ($r^2$ drops from $0.71$ to $0.47$). Daily integration buffers short-term wind spikes ($\text{RMSE} = 1.2\text{–}1.5\text{ mm}\cdot\text{day}^{-1}$).
- **Spatial Resolution Gains:** $30\text{ m}$ Landsat improved daily $r^2$ ($0.56$ vs. $0.47$) and lowered $\text{RMSE}$ ($1.2$ vs. $1.5\text{ mm}\cdot\text{day}^{-1}$) relative to $1\text{ km}$ MODIS while resolving bathymetric and circulation-driven thermal gradients.
- **Machine Learning Benchmark:** None of 11 benchmarked ML models outperformed mechanistic AquaSEBS across the 19 sites. A two-layer TensorFlow Neural Network explained $74\%$ of instantaneous residual variance ($r^2 = 0.74$), with OLS regression confirming **wind speed** ($\text{coef} = +20.10, p < 0.001$) and **relative humidity** ($\text{coef} = -1.56, p < 0.001$) as the primary statistical error predictors, and no significant bias by waterbody area or shape.

## 5. Scientific Caveats and Methodological Limitations

1. **Freshwater Default Without `salinity_gL`:** When `salinity_gL=None` (default, preserving ECOSTRESS Collection 3 `BuildID 0700` consistency), runs assume freshwater ($S = 0$). Supply `salinity_gL` ($\text{g}\cdot\text{L}^{-1}$) over marine or hypersaline waters to apply the Turk (1970) reduction factor $\sigma$.
2. **High-Wind Regime Shift ($u > 7.5\text{ m}\cdot\text{s}^{-1}$):** Instantaneous retrievals during strong wind events ($u > 7.5\text{ m}\cdot\text{s}^{-1}$ daily mean) underestimate aerodynamically driven evaporation. Use the returned `"high_wind_flag"` diagnostic mask to screen or weight instantaneous scenes, or rely on daily integrated totals.
3. **Skin Temperature vs. Bathymetry & Thermal Storage Lag:** ETM estimates instantaneous $W$ from skin temperature ($T_0$) without parameterizing lake depth, stratification, or mixed-layer thickness, which governs diurnal hysteresis in shallow lakes versus deep oceans.
4. **Dew-Point Estimation in Arid Climates:** By default (`dewpoint_method="linear"`), `water_heat_flux()` estimates dew point via $T_d \approx T_a - (100 - 100 \cdot \text{RH}) / 5.0$ for Collection 3 baseline parity. Supply explicit `Td_C` or set `dewpoint_method="magnus_tetens"` in arid environments ($\text{RH} < 50\%$).
5. **Thermal Sensor Accuracy & Cloud Screening:** $W$ is most sensitive to water surface temperature ($T_0$), requiring thermal accuracy within $\pm 0.3\text{--}0.5\text{ K}$ and strict cloud masking.
6. **Shoreline Adjacency, Lateral Advection & Validation Coverage:** Coarse pixels near shorelines or fluctuating reservoir drawdowns can suffer land contamination and lateral dry-air advection; validation sites in Fisher et al. (2023) also lack polar, tropical equatorial, and high-altitude coverage.

## 6. Repository Release and Commit History

| Release Tag | Release Date (UTC) | Release Title | Core Milestone |
| :--- | :--- | :--- | :--- |
| **`v1.0.1`** | 2025-09-30T22:07:17Z | initial release | Initial implementation of Equilibrium Temperature Model (`water_heat_flux`), Verma net radiation + Priestley–Taylor energy balance (`AquaSEBS`), `NASADEM.swb` water masking, and dictionary return schema |
| **`v1.0.2`** | 2025-10-08T17:14:34Z | documentation | Package documentation updates (version bump in `pyproject.toml`) |
| **`v1.1.0`** | 2025-11-11T19:49:53Z | daylight upscaling | Added `upscale_to_daylight` flag, `daylight_ET_from_instantaneous_LE()` integration, and `TicToc` execution timing |
| **`v1.2.0`** | 2026-01-21T21:14:40Z | offline mode | Added `offline_mode` parameter and `MissingOfflineParameter` exception for `Ta_C`, `SWin_Wm2`, and `verma_net_radiation()` |
| **`v1.3.0`** | 2026-01-21T21:55:38Z | `GEOS5FP_INPUTS` constant | Added `GEOS5FP_INPUTS = ["Ta_C"]` to `AquaSEBS/constants.py` |
| **`v1.4.0`** | 2026-11-25 | offline mode, Collection 3 resilience & opt-in science diagnostics | Completed `offline_mode`, lazy `GEOS5FP`, parameter forwarding (`#6`, `#7`), graceful all-`NaN` handling on zero-water scenes (`#14`), `GEOS5FP.wind_speed_mps()` update (`#15`), opt-in `salinity_gL` (`#11`), opt-in `dewpoint_method="magnus_tetens"` (`#12`), diagnostic `"high_wind_flag"` (`#13`), `README.md` LaTeX cleanup (`#16`), and full `pytest` suite (`#8`–`#10`) |

## 7. Contributing, License, and Links

### Contributing

Contributions are welcome via pull requests to [JPL-Evapotranspiration-Algorithms/AquaSEBS](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS). Please ensure `make test` passes before submitting changes.

### License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](LICENSE) file for details.

### Links

- **Homepage:** https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS
- **PyPI:** https://pypi.org/project/AquaSEBS/
- **Issues:** https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS/issues
