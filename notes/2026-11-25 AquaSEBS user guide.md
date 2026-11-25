# AquaSEBS User Guide

**Package:** `AquaSEBS` (v1.3.0)  
**Repository:** [JPL-Evapotranspiration-Algorithms/AquaSEBS](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS)  
**Primary Literature References:**
1. **Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).** *Surface Energy Balance of Fresh and Saline Waters: AquaSEBS*. Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
2. **Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W., Pearson, C., & Huntington, J. L. (2023).** *Remotely sensed terrestrial open water evaporation*. Scientific Reports, 13, 8174. https://doi.org/10.1038/s41598-023-34921-2

---

## 1. Introduction

**AquaSEBS** (Aquatic Surface Energy Balance System) is a Python package for estimating instantaneous and daily open-water evaporation (latent heat flux) and surface energy balance components over lakes, reservoirs, inland seas, and coastal waters using thermal infrared satellite observations and meteorological inputs.

The package implements the operational hybrid formulation validated across 19 global open-water sites in **Fisher et al. (2023)** and deployed in **NASA's ECOSTRESS Collection 2+** (`L3_ET_PT-JPL`) and the **OpenET** platform:
* **Water Heat Flux ($W$):** Computed via the **Equilibrium Temperature Model (ETM)** of Edinger et al. (1968) as adapted for remote sensing by Abdelrady et al. (2016, Eqs. 8–13).
* **Latent Heat Flux ($LE$):** Computed via the energy-limited **Priestley–Taylor (1972)** equation ($\alpha = 1.26$) applied to net available energy $(R_n - W)$.
* **Automated Data & Masking Integration:** Automatically retrieves missing meteorological and shortwave radiative forcings from NASA's **GEOS-5 FP** assimilation system and masks non-water pixels using the **NASADEM Surface Water Body (`NASADEM.swb`)** product.

---

## 2. Installation

### 2.1 Install from PyPI

```bash
pip install AquaSEBS
```

### 2.2 Install from Source (Development Environment)

```bash
git clone https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS.git
cd AquaSEBS

# Create a Conda/Mamba environment (Python 3.11+ recommended)
mamba create -y -n AquaSEBS -c conda-forge python=3.13
mamba activate AquaSEBS

# Install in editable mode with test and build dependencies
pip install -e .[dev]
```

### 2.3 Core Dependencies

* `numpy` — Array math and elementwise operations
* `rasters` — Geospatial `Raster` and `RasterGeometry` handling
* `GEOS5FP` — Automated NASA GEOS-5 FP meteorological and radiation retrieval
* `NASADEM` — Automated static surface water body (`swb`) masking
* `priestley-taylor` ($\ge 1.1.2$) — Psychrometric slope ($\Delta / (\Delta + \gamma)$) and Priestley–Taylor flux solver
* `verma-net-radiation` — Broadband net radiation ($R_n$) calculation
* `daylight-evapotranspiration` — Sinusoidal diurnal integration from instantaneous $LE$ to daylight totals
* `check-distribution` — Input/output numerical range and distribution logging

---

## 3. Core Concepts & Energy Balance Equations

At the water–atmosphere interface, energy partitioning follows:

$$R_n = LE + H + W$$

where:
* $R_n$ ($\text{W}\cdot\text{m}^{-2}$) is broadband net radiation at the water surface.
* $LE$ ($\text{W}\cdot\text{m}^{-2}$) is latent heat flux consumed by evaporation.
* $W$ ($\text{W}\cdot\text{m}^{-2}$, denoted $G_{0w}$ in literature) is water heat storage/exchange flux ($W > 0$ indicates warming/heat storage; $W < 0$ indicates cooling/heat release).
* $H$ ($\text{W}\cdot\text{m}^{-2}$) is sensible heat flux (implicit residual $R_n - LE - W$).

### 3.1 Water Heat Flux (`water_heat_flux`)

Using water surface skin temperature ($T_0 = \text{WST}_C$, $^\circ\text{C}$), near-surface dew-point temperature ($T_d$, $^\circ\text{C}$), wind speed ($u$, $\text{m}\cdot\text{s}^{-1}$), and net absorbed shortwave radiation ($\text{SW}_{\text{net}}$, $\text{W}\cdot\text{m}^{-2}$), AquaSEBS solves the six Equilibrium Temperature Model equations:

1. **Half Temperature Gradient ($T_n$, $^\circ\text{C}$):**
   $$T_n = 0.5\,(\text{WST}_C - T_d)$$
2. **Evaporation Efficiency ($\eta$, dimensionless):**
   $$\eta = 0.35 + 0.015\,\text{WST}_C + 0.0012\,T_n^2$$
3. **Wind Function ($S$, $\text{m}\cdot\text{s}^{-1}$):**
   $$S = 3.3\,u$$
4. **Thermal Exchange Coefficient ($\beta$, $\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$):**
   $$\beta = 4.5 + 0.05\,\text{WST}_C + (\eta + 0.47)\,S$$
5. **Equilibrium Temperature ($T_e$, $^\circ\text{C}$):**
   $$T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}$$
6. **Instantaneous Water Heat Flux ($W_{\text{Wm2}}$, $\text{W}\cdot\text{m}^{-2}$):**
   $$W_{\text{Wm2}} = \beta\,(T_e - \text{WST}_C)$$

### 3.2 Latent Heat Flux (`AquaSEBS`)

Given net radiation $R_n$, water heat flux $W$, air temperature $T_a$ ($^\circ\text{C}$), and psychrometric constant $\gamma$ ($\text{Pa}\cdot\text{K}^{-1}$), instantaneous open-water latent heat flux is computed as:

$$LE_{\text{Wm2}} = \alpha \left(\frac{\Delta}{\Delta + \gamma}\right)(R_n - W_{\text{Wm2}})$$

where $\alpha = 1.26$ (`PT_ALPHA`) is the Priestley–Taylor coefficient over open water and $\epsilon = \Delta / (\Delta + \gamma)$ is the dimensionless psychrometric weighting factor.

---

## 4. Usage Workflows

### 4.1 Standalone Offline Execution (NumPy Arrays / Scalars)

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

### 4.2 Automated Satellite Raster Workflow (with GEOS-5 FP & NASADEM)

When passing a georeferenced `rasters.Raster` for `WST_C` (such as ECOSTRESS L2 LSTE or Landsat ARD Surface Temperature) and an observation timestamp `time_UTC`, `AquaSEBS()` automatically:
1. Infers `geometry` from `WST_C.geometry`.
2. Downloads the `NASADEM.swb` open-water mask and masks non-water pixels to `NaN`.
3. Retrieves missing `Ta_C`, `RH`, `windspeed_mps`, and `SWin_Wm2` from NASA `GEOS5FP`.
4. Computes `Rn_Wm2` via `verma_net_radiation()` (when `albedo` and `emissivity` are provided).

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

### 4.3 Computing Water Heat Flux Directly (`water_heat_flux`)

If you only need the Equilibrium Temperature Model water heat flux and intermediate heat-exchange parameters:

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

### 4.4 Daylight / Daily Evapotranspiration Upscaling

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

---

## 5. API Reference

### 5.1 `AquaSEBS()`

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
    geometry=None,
    time_UTC=None,
    water=None,
    GEOS5FP_connection=None,
    α=PT_ALPHA,
    γ_Pa=GAMMA_PA,
    resampling="cubic",
    upscale_to_daylight=False,
    mask_non_water_pixels=True,
    offline_mode=False
) -> Dict[str, Union[Raster, np.ndarray, float]]
```

#### Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `WST_C` | `Raster \| ndarray` | **Required** | Water surface skin temperature ($^\circ\text{C}$). |
| `emissivity` | `Raster \| ndarray \| None` | `None` | Surface thermal emissivity ($0\text{–}1$). Required if `Rn_Wm2` is omitted. |
| `albedo` | `Raster \| ndarray \| None` | `None` | Broadband surface albedo ($0\text{–}1$). Required if `Rn_Wm2` is omitted. |
| `Ta_C` | `Raster \| ndarray \| None` | `None` | Near-surface air temperature ($^\circ\text{C}$). Auto-fetched from `GEOS5FP` if omitted. |
| `RH` | `Raster \| ndarray \| None` | `None` | Near-surface relative humidity fraction ($0\text{–}1$). Auto-fetched if omitted. |
| `Td_C` | `Raster \| ndarray \| None` | `None` | Near-surface dew-point temperature ($^\circ\text{C}$). Estimated from `Ta_C` and `RH` if omitted. |
| `windspeed_mps` | `Raster \| ndarray \| None` | `None` | Wind speed ($\text{m}\cdot\text{s}^{-1}$). Auto-fetched from `GEOS5FP` if omitted. |
| `SWnet` | `Raster \| ndarray \| None` | `None` | Net shortwave radiation ($\text{W}\cdot\text{m}^{-2}$). Computed from `SWin` and `albedo` if omitted. |
| `Rn_Wm2` | `Raster \| ndarray \| None` | `None` | Broadband net radiation ($\text{W}\cdot\text{m}^{-2}$). Computed via `verma_net_radiation` if omitted. |
| `W_Wm2` | `Raster \| ndarray \| None` | `None` | Water heat flux ($\text{W}\cdot\text{m}^{-2}$). Computed via `water_heat_flux` if omitted. |
| `SWin_Wm2` | `Raster \| ndarray \| None` | `None` | Incoming shortwave radiation ($\text{W}\cdot\text{m}^{-2}$). Auto-fetched if omitted. |
| `geometry` | `RasterGeometry \| None` | `None` | Spatial grid definition; inferred from `WST_C.geometry` if `WST_C` is a `Raster`. |
| `time_UTC` | `datetime \| None` | `None` | Observation UTC timestamp for reanalysis queries and daylight integration. |
| `water` | `Raster \| ndarray \| None` | `None` | Binary/boolean water mask (`True`/`1` = water). Auto-fetched via `NASADEM.swb` if `geometry` is given. |
| `GEOS5FP_connection` | `GEOS5FP \| None` | `None` | Reusable `GEOS5FP` connection object. |
| `α` | `Raster \| ndarray \| float` | `1.26` | Priestley–Taylor coefficient (dimensionless). |
| `γ_Pa` | `Raster \| ndarray \| float` | `GAMMA_PA` | Psychrometric constant ($\text{Pa}\cdot\text{K}^{-1}$). |
| `resampling` | `str` | `"cubic"` | Spatial resampling method for reprojecting meteorological grids. |
| `upscale_to_daylight` | `bool` | `False` | Enable diurnal integration to daylight/daily ET totals. |
| `mask_non_water_pixels` | `bool` | `True` | Mask non-water pixels in `WST_C` to `NaN` using `water`. |
| `offline_mode` | `bool` | `False` | Raise `MissingOfflineParameter` rather than querying `GEOS5FP` when inputs are absent. |

#### Return Dictionary Keys

* **Primary Fluxes:** `LE_Wm2`, `LE_potential_Wm2`, `W_Wm2`, `Rn_Wm2` ($\text{W}\cdot\text{m}^{-2}$)
* **ETM Intermediates (when `W_Wm2` is calculated):** `Tn` ($^\circ\text{C}$), `η` (dimensionless), `S` ($\text{m}\cdot\text{s}^{-1}$), `β` ($\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$), `Te` ($^\circ\text{C}$)
* **Psychrometric & Radiation Intermediates:** `epsilon`, plus `verma_net_radiation` outputs (`SWnet`, `LWin`, `LWout`, etc.) when `Rn_Wm2` is calculated internally
* **Daylight Upscaling Outputs (when `upscale_to_daylight=True`):** `daylight_evapotranspiration` outputs (`LE_daylight_Wm2`, `ET_daylight_mm`, etc.)

---

## 6. Operational Best Practices & Scientific Caveats

1. **Thermal Calibration & Cloud Masking:** Water heat flux $W$ is most sensitive to water surface skin temperature (`WST_C`). Ensure input thermal imagery has radiometric accuracy better than $\pm 0.3\text{--}0.5\text{ K}$ and strict cloud/cloud-shadow screening, as cold cloud edges create artificial temperature gradients $(T_e - T_0)$.
2. **High-Wind Events ($u > 7.5\text{ m}\cdot\text{s}^{-1}$):** Fisher et al. (2023) demonstrated that during strong wind events ($u > 7.5\text{ m}\cdot\text{s}^{-1}$ daily mean), open-water evaporation transitions from radiatively controlled to aerodynamically controlled, causing the Priestley–Taylor formulation to underestimate instantaneous fluxes ($r^2$ drops from $0.71$ to $0.47$). Flag or filter high-wind instantaneous scenes, or rely on daily integrated totals where wind errors average down ($\text{RMSE} \approx 1.2\text{--}1.5\text{ mm}\cdot\text{day}^{-1}$).
3. **Saline and Hypersaline Water Bodies:** Version `1.3.0` computes freshwater potential evaporation ($S = 0$). Over marine ($S \approx 34.7\text{ g}\cdot\text{L}^{-1}$, $\sim 1\%$ reduction) or hypersaline water bodies such as the Great Salt Lake ($S \approx 240\text{ g}\cdot\text{L}^{-1}$, $\sim 27\%$ reduction), apply the Turk (1970) / Abdelrady et al. (2016, Eq. 19) salinity reduction factor post-hoc:
   $$\sigma = 1.025 - 0.0246 \exp(0.00879 \cdot S), \qquad LE_{\text{saline}} = \sigma \cdot LE_{\text{fresh}}$$
4. **Arid Environments ($\text{RH} < 50\%$):** When `Td_C` is omitted, `water_heat_flux()` estimates dew point via $T_d = T_a - (100 - 100 \cdot \text{RH})/5.0$. In arid reservoirs with low relative humidity, supply an explicit psychrometric `Td_C` array/raster for highest accuracy.
5. **Dynamic Water Masks for Fluctuating Reservoirs:** Static water masks (`NASADEM.swb`) may include exposed shoreline sediments during severe reservoir drawdown. For fluctuating reservoirs, supply a scene-specific `water` mask derived from concurrent optical/SWIR water indices (e.g., NDWI/MNDWI) and high-resolution thermal data ($30\text{ m}$ Landsat or $70\text{ m}$ ECOSTRESS).