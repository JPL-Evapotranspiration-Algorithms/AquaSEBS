# AquaSEBS Computational Methodology

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS` (v1.3.0)

**Primary Code Modules:**
- `AquaSEBS/AquaSEBS.py` (`AquaSEBS()`)
- `AquaSEBS/water_heat_flux.py` (`water_heat_flux()`)
- `AquaSEBS/constants.py`

**Primary Literature References:**
1. **Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).** *Surface Energy Balance of Fresh and Saline Waters: AquaSEBS*. Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
2. **Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W., Pearson, C., & Huntington, J. L. (2023).** *Remotely sensed terrestrial open water evaporation*. Scientific Reports, 13, 8174. https://doi.org/10.1038/s41598-023-34921-2

---

## 1. Architectural Overview

The `AquaSEBS` Python package implements a **hybrid remote-sensing surface energy balance workflow** for open water bodies that couples:
1. **Water Masking:** Automated inland and coastal water identification via NASADEM Surface Water Body (`NASADEM.swb`) extents.
2. **Radiative Forcing:** Broadband net radiation ($R_n$) via the Verma et al. (2016) radiation balance (`verma_net_radiation`).
3. **Water Heat Flux ($W$):** Instantaneous water heat storage/release via the **Equilibrium Temperature Model (ETM)** from Abdelrady et al. (2016, Eqs. 8–13) and Fisher et al. (2023, Eqs. 1–6).
4. **Latent Heat Flux ($LE$):** Energy-limited potential open-water evaporation via the **Priestley–Taylor** formulation ($\alpha = 1.26$) from Fisher et al. (2023, Eq. 7), replacing the iterative Monin–Obukhov Similarity (MOS) sensible heat solver in Abdelrady et al. (2016).
5. **Diurnal Integration (Optional):** Upscaling from instantaneous satellite overpass $LE$ to daylight totals via `daylight_ET_from_instantaneous_LE`.

---

## 2. Step-by-Step Execution Pipeline

### Step 1: Spatial Geometry & Water Pixel Masking

1. **Geometry Inference (`AquaSEBS.py`, `water_heat_flux.py`):**
   - If `geometry` is not explicitly provided and `WST_C` is a `rasters.Raster` instance, `geometry` is inferred from `WST_C.geometry`.
   - Water surface temperature (`WST_C`, in $^\circ\text{C}$) is validated for non-`None` presence and checked via `check_distribution()`.
2. **Water Mask Retrieval & Application:**
   - When `water is None` and `geometry is not None`, the binary surface water body extent is retrieved from `NASADEM.swb(geometry=geometry)`.
   - Inside `water_heat_flux()`, if `mask_non_water_pixels=True` (default in `constants.py`), land pixels are masked out:
     $$\text{WST}_C = \text{where}(\text{water},\, \text{WST}_C,\, \text{NaN})$$
   - Because `WST_C` propagates `NaN` values through all subsequent algebraic operations, all downstream ETM intermediates, water heat flux, latent heat flux, and daylight ET outputs inherit the open-water mask.

---

### Step 2: Meteorological & Radiative Forcing Acquisition

All meteorological and radiative inputs accept user-supplied rasters, NumPy arrays, or scalars. When omitted and both `geometry` and `time_UTC` are provided, missing inputs are retrieved dynamically from NASA's **GEOS-5 FP** assimilation system (`GEOS5FP`, default resampling `"cubic"`):

1. **Near-Surface Air Temperature ($T_a$):**
   - Retrieved via `GEOS5FP.Ta_C(time_UTC, geometry, resampling)` if `Ta_C is None` (raises `MissingOfflineParameter` if `offline_mode=True`).
2. **Broadband Net Radiation ($R_n$):**
   - If `Rn_Wm2 is None` and `albedo`, `emissivity`, `geometry`, and `time_UTC` are available:
     - Incoming shortwave radiation ($\text{SW}_{\text{in}}$, $\text{W}\cdot\text{m}^{-2}$) is fetched via `GEOS5FP.SWin(...)` if `SWin_Wm2 is None`.
     - `verma_net_radiation()` computes net radiation and its shortwave/longwave components (`Rn_Wm2`, `SWnet`, `LWin`, `LWout`) from `SWin_Wm2`, `albedo`, `ST_C=WST_C`, `emissivity`, `Ta_C`, and `RH`.
3. **Dew-Point Temperature ($T_d$):**
   - Inside `water_heat_flux()`, if `Td_C is None`, `Ta_C` and fractional relative humidity (`RH`, $0\text{–}1$) are fetched from `GEOS5FP` if needed, and dew point ($^\circ\text{C}$) is estimated via the linear approximation:
     $$T_d = T_a - \frac{100 - 100 \cdot \text{RH}}{5.0}$$
4. **Wind Speed ($u$):**
   - Retrieved via `GEOS5FP.wind_speed(time_UTC, geometry, resampling)` ($\text{m}\cdot\text{s}^{-1}$) if `windspeed_mps is None`.
5. **Net Shortwave Radiation ($\text{SW}_{\text{net}}$):**
   - Inside `water_heat_flux()`, if `SWnet is None`, `SWin` and broadband `ALBEDO` are retrieved from `GEOS5FP` and combined as:
     $$\text{SW}_{\text{net}} = \text{SW}_{\text{in}} \cdot (1 - \alpha_{\text{surf}})$$

---

### Step 3: Water Heat Flux via Equilibrium Temperature Model (`water_heat_flux.py`)

Instantaneous water heat flux ($W$, $\text{W}\cdot\text{m}^{-2}$) is computed sequentially from Abdelrady et al. (2016, Eqs. 8–13):

1. **Half Temperature Gradient ($T_n$, $^\circ\text{C}$):**
   $$T_n = 0.5 \cdot (\text{WST}_C - T_d)$$
2. **Evaporation Efficiency ($\eta$, dimensionless):**
   $$\eta = 0.35 + 0.015 \cdot \text{WST}_C + 0.0012 \cdot T_n^2$$
3. **Scaled Wind Speed Function ($S$, dimensionless):**
   $$S = 3.3 \cdot u$$
4. **Thermal Exchange Coefficient ($\beta$, $\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$):**
   $$\beta = 4.5 + 0.05 \cdot \text{WST}_C + (\eta + 0.47) \cdot S$$
5. **Equilibrium Temperature ($T_e$, $^\circ\text{C}$):**
   $$T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}$$
6. **Water Heat Flux ($W_{\text{Wm2}}$, $\text{W}\cdot\text{m}^{-2}$):**
   $$W_{\text{Wm2}} = \beta \cdot (T_e - \text{WST}_C)$$

All intermediate and input variables (`WST_C`, `Ta_C`, `RH`, `Td_C`, `windspeed_mps`, `SWnet`, `Tn`, `η`, `S`, `β`, `Te`, `W_Wm2`) are returned in a dictionary and merged into the main `AquaSEBS()` output dictionary.

---

### Step 4: Instantaneous Latent Heat Flux (`AquaSEBS.py`)

1. **Psychrometric Weighting Factor ($\epsilon$, dimensionless):**
   - Computed from air temperature (`Ta_C`) and psychrometric constant (`γ_Pa`, default `GAMMA_PA`) via `priestley_taylor.epsilon_from_Ta_C()`:
     $$\epsilon = \frac{\Delta}{\Delta + \gamma}$$
2. **Priestley–Taylor Latent Heat Flux ($\text{LE}_{\text{Wm2}}$, $\text{W}\cdot\text{m}^{-2}$):**
   - Available energy at the water surface is defined as $R_n - W_{\text{Wm2}}$.
   - `priestley_taylor()` computes potential open-water latent heat flux using $\alpha = 1.26$ (`PT_ALPHA` default):
     $$\text{LE}_{\text{Wm2}} = \alpha \cdot \epsilon \cdot (R_n - W_{\text{Wm2}})$$
   - Sensible heat flux ($H$) is not explicitly computed in code (remaining as the energy balance residual $H = R_n - W - \text{LE}$), and salinity scaling is not yet applied in v1.3.0.

---

### Step 5: Optional Diurnal Upscaling (`AquaSEBS.py`)

When `upscale_to_daylight=True` (default `False` in `constants.py`) and `time_UTC` is provided:
- `daylight_ET_from_instantaneous_LE()` is invoked with:
  - `LE_instantaneous_Wm2 = LE_Wm2`
  - `Rn_instantaneous_Wm2 = Rn_Wm2`
  - `G_instantaneous_Wm2 = W_Wm2`
  - `time_UTC` and `geometry`
- The resulting daylight radiation, latent heat flux, and daily evapotranspiration depth ($\text{mm}\cdot\text{day}^{-1}$) products are merged into the returned `results` dictionary.

---

## 3. Literature vs. Codebase Implementation Matrix

| Model Component | Abdelrady et al. (2016) | Fisher et al. (2023) | Implemented Code (`v1.3.0`) |
| :--- | :--- | :--- | :--- |
| **Water Masking** | Manual wind-direction fetch / static mask | MODIS/Landsat land-water mask | Automated `NASADEM.swb` (`rt.where`) |
| **Net Radiation ($R_n$)** | 4-component radiation input | FLiES/BESS + surface energy balance | `verma_net_radiation` (`GEOS5FP` + thermal) |
| **Water Heat Flux ($W$)** | ETM (Eqs. 8–13) | ETM (Eqs. 1–6) | **Exact match** (`water_heat_flux.py`) |
| **Dew Point ($T_d$)** | Measured / ECMWF reanalysis | MODIS `MOD07_L2` / NCEP reanalysis | Linear $T_a, \text{RH}$ rule-of-thumb |
| **Sensible Heat ($H$)** | Iterative MOS ($z_{0m}=0.0002$, $z_{0h}=0.0001$) | Not solved explicitly | Not solved explicitly |
| **Latent Heat ($\text{LE}$)** | SEBS evaporative fraction ($EF_r$) | Priestley–Taylor ($\alpha = 1.26$) | **Exact match** (`priestley_taylor`) |
| **Salinity Correction ($\sigma$)** | Turk (1970) exponential factor | Noted as future sensitivity | Not implemented (tracked in Issue #11) |
| **Daily Upscaling** | Daily $EF$ scaling | Sinusoidal diurnal $R_n$ + $EF$ | Optional `daylight_ET_from_instantaneous_LE` |