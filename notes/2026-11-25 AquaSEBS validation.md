# AquaSEBS Validation Benchmark Summary

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS` (v1.3.0)

**Primary Validation Studies:**
1. **Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W., Pearson, C., & Huntington, J. L. (2023).** *Remotely sensed terrestrial open water evaporation*. Scientific Reports, 13, 8174. https://doi.org/10.1038/s41598-023-34921-2
2. **Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).** *Surface Energy Balance of Fresh and Saline Waters: AquaSEBS*. Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583

---

## 1. Overview

The operational hybrid formulation of **AquaSEBS** (coupling the Equilibrium Temperature Model for water heat flux with the Priestley–Taylor evaporation formulation) was evaluated in **Fisher et al. (2023)** across **19 global open-water sites** spanning **1986–2019** and **seven Köppen–Geiger climate zones**. This benchmark established operational confidence for un-masking open water evaporation in **NASA's ECOSTRESS Collection 2+** (`L3_ET_PT-JPL`) and the **OpenET** platform.

---

## 2. Validation Datasets (Fisher et al., 2023)

### 2.1 In Situ Reference Network (19 Sites)

Ground-truth observations were compiled from four primary networks and syntheses:
* **Great Lakes Evaporation Network (GLEN):** Offshore eddy covariance towers across the Laurentian Great Lakes (`superiorwatersheds.org/GLEN`).
* **US Bureau of Reclamation Open Water Evaporation Network (OWEN):** Western US reservoir stations (`owen.dri.edu`).
* **Zhao & Gao (2019) Synthesis & Literature Compilations:** Global lake and reservoir evaporation records.
* **Primary Field Measurements:** Collison (2019) floating evaporation pan observations at Cochiti Lake, New Mexico.

**Measurement Breakdown:**
* **Methods:** Eddy covariance ($n = 11$), Bowen ratio energy balance / BREB ($n = 5$), bulk mass transfer ($n = 2$), and floating evaporation pan ($n = 1$).
* **Reporting Cadence:** Half-hourly ($n = 9$), hourly ($n = 1$), daily average ($n = 1$), and monthly average ($n = 8$).
* **Climate Zones Represented:** Humid subtropical (`Cfa`, $n=6$), warm-summer humid continental (`Dfb`, $n=5$), hot desert (`BWh`, $n=3$), cold semi-arid (`BSk`, $n=2$), hot semi-arid (`BSh`, $n=1$), cold desert (`BWk`, $n=1$), and hot-summer humid continental (`Dfa`, $n=1$).

| Site | Location | Coordinates (Lat, Lon) | Climate | Method | Sample Rate |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **American Falls** | ID, USA | $42.900^\circ\text{ N}, 112.758^\circ\text{ W}$ | `BSk` | Bulk mass transfer | 30 min |
| **Calm** | FL, USA | $28.142^\circ\text{ N}, 82.582^\circ\text{ W}$ | `Cfa` | BREB | Monthly avg |
| **Cochiti** | NM, USA | $35.617^\circ\text{ N}, 106.315^\circ\text{ W}$ | `BSk` | Floating evaporation pan | 1 h |
| **Five-O** | FL, USA | $30.422^\circ\text{ N}, 85.664^\circ\text{ W}$ | `Cfa` | BREB | Monthly avg |
| **Granite Island** | MI, USA | $46.721^\circ\text{ N}, 87.414^\circ\text{ W}$ | `Dfb` | Eddy covariance | 30 min |
| **Kasumigaura** | Japan | $36.043^\circ\text{ N}, 140.412^\circ\text{ E}$ | `Cfa` | Eddy covariance | Monthly avg |
| **Kinneret** | Israel | $32.817^\circ\text{ N}, 35.600^\circ\text{ E}$ | `BSh` | BREB | Monthly avg |
| **Lahontan** | NV, USA | $39.450^\circ\text{ N}, 119.069^\circ\text{ W}$ | `BWk` | Bulk mass transfer | 30 min |
| **Long Point** | ON, Canada | $42.567^\circ\text{ N}, 80.050^\circ\text{ W}$ | `Dfb` | Eddy covariance | 30 min |
| **Mead** | NV/AZ, USA | $36.083^\circ\text{ N}, 114.781^\circ\text{ W}$ | `BWh` | Eddy covariance | 30 min |
| **Mohave** | NV/AZ, USA | $35.427^\circ\text{ N}, 114.648^\circ\text{ W}$ | `BWh` | Eddy covariance | 30 min |
| **Nasser** | Egypt | $23.954^\circ\text{ N}, 32.878^\circ\text{ E}$ | `BWh` | BREB | Monthly avg |
| **Ross Barnett** | MS, USA | $32.438^\circ\text{ N}, 90.032^\circ\text{ W}$ | `Cfa` | Eddy covariance | Monthly avg |
| **Spectacle Reef** | MI, USA | $45.776^\circ\text{ N}, 84.136^\circ\text{ W}$ | `Dfb` | Eddy covariance | 30 min |
| **Stannard Rock** | MI, USA | $47.184^\circ\text{ N}, 87.225^\circ\text{ W}$ | `Dfb` | Eddy covariance | 30 min |
| **Starr** | FL, USA | $27.957^\circ\text{ N}, 81.588^\circ\text{ W}$ | `Cfa` | BREB | Daily avg |
| **Taihu** | China | $31.382^\circ\text{ N}, 120.196^\circ\text{ E}$ | `Cfa` | Eddy covariance | Monthly avg |
| **White Bear** | MN, USA | $45.077^\circ\text{ N}, 92.983^\circ\text{ W}$ | `Dfa` | Eddy covariance | Monthly avg |
| **White Shoal** | MI, USA | $45.773^\circ\text{ N}, 85.137^\circ\text{ W}$ | `Dfb` | Eddy covariance | 30 min |

### 2.2 Satellite & Reanalysis Inputs

Because ECOSTRESS Collection 2 back-processing was underway during the study, Fisher et al. (2023) evaluated AquaSEBS using two historical satellite records (~10:30 AM local overpass) paired with reanalysis meteorology:

1. **Landsat 5, 7, and 8 (1986–2019, $30\text{ m}$ spatial resolution):**
   * **Surface Temperature ($T_0$):** Landsat ARD ST resampled to $30\text{ m}$ (from $120\text{ m}$ L5, $60\text{ m}$ L7, and $100\text{ m}$ L8 thermal bands), sampled over a $5 \times 5$ pixel ($150\text{ m} \times 150\text{ m}$) offshore window.
   * **Albedo ($\alpha_{\text{surf}}$):** Derived from Landsat ARD Surface Reflectance using Liang (2001) narrowband-to-broadband coefficients.
   * **Meteorology & Radiation:** NCEP-NCAR Reanalysis I (6-hourly, $2.5^\circ$ resolution) for $R_s$, $T_d$, humidity, and wind speed $u$.

2. **Terra/Aqua MODIS (2000–2019, $1\text{ km}$ spatial resolution):**
   * **Scene Filtering:** 11,016 aligned scenes filtered to **686 cloud-free instantaneous scenes** and **52 cloud-free daily/monthly scenes**.
   * **Surface Temperature ($T_0$):** Daily $1\text{ km}$ `MOD11A1` sampled over a $3 \times 3$ pixel offshore window.
   * **Albedo ($\alpha_{\text{surf}}$):** 16-day $500\text{ m}$ `MCD43A3` BRDF albedo.
   * **Atmosphere & Radiation:** `MOD07_L2` for $T_a$ and $T_d$; `MOD04_L2` (aerosol optical thickness) and `MOD06_L2` (cloud optical thickness) driving FLiES/BESS radiative transfer for $R_s$ and $R_n$; NCEP-NCAR Reanalysis I for wind speed $u$.

---

## 3. Quantitative Performance Metrics

### 3.1 Fisher et al. (2023) Global Satellite Validation

| Evaluation Regime | Sensor | $r^2$ | RMSE | Bias | Linear Fit |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Instantaneous (High-wind filtered, $u \le 7.5\text{ m}\cdot\text{s}^{-1}$)** | MODIS ($1\text{ km}$) | $0.71$ | $53.7\text{ W}\cdot\text{m}^{-2}$ ($38\%$ of mean) | $-19.14\text{ W}\cdot\text{m}^{-2}$ ($13\%$ of mean) | $y = 0.8x + 3.0$ |
| **Instantaneous (Unfiltered, all wind speeds)** | MODIS ($1\text{ km}$) | $0.47$ | $84.4\text{–}95.4\text{ W}\cdot\text{m}^{-2}$ ($62\%$ of mean) | $-49.5\text{ W}\cdot\text{m}^{-2}$ ($36\%$ of mean) | $y = 0.6x + 13.0$ |
| **Daily Evaporation** | MODIS ($1\text{ km}$) | $0.47$ | $1.5\text{ mm}\cdot\text{day}^{-1}$ ($41\%$ of mean) | $-0.03\text{ to }+0.19\text{ mm}\cdot\text{day}^{-1}$ ($1\%$ of mean) | $y = 0.96x + 0.4$ |
| **Daily Evaporation** | Landsat ($30\text{ m}$) | $0.56$ | $1.2\text{ mm}\cdot\text{day}^{-1}$ ($38\%$ of mean) | $-0.82\text{ mm}\cdot\text{day}^{-1}$ ($26\%$ of mean) | $y = 0.78x + 0.1$ |

*(Note: In Fisher et al., 2023, prose and figure captions report $84.4$ vs. $95.4\text{ W}\cdot\text{m}^{-2}$ for unfiltered instantaneous RMSE and $+0.19$ vs. $-0.03\text{ mm}\cdot\text{day}^{-1}$ for MODIS daily bias.)*

### 3.2 Abdelrady et al. (2016) Reanalysis & Field Benchmarks

* **Water Heat Flux ($G_{0w}$) vs. ECMWF ERA-Interim:** Across five seas/oceans (Indian, Atlantic, Mediterranean, North Sea, Red Sea; $n \approx 1005$), ETM achieved $\text{rRMSE} = 4.3\%\text{–}8.5\%$ ($\text{RMSE} = 26.5\text{–}68.5\text{ W}\cdot\text{m}^{-2}$).
* **Freshwater Lakes:**
  * **Lake Tana (in situ, $n = 23$):** $\text{LE}$ $r^2 = 0.98$, $\text{RMSE} = 35.6\text{ W}\cdot\text{m}^{-2}$ ($\text{rRMSE} = 10.3\%$); sensible heat $H$ $\text{RMSE} = 4.8\text{ W}\cdot\text{m}^{-2}$ ($\text{rRMSE} = 7.9\%$).
  * **Lake Victoria (vs. ECMWF, $n = 1005$):** Daily evaporation $\text{RMSE} = 1.5\text{ mm}\cdot\text{day}^{-1}$ ($\text{rRMSE} = 7.4\%$), water heat flux $\text{RMSE} = 19.9\text{ W}\cdot\text{m}^{-2}$ ($\text{rRMSE} = 6.7\%$).
* **Saline & Hypersaline Waters:** Incorporating the Turk (1970) salinity factor reduced Great Salt Lake ($S = 240\text{ g}\cdot\text{L}^{-1}$) evaporation by $26.5\%$, improving 3-hourly $\text{RMSE}$ from $0.63$ to $0.38\text{ mm}\cdot(3\text{ h})^{-1}$.

---

## 4. Key Findings & Error Attribution

### 4.1 High-Wind Regime Shift ($u > 7.5\text{ m}\cdot\text{s}^{-1}$)
* Open-water evaporation shifts from **radiatively controlled** under calm/moderate winds to **atmospherically (aerodynamically) controlled** during high-wind events (mean daily $u > 7.5\text{ m}\cdot\text{s}^{-1}$).
* Because the operational Priestley–Taylor formulation is energy-limited ($\alpha = 1.26$), unscreened high-wind events cause systematic instantaneous underestimation ($r^2$ drops from $0.71$ to $0.47$; bias worsens from $-19.1$ to $-49.5\text{ W}\cdot\text{m}^{-2}$).
* Temporal integration to daily sums dampens short-term wind spikes, keeping daily $\text{RMSE}$ between $1.2$ and $1.5\text{ mm}\cdot\text{day}^{-1}$ with near-zero bias on MODIS.

### 4.2 Spatial Resolution Gains (30 m Landsat vs. 1 km MODIS)
* Higher spatial resolution ($30\text{ m}$ Landsat) improved daily $r^2$ ($0.56$ vs. $0.47$) and lowered $\text{RMSE}$ ($1.2$ vs. $1.5\text{ mm}\cdot\text{day}^{-1}$), while resolving within-lake spatial gradients driven by bathymetry, shoreline proximity, and circulation currents (e.g., Lahontan Reservoir, American Falls Reservoir, and White Shoal).
* Three outlier scenes at Lake Five-O lowered the Landsat $r^2$ from $0.71$ to $0.56$.

### 4.3 Machine Learning Benchmark & Statistical Predictors of Error
* **Process Model Robustness:** None of 11 benchmarked machine learning models (OLS, Ridge, LASSO, Elastic Net, MLP, TensorFlow Neural Network, Decision Tree, Random Forest, SVM, KNN, Gradient Boosting) outperformed mechanistic AquaSEBS across the 19 sites, confirming that AquaSEBS captures the primary physical variance without site-specific overfitting.
* **Residual Predictability:** A two-layer TensorFlow Neural Network trained to predict instantaneous model residuals ($\text{In Situ} - \text{Satellite}$) explained $74\%$ of the error variance ($r^2 = 0.74$).
* **OLS Error Predictors:** Linear regression of instantaneous residuals against ancillary variables showed:
  * **Wind speed ($\text{m}\cdot\text{s}^{-1}$):** Dominant error driver ($\text{coef} = +20.10$, $t = 11.13$, $p < 0.001$).
  * **Relative humidity ($\%$):** Secondary error driver ($\text{coef} = -1.56$, $t = -3.92$, $p < 0.001$).
  * **Wind direction ($^\circ$):** Minor effect ($\text{coef} = +0.14$, $t = 2.40$, $p = 0.02$).
  * **Waterbody area, shape, air temperature, and net radiation:** No statistically significant bias ($p > 0.19$).

### 4.4 Attributed Sources of Model–Data Mismatch
Fisher et al. (2023) attribute residual discrepancies across five interacting factors:
1. **Model formulation limits** during high-wind aerodynamic regimes and unmodeled salinity/bathymetric storage effects.
2. **In situ measurement uncertainty** and energy-balance closure errors across eddy covariance, BREB, bulk transfer, and floating pan instrumentation.
3. **Coarse meteorological forcing errors** (e.g., $2.5^\circ$ NCEP-NCAR reanalysis wind speed).
4. **Spatiotemporal scale mismatch** between dynamic upwind tower footprints and fixed satellite pixels or instantaneous overpass times (~10:30 AM).
5. **Lateral advection** of warm, dry air across sharp wet/dry shoreline boundaries in arid reservoirs.