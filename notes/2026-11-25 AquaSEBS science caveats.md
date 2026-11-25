# AquaSEBS Scientific Caveats and Methodological Limitations

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS` (v1.3.0)  
**Primary Literature References:**
1. Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016). *Surface Energy Balance of Fresh and Saline Waters: AquaSEBS*. Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
2. Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W., Pearson, C., & Huntington, J. L. (2023). *Remotely sensed terrestrial open water evaporation*. Scientific Reports, 13, 8174. https://doi.org/10.1038/s41598-023-34921-2

---

## 1. Overview

The `AquaSEBS` Python package implements a hybrid satellite remote sensing framework for estimating open-water surface energy balance components. It couples the **Equilibrium Temperature Model (ETM)** parameterization for instantaneous water heat flux ($W$, denoted $G_{0w}$ in Abdelrady et al., 2016) with the **Priestley–Taylor** energy-limited latent heat flux formulation ($LE$) evaluated across 19 global validation sites in Fisher et al. (2023).

While validated for operational use in missions such as NASA ECOSTRESS and OpenET, users and downstream modelers should account for the following physical, mathematical, and empirical caveats when interpreting outputs.

---

## 2. Core Methodological & Physical Caveats

### 2.1 Absence of Water Salinity Correction in Code Implementation
* **Physical Mechanism:** Dissolved salts lower the saturation vapor pressure over water surfaces in accordance with Raoult's law, reducing evaporation rates and reallocating available energy toward water warming and sensible heat flux.
* **Literature Formulation:** Abdelrady et al. (2016, Eq. 19) incorporate the empirical Turk (1970) salinity reduction factor:
  $$\sigma = 1.025 - 0.0246 \exp(0.00879 \cdot S)$$
  where $S$ is salinity ($\text{g}\cdot\text{L}^{-1}$). While negligible in fresh water ($<1\%$ reduction for $S < 20\text{ g}\cdot\text{L}^{-1}$), this reduction reaches $\sim 1\%$ over open ocean ($34.7\text{ g}\cdot\text{L}^{-1}$) and $\sim 27\%$ over hypersaline water bodies such as the Great Salt Lake ($240\text{ g}\cdot\text{L}^{-1}$).
* **Current Code Behavior:** Neither `AquaSEBS()` nor `water_heat_flux()` currently accepts a salinity input or applies $\sigma$. All calculations implicitly assume freshwater conditions ($S = 0$), which will systematically overestimate evaporation over brackish, marine, and hypersaline water bodies.

### 2.2 High-Wind Regime Shift ($u > 7.5\text{ m}\cdot\text{s}^{-1}$)
* **Physical Mechanism:** Open-water evaporation transitions between two distinct physical regimes:
  1. **Radiatively controlled** under calm-to-moderate winds, where available net energy ($R_n - W$) governs latent heat flux.
  2. **Atmospherically (aerodynamically) controlled** during strong wind events, where turbulent mechanical mixing and vapor pressure deficit drive rapid evaporation even when net radiation is low.
* **Model Limitation:** Unlike the original Abdelrady et al. (2016) formulation—which iteratively solves Monin–Obukhov Similarity (MOS) equations with calibrated water roughness lengths ($z_{0m} = 0.0002\text{ m}$, $z_{0h} = 0.0001\text{ m}$)—this implementation uses the energy-limited Priestley–Taylor equation ($\alpha = 1.26$) for $LE$ and only incorporates wind speed inside the ETM water heat flux parameterization ($S(W) = 3.3u$).
* **Empirical Impact:** Fisher et al. (2023) demonstrated that unscreened high-wind events ($u > 7.5\text{ m}\cdot\text{s}^{-1}$ daily mean) severely degrade instantaneous retrieval accuracy:
  * **High-wind filtered ($u \le 7.5\text{ m}\cdot\text{s}^{-1}$):** $r^2 = 0.71$, $\text{RMSE} = 53.7\text{ W}\cdot\text{m}^{-2}$ ($38\%$ of mean), $\text{Bias} = -19.1\text{ W}\cdot\text{m}^{-2}$ ($13\%$).
  * **Unfiltered (all wind speeds):** $r^2 = 0.47$, $\text{RMSE} = 84.4\text{ W}\cdot\text{m}^{-2}$ ($62\%$ of mean), $\text{Bias} = -49.5\text{ W}\cdot\text{m}^{-2}$ ($36\%$).
  * Wind speed is the single largest statistical predictor of instantaneous AquaSEBS error ($\text{coef} = 20.10$, $p < 0.001$). Temporal integration to daily totals dampens this sensitivity ($\text{RMSE} \approx 1.2\text{--}1.5\text{ mm}\cdot\text{day}^{-1}$).

### 2.3 Skin Temperature Approximation vs. Bathymetry and Thermal Storage Lag
* **Physical Mechanism:** Water bodies absorb shortwave radiation across depth rather than at an opaque surface, and vertical convective mixing redistributes heat across a mixed layer bounded by lake bathymetry.
* **Model Limitation:** Satellite thermal infrared sensors observe only the top skin temperature ($T_0$, top $\sim 10\text{--}100\,\mu\text{m}$). The Equilibrium Temperature Model estimates instantaneous $W$ purely from $T_0$ and surface meteorology without accounting for water depth, stratification, thermal inertia, or horizontal heat advection by currents.
* **Empirical Impact:**
  * Deep oceans maintain nearly constant diurnal latent heat flux by adjusting mixed-layer depth, whereas shallow inland lakes exhibit pronounced diurnal hysteresis that is not explicitly parameterized.
  * Abdelrady et al. (2016) note that errors in $W$ are the dominant source of propagated uncertainty in overall energy balance closure, and that AquaSEBS tends to slightly underestimate nocturnal fluxes when $R_n < 0$.

### 2.4 Linear Dew-Point Temperature Approximation
* **Model Formulation:** When dew-point temperature (`Td_C`) is not directly supplied to `water_heat_flux()`, it is estimated from air temperature (`Ta_C`) and fractional relative humidity (`RH`) using a simplified linear rule-of-thumb:
  $$T_d \approx T_a - \frac{100 - 100 \cdot \text{RH}}{5.0}$$
* **Caveat:** While reasonably accurate in humid boundary layers ($\text{RH} > 50\%$), this linear approximation diverges from exact Magnus–Tetens or psychrometric saturation vapor pressure inversions in arid/semi-arid environments ($\text{RH} < 50\%$). Because $T_d$ enters $T_n$, $\eta$, $\beta$, and $T_e$, dew-point errors propagate directly into $W$.

### 2.5 Sensitivity to Water Surface Temperature ($T_0$) Accuracy
* **Sensitivity:** Both Abdelrady et al. (2016) and Fisher et al. (2023) identify water surface temperature ($T_0$) as the most sensitive input variable in the ETM formulation.
* **Operational Requirement:** Reliable retrievals require thermal infrared $T_0$ accuracy within $\pm 0.3\text{--}0.5\text{ K}$ and rigorous cloud/cloud-shadow masking, as undetected thin cirrus or cold cloud edges introduce severe artificial gradients between $T_0$ and $T_e$.

---

## 3. Spatial, Temporal, and Validation Caveats

### 3.1 Shoreline Adjacency and Horizontal Advection
* **Mixed Pixels:** Coarse thermal pixels (e.g., 1 km MODIS) near shorelines or over narrow reservoirs suffer from land-temperature contamination. Even with static water masks (e.g., NASADEM SWB), seasonal reservoir drawdown can expose shoreline sediments that invalidate static masks.
* **Lateral Advection:** Across sharp wet/dry boundaries (particularly in arid reservoirs), warm, dry air advected from surrounding land enhances actual evaporation near the upwind shore in ways not captured by 1D vertical energy balance assumptions.

### 3.2 Diurnal Upscaling Assumptions
* **Constant Evaporative Fraction:** When `upscale_to_daylight=True`, instantaneous $LE$ at satellite overpass time is scaled to daily totals assuming self-preservation of the evaporative fraction ($EF$) and a sinusoidal diurnal radiation curve.
* **Caveat:** Over water bodies with strong diurnal wind cycles (e.g., afternoon lake breezes or frontal passages) or shallow thermal lag, instantaneous overpass $EF$ (typically $\sim\text{10:30 AM}$ for Landsat/MODIS) may not represent full-day integrated energy partitioning.

### 3.3 Geographic & Climatic Representativeness of Validation
* Fisher et al. (2023) validated AquaSEBS across 19 global sites spanning 7 Köppen–Geiger climate zones, representing one of the largest open-water validation benchmarks to date.
* However, the validation dataset lacks representation in **high-latitude/polar regions**, **tropical equatorial lakes**, and **high-altitude water bodies**, and relies on heterogeneous in situ reference methods (eddy covariance, Bowen ratio, bulk mass transfer, floating pans) that carry their own measurement and energy-closure uncertainties.

---

## 4. Summary of Operational Guidance

| Condition / Application | Expected Behavior & Recommendation |
| :--- | :--- |
| **Freshwater lakes & reservoirs ($u \le 7.5\text{ m}\cdot\text{s}^{-1}$)** | Optimal operating regime; expect daily RMSE $\approx 1.2\text{--}1.5\text{ mm}\cdot\text{day}^{-1}$ and instantaneous RMSE $\approx 50\text{--}55\text{ W}\cdot\text{m}^{-2}$. |
| **High-wind conditions ($u > 7.5\text{ m}\cdot\text{s}^{-1}$)** | Flag instantaneous retrievals as low-confidence; rely on daily integrated estimates where wind-driven errors average down. |
| **Saline / Marine / Hypersaline waters** | Apply the Turk (1970) salinity reduction factor $\sigma(S)$ post-hoc until native salinity support is added to the API. |
| **Arid environments ($\text{RH} < 50\%$)** | Supply explicit psychrometric `Td_C` rather than relying on the internal linear `Ta_C`/`RH` approximation. |
| **Small or fluctuating reservoirs** | Use high-resolution thermal data (e.g., 30 m Landsat, 70 m ECOSTRESS) with dynamic water masks to avoid shoreline land contamination. |