# AquaSEBS Scientific Background and Theoretical Framework

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS` (v1.3.0)

**Primary Literature References:**
1. **Abdelrady, A., Timmermans, J., Vekerdy, Z., & Salama, M. S. (2016).** *Surface Energy Balance of Fresh and Saline Waters: AquaSEBS*. Remote Sensing, 8(7), 583. https://doi.org/10.3390/rs8070583
2. **Fisher, J. B., Dohlen, M. B., Halverson, G. H., Collison, J. W., Pearson, C., & Huntington, J. L. (2023).** *Remotely sensed terrestrial open water evaporation*. Scientific Reports, 13, 8174. https://doi.org/10.1038/s41598-023-34921-2

---

## 1. Motivation & Operational Context

Evaporation from open water bodies—spanning oceans, natural lakes, reservoirs, and inland seas—plays a governing role in global and regional hydrological cycles. Oceans account for roughly $86\%$ of global evaporation feeding the atmospheric water cycle, while terrestrial reservoirs and lakes are critical freshwater resources increasingly stressed by climate variability and human abstraction.

Measuring open-water evaporation in situ is logistically difficult and spatially unrepresentative:
* **Instrumentation challenges:** Eddy covariance towers require stable platforms and strict fetch/footprint constraints; Bowen ratio systems lose sub-daily temporal fidelity; and bulk aerodynamic/mass-transfer formulations depend heavily on site-specific empirical wind coefficients.
* **Spatial heterogeneity:** Evaporation rates vary markedly across a single water body due to underlying bathymetry, subsurface heat storage and circulation, shoreline geometry, and horizontal advection of warm, dry air from surrounding land.

Satellite thermal infrared remote sensing resolves spatial gradients in water surface temperature ($T_0$) at fine scales ($30\text{ m}$ Landsat, $70\text{ m}$ ECOSTRESS, $1\text{ km}$ MODIS). **AquaSEBS** (Aquatic Surface Energy Balance System) was developed to adapt terrestrial remote-sensing energy balance models to open water surfaces and now serves as the operational open-water evaporation algorithm for **NASA's ECOSTRESS mission** (`L3_ET_PT-JPL` Collection 2+) and the **OpenET** platform.

---

## 2. Evolution of the AquaSEBS Framework

### 2.1 Original Formulation: Abdelrady et al. (2016)
Abdelrady et al. (2016) adapted Su's (2002) one-source **Surface Energy Balance System (SEBS)** for fresh and saline water bodies through five modifications:
1. Replacing constant-fraction ground heat flux approximations ($G_0 = 0.5 R_n$) with Edinger et al.'s (1968) **Equilibrium Temperature Model (ETM)** for instantaneous water heat flux ($G_{0w}$).
2. Calibrating the aerodynamic roughness length for momentum transfer ($z_{0m} = 0.0002\text{ m}$, with zero-plane displacement $d_0 = 0$) over water surfaces.
3. Inverting Monin–Obukhov Similarity (MOS) equations against eddy covariance measurements over Lake IJsselmeer (The Netherlands) to determine the roughness length for heat transfer ($z_{0h} = 0.0001\text{ m}$, corresponding to $kB^{-1} \approx 0.3$).
4. Defining water-specific wet (open freshwater) and dry (salt-saturated water at $\text{RH} > 70\%$) boundary limits for evaporative fraction ($EF$) upscaling.
5. Incorporating Turk's (1970) empirical salinity reduction factor ($\sigma$) to account for Raoult's law vapor pressure depression over saline and hypersaline waters.

### 2.2 Operational Remote-Sensing Hybrid: Fisher et al. (2023)
To enable global, operational satellite retrievals across millions of water bodies without requiring site-calibrated aerodynamic resistance iterations, Fisher et al. (2023) coupled:
* **Abdelrady et al.'s (2016) exact 6-equation Equilibrium Temperature Model (ETM)** to compute instantaneous water heat flux ($W$, denoted $G_{0w}$ in the literature).
* **The Priestley–Taylor (1972) energy-limited formulation** ($\alpha = 1.26$) to partition available energy ($R_n - W$) into instantaneous latent heat flux ($LE$).
* **Sinusoidal diurnal radiation scaling** (Bisht et al., 2005; Fisher, 2018) to integrate instantaneous satellite overpass retrievals into daily evapotranspiration totals.

This hybrid architecture is what is implemented in the `AquaSEBS` Python package.

---

## 3. Mathematical Formulation

### 3.1 Surface Energy Balance
At the water–atmosphere interface, the instantaneous energy balance is governed by:

$$R_n = LE + H + W$$

where:
* $R_n$ is net broadband radiation ($\text{W}\cdot\text{m}^{-2}$),
* $LE$ (or $\lambda E$) is latent heat flux consumed by evaporation ($\text{W}\cdot\text{m}^{-2}$),
* $H$ is sensible heat flux exchanged with the atmosphere ($\text{W}\cdot\text{m}^{-2}$), and
* $W$ (or $G_{0w}$) is water heat flux stored in or released from the water column ($\text{W}\cdot\text{m}^{-2}$).

---

### 3.2 Net Radiation ($R_n$)
When not supplied directly, $R_n$ is computed via the Verma et al. (2016) surface radiation balance:

$$R_n = \text{SW}_{\text{in}}(1 - \alpha_{\text{surf}}) + \text{LW}_{\text{in}} - \text{LW}_{\text{out}}$$

where $\text{SW}_{\text{in}}$ is incoming shortwave radiation ($\text{W}\cdot\text{m}^{-2}$), $\alpha_{\text{surf}}$ is broadband surface albedo (dimensionless), and $\text{LW}_{\text{in}}, \text{LW}_{\text{out}}$ are downwelling and upwelling longwave radiation fluxes ($\text{W}\cdot\text{m}^{-2}$) derived from water surface temperature ($T_0$), surface emissivity ($\varepsilon$), near-surface air temperature ($T_a$), and relative humidity ($\text{RH}$).

---

### 3.3 Water Heat Flux via Equilibrium Temperature Model ($W$)
In situ water heat flux is defined by the depth-integrated heat storage rate:

$$W = \rho_w c_w \int_{0}^{Z} \frac{\Delta T}{\Delta t}\, dZ$$

Because thermal infrared satellites observe only the top skin temperature ($T_0$, upper $\sim 10\text{--}100\,\mu\text{m}$) rather than vertical temperature profiles, AquaSEBS applies the **Equilibrium Temperature Model (ETM)** (Edinger et al., 1968; Abdelrady et al., 2016, Eqs. 8–13; Fisher et al., 2023, Eqs. 1–6).

The **equilibrium temperature** ($T_e$) is defined as the hypothetical water surface temperature at which net heat exchange between the water surface and the atmosphere equals zero. It is solved in six sequential steps:

1. **Half Temperature Gradient ($T_n$):**
   $$T_n = 0.5\,(T_0 - T_d)$$
   where $T_0$ is water surface skin temperature ($^\circ\text{C}$) and $T_d$ is near-surface air dew-point temperature ($^\circ\text{C}$).

2. **Evaporation Efficiency Factor ($\eta$):**
   $$\eta = 0.35 + 0.015\,T_0 + 0.0012\,T_n^2$$
   accounting for baseline neutral exchange ($0.35$), linear thermal dependence ($0.015\,T_0$), and non-linear vapor gradient enhancement ($0.0012\,T_n^2$).

3. **Scaled Wind Speed Function ($S(W)$):**
   $$S(W) = 3.3\,u$$
   where $u$ is near-surface wind speed ($\text{m}\cdot\text{s}^{-1}$) and $3.3$ is an empirical wind-mixing coefficient over open water.

4. **Thermal Exchange Coefficient ($\beta$):**
   $$\beta = 4.5 + 0.05\,T_0 + (\eta + 0.47)\,S(W)$$
   representing total interfacial heat conductance ($\text{W}\cdot\text{m}^{-2}\cdot^\circ\text{C}^{-1}$) combining radiative/conductive baseline transfer ($4.5 + 0.05\,T_0$) with turbulent wind-driven latent and sensible exchange.

5. **Equilibrium Temperature ($T_e$):**
   $$T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}$$
   where $\text{SW}_{\text{net}} = \text{SW}_{\text{in}}(1 - \alpha_{\text{surf}})$ is net absorbed shortwave radiation ($\text{W}\cdot\text{m}^{-2}$).

6. **Instantaneous Water Heat Flux ($W$):**
   $$W = \beta\,(T_e - T_0)$$
   Positive $W$ ($T_e > T_0$) indicates net energy storage (water warming); negative $W$ ($T_e < T_0$) indicates heat release from the water column to the atmosphere (water cooling).

---

### 3.4 Latent Heat Flux via Priestley–Taylor ($LE$)
Over open water surfaces where moisture supply is unlimited and surface resistance is zero, instantaneous latent heat flux is computed using the Priestley–Taylor (1972) formulation (Fisher et al., 2023, Eq. 7):

$$LE = \alpha \left(\frac{\Delta}{\Delta + \gamma}\right)(R_n - W)$$

where:
* $\alpha = 1.26$ is the Priestley–Taylor advection-aridity coefficient (dimensionless),
* $\Delta$ is the slope of the saturation vapor pressure–temperature curve at air temperature $T_a$ ($\text{Pa}\cdot^\circ\text{C}^{-1}$ or $\text{kPa}\cdot^\circ\text{C}^{-1}$),
* $\gamma$ is the psychrometric constant ($\approx 66\text{ Pa}\cdot^\circ\text{C}^{-1}$ or $0.066\text{ kPa}\cdot^\circ\text{C}^{-1}$), and
* $\epsilon = \frac{\Delta}{\Delta + \gamma}$ is the dimensionless psychrometric energy-partitioning factor.

---

### 3.5 Salinity Reduction Factor ($\sigma$)
Dissolved salts reduce the chemical activity and saturation vapor pressure of water according to Raoult's law, suppressing evaporation and reallocating available energy into water heat storage and sensible heat flux. Abdelrady et al. (2016, Eqs. 19–20) parameterize this effect using Turk's (1970) empirical brine evaporation curve:

$$\sigma = 1.025 - 0.0246 \exp(0.00879 \cdot S)$$

$$E_s = \sigma \cdot E_{\text{fresh}}$$

where $S$ is water salinity ($\text{g}\cdot\text{L}^{-1}$) and $\sigma$ is the dimensionless salinity reduction factor:
* **Fresh to brackish water ($0\text{–}20\text{ g}\cdot\text{L}^{-1}$):** $< 1\%$ reduction ($\sigma \approx 0.99\text{–}1.00$).
* **Mean open ocean ($34.7\text{ g}\cdot\text{L}^{-1}$):** $\sim 1\%$ reduction ($\sigma \approx 0.99$).
* **Hypersaline lakes (e.g., Great Salt Lake, $240\text{ g}\cdot\text{L}^{-1}$):** $\sim 22\%\text{–}27\%$ reduction ($\sigma \approx 0.78$).
* **Extreme brine ($300\text{ g}\cdot\text{L}^{-1}$):** $\sim 32\%$ reduction ($\sigma \approx 0.68$).

*(Note: Native salinity scaling is tracked as an enhancement in Issue #11 and is not yet applied inside `AquaSEBS()` v1.3.0.)*

---

### 3.6 Diurnal Upscaling to Daily Evaporation
Over water surfaces, the evaporative fraction ($EF$) exhibits high daytime stability ($\text{STD} \approx 0.05$ in Abdelrady et al., 2016) due to thermal buffering by the water column. When `upscale_to_daylight=True`, AquaSEBS integrates instantaneous overpass $LE$ to daily totals following the ECOSTRESS Collection 1 / Bisht et al. (2005) sinusoidal diurnal radiation framework, preserving the instantaneous ratio between latent heat flux and available energy ($R_n - W$) across daylight hours.

---

## 4. Empirical Validation & Physical Behavior

### 4.1 Summary of Literature Benchmarks

| Study | Scope & Data Sources | Key Findings & Error Metrics |
| :--- | :--- | :--- |
| **Abdelrady et al. (2016)** | Eddy covariance at Lake IJsselmeer (NL); in situ data at Lake Tana (ET); ECMWF ERA-Interim across 5 seas/oceans, Lake Victoria, and Great Salt Lake | • Water heat flux $G_{0w}$ vs. ECMWF: $\text{rRMSE} = 4.3\%\text{–}8.5\%$ ($\text{RMSE} = 26.5\text{–}68.5\text{ W}\cdot\text{m}^{-2}$)<br>• Lake Tana $LE$: $R^2 = 0.98$, $\text{RMSE} = 35.6\text{ W}\cdot\text{m}^{-2}$ ($\text{rRMSE} = 10.3\%$)<br>• Lake Victoria daily $E$: $\text{RMSE} = 1.5\text{ mm}\cdot\text{day}^{-1}$ ($\text{rRMSE} = 7.4\%$)<br>• Great Salt Lake salinity correction reduced evaporation by $26.5\%$ |
| **Fisher et al. (2023)** | **19 global in situ sites** across 7 Köppen–Geiger climate zones (GLEN, OWEN, flux towers); 686 cloud-free MODIS ($1\text{ km}$) scenes & Landsat ($30\text{ m}$) scenes (1986–2019) | • **Instantaneous MODIS ($u \le 7.5\text{ m}\cdot\text{s}^{-1}$):** $r^2 = 0.71$, $\text{RMSE} = 53.7\text{ W}\cdot\text{m}^{-2}$ ($38\%$), $\text{Bias} = -19.1\text{ W}\cdot\text{m}^{-2}$ ($13\%$)<br>• **Instantaneous MODIS (unfiltered winds):** $r^2 = 0.47$, $\text{RMSE} = 84.4\text{ W}\cdot\text{m}^{-2}$ ($62\%$), $\text{Bias} = -49.5\text{ W}\cdot\text{m}^{-2}$ ($36\%$)<br>• **Daily MODIS:** $\text{RMSE} = 1.5\text{ mm}\cdot\text{day}^{-1}$, $\text{Bias} = +0.19\text{ mm}\cdot\text{day}^{-1}$ ($1\%$)<br>• **Daily Landsat ($30\text{ m}$):** $r^2 = 0.56$, $\text{RMSE} = 1.2\text{ mm}\cdot\text{day}^{-1}$<br>• **ML Benchmark:** None of 11 ML models beat AquaSEBS; neural network predicted residuals with $r^2 = 0.74$ |

### 4.2 Governing Physical Controls & Regime Transitions
1. **Radiative vs. Aerodynamic Control:** At daily-to-seasonal timescales and moderate wind speeds ($u \le 7.5\text{ m}\cdot\text{s}^{-1}$), open-water evaporation is primarily **radiatively controlled**, aligning closely with the Priestley–Taylor assumption. During short-term high-wind events ($u > 7.5\text{ m}\cdot\text{s}^{-1}$), evaporation shifts to an **atmospherically/aerodynamically controlled** regime where turbulent vapor transport dominates instantaneous fluxes.
2. **Depth & Thermal Inertia:** Deep oceans absorb and redistribute solar energy through changes in convective mixed-layer thickness, maintaining nearly flat diurnal $LE$ curves. Shallow lakes have depth-limited mixing layers, producing strong diurnal hysteresis in water heat storage and turbulent fluxes.
3. **Input Sensitivity Hierarchy:** Water heat flux ($W$) is the primary source of propagated energy-balance uncertainty and is most sensitive to **water surface temperature ($T_0$)** (requiring thermal accuracy $< 0.3\text{--}0.5\text{ K}$), followed by **dew-point temperature ($T_d$)** and **wind speed ($u$)**.

---

## 5. Additional Supporting References

* **Bisht, G., Venturini, V., Islam, S., & Jiang, L. (2005).** Estimation of the net radiation using MODIS. *Remote Sensing of Environment*, 97(1), 52–67.
* **Edinger, J. E., Duttweiler, D. W., & Geyer, J. C. (1968).** The response of water temperatures to meteorological conditions. *Water Resources Research*, 4(5), 1137–1143.
* **Fisher, J. B. (2018).** *Level-3 Evapotranspiration L3 (ET_PT-JPL) Algorithm Theoretical Basis Document (ECOSTRESS)*. Jet Propulsion Laboratory, California Institute of Technology.
* **Melton, F. S., et al. (2022).** OpenET: Filling the biggest data gap in water management for the Western United States. *JAWRA*, 58(6), 971–994.
* **Priestley, C. H. B., & Taylor, R. J. (1972).** On the assessment of surface heat flux and evaporation using large-scale parameters. *Monthly Weather Review*, 100(2), 81–92.
* **Su, Z. (2002).** The Surface Energy Balance System (SEBS) for estimation of turbulent heat fluxes. *Hydrology and Earth System Sciences*, 6(1), 85–99.
* **Turk, L. J. (1970).** Evaporation of brine: A field study on the Bonneville Salt Flats, Utah. *Water Resources Research*, 6(4), 1209–1215.
* **Verma, M., et al. (2016).** Global surface net-radiation at 5 km from MODIS Terra. *Remote Sensing*, 8(9), 739.