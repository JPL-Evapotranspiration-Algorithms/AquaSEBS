# AquaSEBS Comprehensive Development Plan & Release Roadmap

**Package:** `AquaSEBS` (Aquatic Surface Energy Balance System)  
**Repository:** [JPL-Evapotranspiration-Algorithms/AquaSEBS](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS)  
**Author & Maintainer:** Gregory H. Halverson ([gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)), NASA Jet Propulsion Laboratory  
**Current Target Release:** `v1.4.0`  
**Operational Baseline Compatibility:** NASA ECOSTRESS Collection 3 (`BuildID 0700`)  
**Default Development Runtime:** Python `3.14` (CI support matrix: Python `3.11` – `3.14`)  

---

## 1. Executive Summary

`AquaSEBS` implements the open-water surface energy balance and evaporation model combining the Equilibrium Temperature Model (ETM) of Edinger et al. (1968) and Abdelrady et al. (2016) with the energy-limited Priestley–Taylor (1972) evaporation parameterization validated globally across 19 open-water sites by Fisher et al. (2023) for NASA's **ECOSTRESS** (`L3_ET_PT-JPL` Collection 2+) and **OpenET** systems.

This document establishes the complete engineering and scientific development plan across:
1. **Completed Historical Releases (`v1.0.1` – `v1.3.0`)**: Foundational physics implementation, diurnal upscaling, initial `offline_mode` hooks, and input configuration constants.
2. **Completed `v1.4.0` Engineering Foundations (`#6`–`#10`)**: Full `offline_mode` enforcement, lazy `GEOS5FP` initialization, parameter forwarding (`RH`, `SWin_Wm2`, `albedo`), top-level namespace cleanup, comprehensive `pytest` unit test coverage, Python `3.14` standardization, and automated version management (`make set-version`).
3. **Pending `v1.4.0` Operational & Scientific Deliverables (`#11`–`#16`)**:
   * **Operational Production Fixes (`#14`, `#15`)**: Graceful all-`NaN` handling when scenes contain zero valid water pixels (eliminating `BlankOutputError` aborts in `JET3` Issue `#92`) and updating `GEOS5FP.wind_speed()` to `GEOS5FP.wind_speed_mps()`.
   * **Opt-In / Non-Breaking Scientific Enhancements (`#11`, `#12`, `#13`)**: Native Turk (1970) salinity correction (`salinity_gL=None` by default), opt-in Magnus–Tetens dew-point calculation (preserving the linear approximation by default so ECOSTRESS Collection 3 `BuildID 0700` baseline outputs are unchanged), and non-intrusive high-wind diagnostic flagging (`high_wind_flag`).
   * **Documentation & Release Finalization (`#16`)**: Escaped LaTeX cleanup in `README.md`, full test/doctest verification, closure of `#11`–`#16`, and PyPI/GitHub release of `v1.4.0`.

---

## 2. ECOSTRESS Collection 3 (`BuildID 0700`) Output Consistency Governance

Under NASA's ECOSTRESS Science Data System (SDS) operational governance, all software updates deployed into Collection 3 forward processing (`ECOv003-L3T-L4T-JET` $\rightarrow$ `JET3` $\rightarrow$ `AquaSEBS`, static **`BuildID 0700`**) **must produce primary science outputs numerically identical to the established late-2025 Collection 3 baseline**.

To satisfy this requirement in `v1.4.0`:
1. **Strict Numerical Invariance on Valid Pixels:** For every pixel where `(water == 1) & np.isfinite(WST_C)` is non-empty under default function arguments, `W_Wm2`, `LE_Wm2`, `LE_potential_Wm2`, and daylight-upscaled ET totals **must remain bit-for-bit identical** to `v1.3.0`.
2. **All New Physics Capabilities Must Be Opt-In or Diagnostic-Only:**
   * **Salinity Correction (`#11`):** `salinity_gL` defaults to `None` (freshwater baseline, $\sigma = 1.0$). Salinity scaling is applied only when a caller explicitly supplies `salinity_gL`.
   * **Dew-Point Formulation (`#12`):** Retain the existing linear psychrometric rule-of-thumb ($T_d = T_a - (100 - 100 \cdot \text{RH}) / 5.0$) as the default (`dewpoint_method="linear"` or `use_magnus_tetens=False`) so `Td_C`, `Tn`, `β`, `Te`, and `W_Wm2` do not shift in Collection 3 production; expose the Magnus–Tetens formulation (`dewpoint_method="magnus_tetens"`) as an explicit opt-in parameter.
   * **High-Wind Regime Flagging (`#13`):** Expose `"high_wind_flag"` (`windspeed_mps > 7.5`) as an additional boolean metadata key in the returned dictionary without altering or masking `LE_Wm2` values.
3. **Graceful Zero-Water Handling (`#14`):** Returning shape- and `RasterGeometry`-matched `NaN` arrays when zero valid water pixels exist restores terrestrial `L3T_JET`, `L3T_ETAUX`, `L4T_ESI`, and `L4T_WUE` granules on land-only and cloud-masked coastal tiles without altering a single terrestrial pixel value.

---

## 3. Historical Release Progression (`v1.0.1` – `v1.3.0`)

### 3.1 `v1.0.1` — Initial Release (2025-09-30)
* **Release Title:** `initial release`
* **Key Deliverables:**
  * Implemented `water_heat_flux()` (`AquaSEBS/water_heat_flux.py`) solving the 6-equation Equilibrium Temperature Model (ETM) for instantaneous water heat storage/release flux ($W$ / $G_{0w}$):
    $$T_n = 0.5(T_0 - T_d), \quad \eta = 0.35 + 0.015 T_0 + 0.0012 T_n^2, \quad S = 3.3 u$$
    $$\beta = 4.5 + 0.05 T_0 + (\eta + 0.47) S, \quad T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}, \quad W = \beta(T_e - T_0)$$
  * Implemented `AquaSEBS()` (`AquaSEBS/AquaSEBS.py`) coupling `verma_net_radiation` for broadband net radiation ($R_n$), `water_heat_flux()` for $W$, and `priestley_taylor` ($\alpha = 1.26$) for instantaneous open-water latent heat flux ($\text{LE}$).
  * Integrated automated static open-water masking via `NASADEM.swb` and automated meteorological/radiative forcing retrieval via `GEOS5FP`.
  * Established dictionary return schemas across both core functions and basic package import smoke tests (`tests/test_import_AquaSEBS.py`).

### 3.2 `v1.0.2` — Package Documentation Update (2025-10-08)
* **Release Title:** `documentation`
* **Key Deliverables:**
  * Updated package documentation and bumped metadata version in `pyproject.toml`.

### 3.3 `v1.1.0` — Daylight Evapotranspiration Upscaling (2025-11-11)
* **Release Title:** `daylight upscaling`
* **Key Deliverables:**
  * Added `upscale_to_daylight: bool = False` (configured via `UPSCALE_TO_DAYLIGHT` in `AquaSEBS/constants.py`).
  * Integrated `daylight_ET_from_instantaneous_LE()` from `daylight-evapotranspiration` to integrate instantaneous overpass $\text{LE}$ to daily/daylight totals using sinusoidal diurnal radiation scaling while preserving instantaneous evaporative fraction $\text{EF} = \text{LE} / (R_n - W)$.
  * Added `pytictoc.TicToc` execution timing and structured logging around diurnal upscaling.

### 3.4 `v1.2.0` — Initial Offline Mode Support (2026-01-21)
* **Release Title:** `offline mode`
* **Key Deliverables:**
  * Introduced `MissingOfflineParameter` exception class in `AquaSEBS/exceptions.py`.
  * Added `offline_mode: bool = False` parameter to `AquaSEBS()` and forwarded it to `verma_net_radiation()`, raising `MissingOfflineParameter` when `Ta_C` or `SWin_Wm2` were omitted in offline runs.

### 3.5 `v1.3.0` — `GEOS5FP_INPUTS` Constant (2026-01-21)
* **Release Title:** `GEOS5FP_INPUTS` constant
* **Key Deliverables:**
  * Added `GEOS5FP_INPUTS = ["Ta_C"]` to `AquaSEBS/constants.py` as an initial step toward tracking required assimilation inputs.

---

## 4. `v1.4.0` Development Plan (All Known Issues Consolidated)

All open issues (`#11`–`#16`) and post-`v1.3.0` engineering improvements are consolidated into **`v1.4.0`**.

### 4.1 Completed & Closed Engineering Milestones in `v1.4.0` (`#6`–`#10`)

| Workstream | Status | GitHub Issue | Summary of Changes |
| :--- | :---: | :---: | :--- |
| **Complete `offline_mode` & Lazy `GEOS5FP`** | Closed | `#6` | Added `offline_mode: bool = False` to `water_heat_flux()`, deferred `GEOS5FP()` instantiation via `_get_geos5fp()` helper in both `AquaSEBS()` and `water_heat_flux()`, gated `NASADEM.swb()` on `mask_non_water_pixels` and `offline_mode`, and raised `MissingOfflineParameter` when required inputs are omitted offline. |
| **Parameter Forwarding (`AquaSEBS` $\rightarrow$ `water_heat_flux`)** | Closed | `#7` | Added `SWin_Wm2` parameter to `water_heat_flux()` and forwarded `albedo`, `RH`, `SWin_Wm2`, and `offline_mode` from `AquaSEBS()` to eliminate duplicate `GEOS5FP` fetches and offline failures. |
| **Return Schema & Docstring Reconciliation** | Closed | `#8` | Standardized Unicode return keys (`"η"`, `"β"`) and documented `"G_Wm2"` across `water_heat_flux.py`, `AquaSEBS.py`, and `README.md`; fixed `water_heat_flux()` doctest output (`294.75` $\rightarrow$ `331.15`); removed obsolete `TODO` comments. |
| **Package Metadata & Dependency Hygiene** | Closed | `#9` | Updated `pyproject.toml` description for AquaSEBS, declared direct dependencies (`verma-net-radiation`, `pytictoc`), removed unused `ECOv002-*` entries, and deduplicated `pytest` in `[project.optional-dependencies] dev`. |
| **Comprehensive Unit Test Suite & CI Modernization** | Closed | `#10` | Added `tests/test_water_heat_flux.py`, `tests/test_AquaSEBS.py`, `tests/test_masking_and_upscaling.py`, and `tests/test_offline_mode.py`; standardized on Python `3.14` default and CI matrix `["3.11", "3.12", "3.13", "3.14"]`. |
| **Automated Version Synchronization & Exports** | Completed | — | Added `make set-version` target to `makefile`; replaced wildcard imports in `AquaSEBS/__init__.py` with explicit exports and `__all__`; removed dead `GEOS5FP_INPUTS` constant. |

---

### 4.2 Pending `v1.4.0` Implementation Tasks (`#11`–`#16`)

#### Task 1 (`#14`): Graceful Zero-Water-Pixel Handling (No `BlankOutputError` on Land-Only or Cloud-Masked Water Tiles)
* **Problem:** In ECOSTRESS Collection 3 DAG `EcostressDAG_Col3.47400.0802.20261115T235815542` (tile `52KBF`, `orbit=47400`, `scene=030`), `check_distribution(WST_C, "WST_C")` succeeds prior to masking because land temperatures are valid, but `WST_C = rt.where(water, WST_C, np.nan)` converts the entire raster to `NaN` when zero open-water pixels are unmasked. Subsequent `check_distribution()` calls on `Tn`, `η`, `S`, `β`, `Te`, `W_Wm2`, `epsilon`, and `LE_Wm2` raise `BlankOutputError`, aborting `JET3` (`JET3` Issue `#92`).
* **Implementation (`AquaSEBS/water_heat_flux.py` & `AquaSEBS/AquaSEBS.py`):**
  * After masking `WST_C` in `water_heat_flux()`, evaluate `has_valid_water = bool(np.any(np.isfinite(np.asarray(WST_C))))`.
  * When `not has_valid_water`: log a warning and return the standard output dictionary populated with shape- and `RasterGeometry`-matched `NaN` arrays (`rt.full_like(WST_C, np.nan)` for `Raster`, `np.full_like(..., np.nan)` for `ndarray`) without calling `check_distribution()` on derived all-`NaN` arrays.
  * In `AquaSEBS()`, check `has_valid_water = bool(np.any(np.isfinite(np.asarray(W_Wm2))))`. When `False`, bypass `check_distribution()` on `W_Wm2`/`LE_Wm2` and skip `daylight_ET_from_instantaneous_LE()`, populating `LE_Wm2`, `LE_potential_Wm2`, `G_Wm2`, and (when `upscale_to_daylight=True`) daylight outputs with `NaN` arrays.
* **Collection 3 Consistency:** Zero numerical change whenever at least one finite water pixel exists.
* **Unit Tests:** Add `test_all_land_or_cloud_masked_water_returns_nan_without_blank_output_error()` in `tests/test_masking_and_upscaling.py` covering both `water = np.zeros(...)` and `water = np.ones(...)` with all-`NaN` `WST_C` under `upscale_to_daylight=True` and `False`.

#### Task 2 (`#15`): Update `GEOS5FP.wind_speed()` to `GEOS5FP.wind_speed_mps()`
* **Problem:** `AquaSEBS/water_heat_flux.py` calls `_get_geos5fp().wind_speed(...)`, which was removed in `GEOS5FP >= 2.24.0` in favor of `.wind_speed_mps(...)`.
* **Implementation (`AquaSEBS/water_heat_flux.py` & `tests/test_offline_mode.py`):**
  * Replace `_get_geos5fp().wind_speed(...)` with `_get_geos5fp().wind_speed_mps(...)` (with `getattr(..., "wind_speed_mps", None) or getattr(..., "wind_speed")` fallback if needed for legacy `GEOS5FP` compatibility).
  * Update `FailingGEOS5FP` in `tests/test_offline_mode.py` to define `wind_speed_mps`.
* **Collection 3 Consistency:** `JET3` already passes `windspeed_mps` retrieved via `GEOS5FP.wind_speed_mps(...)`; standalone online calls now resolve identical wind-speed grids.

#### Task 3 (`#11`): Opt-In Turk (1970) Water Salinity Correction Factor (`salinity_gL=None`)
* **Scientific Reference:** Abdelrady et al. (2016), Section 2.2.5, Eqs. 19–20; Turk (1970).
* **Implementation (`AquaSEBS/AquaSEBS.py`):**
  * Add optional `salinity_gL: Union[Raster, np.ndarray, float] = None` parameter to `AquaSEBS()`.
  * **Default (`salinity_gL is None`):** Perform no salinity scaling (freshwater potential evaporation, exact Collection 3 `BuildID 0700` behavior).
  * **When `salinity_gL is not None`:** Compute Turk's (1970) empirical salinity reduction factor:
    $$\sigma = \text{clip}\!\left(1.025 - 0.0246 \exp(0.00879 \cdot S),\; 0.0,\; 1.0\right)$$
    Scale `LE_Wm2` and `LE_potential_Wm2` by $\sigma$ prior to optional daylight upscaling (so daylight ET automatically inherits the salinity scaling), and include `"salinity_reduction_factor": σ` in the returned dictionary.
* **Unit Tests:** Add tests in `tests/test_AquaSEBS.py` verifying: (a) `salinity_gL=None` produces identical results to `v1.3.0`, (b) freshwater (`0 g/L`) yields $\sigma \approx 1.0$, (c) open ocean (`34.7 g/L`) reduces `LE_Wm2` by $\sim 0.8\%$, and (d) Great Salt Lake (`240 g/L`) reduces `LE_Wm2` by $\sim 22.8\%$.

#### Task 4 (`#12`): Opt-In Magnus–Tetens Saturation Vapor Pressure Dew-Point Formula (`dewpoint_method="linear"`)
* **Scientific Reference:** WMO / Alduchov & Eskridge (1996) Magnus–Tetens psychrometric formulation.
* **Implementation (`AquaSEBS/constants.py`, `AquaSEBS/water_heat_flux.py`, `AquaSEBS/AquaSEBS.py`):**
  * Define `DEWPOINT_METHOD = "linear"` in `AquaSEBS/constants.py` and add `dewpoint_method: str = DEWPOINT_METHOD` to `water_heat_flux()` and `AquaSEBS()`.
  * **Default (`dewpoint_method="linear"`):** Preserve the existing linear approximation $T_d = T_a - (100 - 100 \cdot \text{RH}) / 5.0$ so Collection 3 `BuildID 0700` outputs remain 100% invariant.
  * **Opt-In (`dewpoint_method="magnus_tetens"`):** Compute dew-point temperature using Magnus–Tetens ($a = 17.625$, $b = 243.04^\circ\text{C}$):
    $$\gamma_m = \ln(\text{clip}(\text{RH}, 10^{-6}, 1.0)) + \frac{a \cdot T_a}{b + T_a}, \quad T_d = \frac{b \cdot \gamma_m}{a - \gamma_m}$$
  * Raise `ValueError` if an unrecognized `dewpoint_method` is passed.
* **Unit Tests:** Add tests in `tests/test_water_heat_flux.py` verifying default `"linear"` parity, `"magnus_tetens"` psychrometric accuracy across humid ($\text{RH}=0.80$), moderate ($\text{RH}=0.50$), and arid ($\text{RH}=0.15$) regimes, and safe boundary clamping at $\text{RH}=0.0$.

#### Task 5 (`#13`): Non-Intrusive High-Wind Regime Diagnostic Flagging (`high_wind_flag`)
* **Scientific Reference:** Fisher et al. (2023), Abstract, Results, and Table 2 ($u > 7.5\text{ m}\cdot\text{s}^{-1}$ regime shift).
* **Implementation (`AquaSEBS/constants.py`, `AquaSEBS/water_heat_flux.py`, `AquaSEBS/AquaSEBS.py`):**
  * Define `HIGH_WIND_THRESHOLD_MPS = 7.5` in `AquaSEBS/constants.py`.
  * Whenever `windspeed_mps` is available in `water_heat_flux()` / `AquaSEBS()`, compute `high_wind_flag = windspeed_mps > HIGH_WIND_THRESHOLD_MPS` (preserving `Raster` geometry when inputs are `Raster`s) and include `"high_wind_flag"` in the output dictionary.
  * **Collection 3 Consistency:** Diagnostic metadata only; does not modify or mask `W_Wm2`, `LE_Wm2`, or daylight ET outputs.
* **Unit Tests:** Add tests in `tests/test_water_heat_flux.py` and `tests/test_AquaSEBS.py` verifying `"high_wind_flag"` is `False` for $u \le 7.5\text{ m}\cdot\text{s}^{-1}$ and `True` for $u > 7.5\text{ m}\cdot\text{s}^{-1}$.

#### Task 6 (`#16`): Fix Escaped LaTeX Control Characters & Document New `v1.4.0` Options in `README.md`
* **Implementation (`README.md`):**
  * Restore clean LaTeX math formatting (`$\alpha = 1.26$`, `$\beta$`, `\text{...}`, `\rho_w`, `\frac{...}`, `\approx`, `\right)`) across `README.md` Sections 1.1, 1.4, 2.2, 2.3, 3.1, 3.2, 4.2, 4.3, 5, and 6, and fix `Zolt{\'a}n` in the BibTeX citation block.
  * Document `salinity_gL`, `dewpoint_method`, `"high_wind_flag"`, and zero-water `NaN` behavior in the API Reference and Release History sections of `README.md`.

---

### 4.3 Release Verification & Execution Checklist (`v1.4.0`)

1. **Remove Temporary Local Scripts:** Ensure `create_v1_4_0_issues.fish` is deleted so `git status` shows only tracked source/test/doc changes.
2. **Verify Editable Install & Full Pytest + Doctest Suite on Python 3.14:**
   ```fish
   pip install -e . --no-deps
   make test
   pytest --doctest-modules AquaSEBS