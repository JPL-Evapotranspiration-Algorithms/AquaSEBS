# AquaSEBS Comprehensive Development Plan & Release Roadmap

**Package:** `AquaSEBS` (Aquatic Surface Energy Balance System)  
**Repository:** [JPL-Evapotranspiration-Algorithms/AquaSEBS](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS)  
**Author & Maintainer:** Gregory H. Halverson ([gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)), NASA Jet Propulsion Laboratory  
**Current Target Release:** `v1.4.0`  
**Default Development Runtime:** Python `3.14` (CI support matrix: Python `3.11` – `3.14`)  

---

## 1. Executive Summary

`AquaSEBS` implements the open-water surface energy balance and evaporation model combining the Equilibrium Temperature Model (ETM) of Edinger et al. (1968) and Abdelrady et al. (2016) with the energy-limited Priestley–Taylor (1972) evaporation parameterization validated globally across 19 open-water sites by Fisher et al. (2023) for NASA's **ECOSTRESS** (`L3_ET_PT-JPL` Collection 2+) and **OpenET** systems.

This document establishes the complete engineering and scientific development plan across:
1. **Completed Historical Releases (`v1.0.1` – `v1.3.0`)**: Foundational physics implementation, diurnal upscaling, initial `offline_mode` hooks, and input configuration constants.
2. **Current Release Cycle (`v1.4.0` — Technical & Infrastructure Consolidation)**: Full `offline_mode` enforcement, lazy `GEOS5FP` initialization, complete parameter forwarding (`RH`, `SWin_Wm2`, `albedo`), top-level namespace cleanup, comprehensive `pytest` unit test coverage, Python `3.14` environment standardization, and automated version management (`make set-version`).
3. **Pending Release Finalization Checklist (`v1.4.0`)**: Remaining verification, formatting cleanup, issue closure (`#6`, `#7`), and GitHub/PyPI release execution.
4. **Shelved Scientific Roadmap (`v1.5.0+`)**: Well-documented, scoped scientific enhancements (`#11`, `#12`, `#13`) deferred from `v1.4.0` to preserve numerical stability and clear release boundaries.

---

## 2. Historical Release Progression (`v1.0.1` – `v1.3.0`)

### 2.1 `v1.0.1` — Initial Release (2025-09-30)
* **Release Title:** `initial release`
* **Key Deliverables:**
  * Implemented `water_heat_flux()` (`AquaSEBS/water_heat_flux.py`) solving the 6-equation Equilibrium Temperature Model (ETM) for instantaneous water heat storage/release flux ($W$ / $G_{0w}$):
    $$T_n = 0.5(T_0 - T_d), \quad \eta = 0.35 + 0.015 T_0 + 0.0012 T_n^2, \quad S = 3.3 u$$
    $$\beta = 4.5 + 0.05 T_0 + (\eta + 0.47) S, \quad T_e = T_d + \frac{\text{SW}_{\text{net}}}{\beta}, \quad W = \beta(T_e - T_0)$$
  * Implemented `AquaSEBS()` (`AquaSEBS/AquaSEBS.py`) coupling `verma_net_radiation` for broadband net radiation ($R_n$), `water_heat_flux()` for $W$, and `priestley_taylor` ($\alpha = 1.26$) for instantaneous open-water latent heat flux ($\text{LE}$).
  * Integrated automated static open-water masking via `NASADEM.swb` and automated meteorological/radiative forcing retrieval via `GEOS5FP`.
  * Established dictionary return schemas across both core functions and basic package import smoke tests (`tests/test_import_AquaSEBS.py`).

### 2.2 `v1.0.2` — Package Documentation Update (2025-10-08)
* **Release Title:** `documentation`
* **Key Deliverables:**
  * Updated package documentation and bumped metadata version in `pyproject.toml`.

### 2.3 `v1.1.0` — Daylight Evapotranspiration Upscaling (2025-11-11)
* **Release Title:** `daylight upscaling`
* **Key Deliverables:**
  * Added `upscale_to_daylight: bool = False` (configured via `UPSCALE_TO_DAYLIGHT` in `AquaSEBS/constants.py`).
  * Integrated `daylight_ET_from_instantaneous_LE()` from `daylight-evapotranspiration` to integrate instantaneous overpass $\text{LE}$ to daily/daylight totals using sinusoidal diurnal radiation scaling while preserving instantaneous evaporative fraction $\text{EF} = \text{LE} / (R_n - W)$.
  * Added `pytictoc.TicToc` execution timing and structured logging around diurnal upscaling.

### 2.4 `v1.2.0` — Initial Offline Mode Support (2026-01-21)
* **Release Title:** `offline mode`
* **Key Deliverables:**
  * Introduced `MissingOfflineParameter` exception class in `AquaSEBS/exceptions.py`.
  * Added `offline_mode: bool = False` parameter to `AquaSEBS()` and forwarded it to `verma_net_radiation()`, raising `MissingOfflineParameter` when `Ta_C` or `SWin_Wm2` were omitted in offline runs.

### 2.5 `v1.3.0` — `GEOS5FP_INPUTS` Constant (2026-01-21)
* **Release Title:** `GEOS5FP_INPUTS` constant
* **Key Deliverables:**
  * Added `GEOS5FP_INPUTS = ["Ta_C"]` to `AquaSEBS/constants.py` as an initial step toward tracking required assimilation inputs.

---

## 3. `v1.4.0` Development Plan (Consolidating All Post-`v1.3.0` Work)

All engineering work since `v1.3.0` targets **`v1.4.0`** directly (superseding the intermediate `1.3.1` metadata bump). Scientific enhancements (`#11`, `#12`, `#13`) are intentionally excluded from `v1.4.0` so that this release focuses strictly on runtime reliability, offline execution guarantees, parameter forwarding, test coverage, and tooling consistency.

### 3.1 Completed Engineering Milestones in `v1.4.0`

| Workstream | Status | GitHub Issue | Summary of Changes |
| :--- | :---: | :---: | :--- |
| **Package Metadata & Dependency Hygiene** | Completed | `#9` (Closed) | Updated `pyproject.toml` description for AquaSEBS, declared direct dependencies (`verma-net-radiation`, `pytictoc`), removed unused `ECOv002-*` entries, and deduplicated `pytest` in `[project.optional-dependencies] dev`. |
| **Return Schema & Docstring Reconciliation** | Completed | `#8` (Closed) | Standardized Unicode return keys (`"η"`, `"β"`) and documented `"G_Wm2"` across `water_heat_flux.py`, `AquaSEBS.py`, and `README.md`; fixed `water_heat_flux()` doctest output (`294.75` $\rightarrow$ `331.15`); removed obsolete `TODO` comments. |
| **Comprehensive Unit Test Suite** | Completed | `#10` (Closed) | Added `tests/test_water_heat_flux.py`, `tests/test_AquaSEBS.py`, `tests/test_masking_and_upscaling.py`, and `tests/test_offline_mode.py` verifying analytical ETM equations, sign behavior, Priestley–Taylor energy balance, `verma_net_radiation` integration, water masking, and daylight upscaling. |
| **Automated Version Synchronization** | Completed | — | Added `make set-version` target to `makefile` updating `pyproject.toml`, `CITATION.cff` (`version`, `preferred-citation.version`, `date-released`), and `README.md` in one command. |
| **Python 3.14 & CI Modernization** | Completed | `#10` (Closed) | Standardized default development environment on Python `3.14` (`makefile`, `README.md`, `python-publish.yml`), updated GitHub Actions to `actions/checkout@v4` and `actions/setup-python@v5`, and configured CI test matrix across Python `["3.11", "3.12", "3.13", "3.14"]`. |
| **Complete `offline_mode` & Lazy `GEOS5FP`** | Implemented | `#6` (Open) | Added `offline_mode: bool = False` to `water_heat_flux()`, deferred `GEOS5FP()` instantiation via `_get_geos5fp()` helper in both `AquaSEBS()` and `water_heat_flux()`, gated `NASADEM.swb()` on `mask_non_water_pixels` and `offline_mode`, and raised `MissingOfflineParameter` when required inputs are omitted offline. |
| **Parameter Forwarding (`AquaSEBS` $\rightarrow$ `water_heat_flux`)** | Implemented | `#7` (Open) | Added `SWin_Wm2` parameter to `water_heat_flux()` and forwarded `albedo`, `RH`, `SWin_Wm2`, and `offline_mode` from `AquaSEBS()` to eliminate duplicate `GEOS5FP` fetches and offline failures. |
| **Top-Level Namespace & Constants Cleanup** | Implemented | — | Replaced wildcard imports in `AquaSEBS/__init__.py` with explicit exports and `__all__`, removed dead `GEOS5FP_INPUTS` from `AquaSEBS/constants.py`, and removed temporary root utility scripts. |

---

### 3.2 In-Progress & Pending Tasks to Finalize `v1.4.0`

The following items remain to be completed before tagging and publishing `v1.4.0`:

#### Task 1: Fix Escaped LaTeX Formatting Typo in `README.md`
* **Observation:** During the `README.md` update, several LaTeX macros had literal control characters or missing backslashes (e.g., `$ lpha = 1.26$`, `$ eta$`, `$	ext{W}$`, `$ho_w$`, `$rac{...}$`, `Zolt{'a}n`).
* **Action:** Restore clean LaTeX math formatting (`$\alpha = 1.26$`, `$\beta$`, `\text{...}`, `\rho_w`, `\frac{...}`, `\approx`, `\right)`) across `README.md` Sections 1.1, 1.4, 2.2, 2.3, 3.1, 3.2, 4.2, 4.3, 5, and 6, as well as `Zolt{\'a}n` in the BibTeX citation block.

#### Task 2: Verify Editable Install & Full Pytest Suite on Python 3.14
* **Action:**
  ```fish
  pip install -e . --no-deps
  make test
  pytest --doctest-modules AquaSEBS
  ```
* **Acceptance Criteria:**
  * `AquaSEBS.__version__ == "1.4.0"`.
  * All unit tests in `tests/` (`test_import_AquaSEBS.py`, `test_water_heat_flux.py`, `test_AquaSEBS.py`, `test_masking_and_upscaling.py`, `test_offline_mode.py`) pass with zero errors.
  * Doctests in `AquaSEBS/AquaSEBS.py` and `AquaSEBS/water_heat_flux.py` pass cleanly.

#### Task 3: Close Technical GitHub Issues `#6` and `#7`
* **Action:**
  ```fish
  gh issue close 6 --comment "Resolved in v1.4.0: Completed offline_mode support across AquaSEBS() and water_heat_flux(), gated NASADEM.swb() on mask_non_water_pixels and offline_mode, and deferred GEOS5FP() instantiation until a network retrieval is needed."
  gh issue close 7 --comment "Resolved in v1.4.0: Added SWin_Wm2 support to water_heat_flux() and forwarded albedo, RH, SWin_Wm2, and offline_mode from AquaSEBS()."
  ```
* **Acceptance Criteria:** `make show-issues` lists only the three shelved scientific issues (`#11`, `#12`, `#13`).

#### Task 4: Commit, Tag, and Publish `v1.4.0`
* **Action:**
  1. Review `git diff` and `git status` to confirm no untracked artifacts remain.
  2. Commit all `v1.4.0` changes to `main`.
  3. Create GitHub Release `v1.4.0` (triggering `.github/workflows/python-publish.yml` to build and publish `AquaSEBS 1.4.0` to PyPI on Python `3.14`), or publish via `make dist`.

---

## 4. Shelved Scientific Roadmap (`v1.5.0+`)

The following three scientific enhancements are intentionally shelved from `v1.4.0` and tracked in GitHub Issues for future scientific release cycles:

### 4.1 Issue #11 — Native Water Salinity Correction Factor (`salinity_gL`)
* **Scientific Reference:** Abdelrady et al. (2016), Section 2.2.5, Eqs. 19–20; Turk (1970).
* **Problem Statement:** `AquaSEBS` currently models freshwater potential evaporation ($S = 0\text{ g}\cdot\text{L}^{-1}$), overestimating evaporation by $\sim 1\%$ over open ocean ($S \approx 34.7\text{ g}\cdot\text{L}^{-1}$) and up to $\sim 27\%$ over hypersaline lakes such as the Great Salt Lake ($S \approx 240\text{ g}\cdot\text{L}^{-1}$).
* **Planned Design:**
  * Add optional `salinity_gL: Union[Raster, np.ndarray, float] = None` parameter to `AquaSEBS()`.
  * Compute Turk's (1970) empirical salinity reduction factor:
    $$\sigma = \text{clip}\!\left(1.025 - 0.0246 \exp(0.00879 \cdot S),\; 0.0,\; 1.0\right)$$
  * Scale `LE_Wm2` (and downstream daylight ET outputs when `upscale_to_daylight=True`) by $\sigma$ when `salinity_gL` is provided, and return `"salinity_reduction_factor"` in the output dictionary.
  * Add analytical and hypersaline benchmark unit tests based on Table 1 and Table 8 of Abdelrady et al. (2016).

### 4.2 Issue #12 — Magnus–Tetens Saturation Vapor Pressure Dew-Point Formula
* **Scientific Reference:** WMO / Alduchov & Eskridge (1996) Magnus–Tetens psychrometric formulation.
* **Problem Statement:** When `Td_C` is omitted, `water_heat_flux()` currently uses the linear rule-of-thumb $T_d \approx T_a - (100 - 100 \cdot \text{RH}) / 5.0$, which diverges significantly when $\text{RH} < 50\%$ (common across arid Western U.S. reservoirs and dryland inland seas).
* **Planned Design:**
  * Replace the linear approximation in `water_heat_flux()` with the Magnus–Tetens formulation ($a = 17.625$, $b = 243.04^\circ\text{C}$):
    $$\gamma_m = \ln(\text{clip}(\text{RH}, 10^{-6}, 1.0)) + \frac{a \cdot T_a}{b + T_a}, \quad T_d = \frac{b \cdot \gamma_m}{a - \gamma_m}$$
  * Add unit tests verifying accuracy across humid ($\text{RH} = 80\%$), moderate ($\text{RH} = 50\%$), and arid ($\text{RH} = 15\%$) regimes as well as boundary clamping at $\text{RH} \rightarrow 0$.

### 4.3 Issue #13 — High-Wind Regime Diagnostic Flagging ($u > 7.5\text{ m}\cdot\text{s}^{-1}$)
* **Scientific Reference:** Fisher et al. (2023), Abstract, Results, and Table 2.
* **Problem Statement:** Across 19 global validation sites, Fisher et al. (2023) showed that instantaneous AquaSEBS retrievals degrade when daily mean wind speed exceeds $7.5\text{ m}\cdot\text{s}^{-1}$ ($r^2$ drops from $0.71$ to $0.47$, RMSE increases from $38\%$ to $62\%$ of mean) as evaporation shifts from radiatively controlled to aerodynamically controlled.
* **Planned Design:**
  * Define `HIGH_WIND_THRESHOLD_MPS = 7.5` in `AquaSEBS/constants.py`.
  * Expose a boolean diagnostic mask `"high_wind_flag"` (`windspeed_mps > HIGH_WIND_THRESHOLD_MPS`) in the output dictionary whenever `windspeed_mps` is resolved, allowing downstream consumers (e.g., ECOSTRESS / OpenET quality-control pipelines) to filter or weight instantaneous scenes.

---

## 5. Summary Timeline & Milestone Matrix

| Version | Date | Focus Area | Key Outcome |
| :--- | :--- | :--- | :--- |
| `v1.0.1` | 2025-09-30 | Core Physics | Initial ETM + Priestley–Taylor + `GEOS5FP` + `NASADEM.swb` release |
| `v1.0.2` | 2025-10-08 | Documentation | Initial package metadata & README updates |
| `v1.1.0` | 2025-11-11 | Diurnal Scaling | Added `upscale_to_daylight` and `daylight_ET_from_instantaneous_LE` |
| `v1.2.0` | 2026-01-21 | Offline Guardrails | Added `offline_mode` and `MissingOfflineParameter` to `AquaSEBS()` |
| `v1.3.0` | 2026-01-21 | Configuration | Added `GEOS5FP_INPUTS` constant |
| **`v1.4.0`** | **2026-11-25** | **Technical & CI Consolidation** | **Complete `offline_mode`, lazy `GEOS5FP`, parameter forwarding (`#6`, `#7`), clean exports, full `pytest` suite (`#8`–`#10`), Python `3.14` default & `3.11`–`3.14` CI** |
| `v1.5.0` *(Planned)* | TBD | Scientific Enhancements | Salinity correction (`#11`), Magnus–Tetens $T_d$ (`#12`), high-wind flag (`#13`) |