# AquaSEBS Repository Release and Commit History

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS`  
**Author:** Gregory H. Halverson (`gregory.h.halverson@jpl.nasa.gov`)  
**Current Version:** `v1.3.0` (plus uncommitted package hygiene updates)

---

## 1. Release Timeline Summary

| Release Tag | Release Date (UTC) | Release Title | Core Milestone |
| :--- | :--- | :--- | :--- |
| **`v1.0.1`** | 2025-09-30T22:07:17Z | initial release | Initial implementation of Equilibrium Temperature Model (`water_heat_flux`), Verma net radiation + Priestley–Taylor energy balance (`AquaSEBS`), `NASADEM.swb` water masking, and dictionary return schema |
| **`v1.0.2`** | 2025-10-08T17:14:34Z | documentation | Package documentation updates (version bump in `pyproject.toml`) |
| **`v1.1.0`** | 2025-11-11T19:49:53Z | daylight upscaling | Added `upscale_to_daylight` flag, `daylight_ET_from_instantaneous_LE()` integration, and `TicToc` execution timing |
| **`v1.2.0`** | 2026-01-21T21:14:40Z | offline mode | Added `offline_mode` parameter and `MissingOfflineParameter` exception for `Ta_C`, `SWin_Wm2`, and `verma_net_radiation()` |
| **`v1.3.0`** | 2026-01-21T21:55:38Z | `GEOS5FP_INPUTS` constant | Added `GEOS5FP_INPUTS = ["Ta_C"]` to `AquaSEBS/constants.py` |

---

## 2. Detailed Commit Progression

### 2.1 Initial Development & First Release (`v1.0.1`)

* **`2d94ffc` (2025-09-29) — `initial commit`**
  * Created the initial package layout (`AquaSEBS/__init__.py`, `AquaSEBS.py`, `constants.py`, `version.py`, `water_heat_flux.py`, `pyproject.toml`, and import tests in `tests/`).
  * `AquaSEBS.py` began as a PT-JPL prototype containing `LE_soil`, `Rn_canopy`, `PET_water`, and `PET_land` stubs, while `pyproject.toml` carried over the PT-JPL description.
  * `water_heat_flux.py` implemented the ETM sequence returning scalar/array `W_Wm2`, with an initial algebraic typo in Eq. 12 (`Te = Td_C + SWnet * β`).

* **`86c1dc6` (2025-09-30) — `water heat flux working`**
  * Fixed the equilibrium temperature formula in `water_heat_flux.py` to `Te = Td_C + (SWnet / β)`.
  * Expanded docstrings and inline equations with explicit citations to Abdelrady et al. (2016, Section 2.2.1, Eqs. 8–13) and Fisher et al. (2023, Methods).

* **`97a225d` (2025-09-30) — `water surface energy balance`**
  * Replaced the PT-JPL prototype in `AquaSEBS.py` with the open-water energy balance pipeline: optional net radiation retrieval via `verma_net_radiation()`, water heat flux calculation via `water_heat_flux()`, and Priestley–Taylor latent heat calculation (`LE_Wm2 = PT_ALPHA * epsilon * (Rn_Wm2 - W_Wm2)`).

* **`1e01c24` (2025-09-30) — `water surface latent heat flux working`**
  * Added automatic `geometry` inference from `WST_C` when `WST_C` is a `Raster`, default `GEOS5FP()` connection initialization, automated `Ta_C` retrieval, distribution checks (`check_distribution`), and corrected `ST_C=WST_C` in `verma_net_radiation()`.

* **`c00aba5` (2025-09-30) — `removing lower and upper bounds`**
  * Removed the empirical `rt.clip(W_Wm2, 0, W_MAX_PROPORTION * Rn_Wm2)` constraint in `AquaSEBS.py` so water heat flux relies directly on the unclipped ETM physics.

* **`c741599` (2025-09-30) — `using priestley-taylor package`**
  * Delegated latent heat flux computation to `priestley_taylor.priestley_taylor()` (`priestley-taylor>=1.1.2`), added `α` and `γ_Pa` parameters, guarded `Ta_C` and `Rn_Wm2` lookups on `geometry` and `time_UTC`, and added `NASADEM` to `pyproject.toml`.

* **`cdc0a79` & `8b04d74` (2025-09-30) — `v1.0.0` / `v1.0.1 initial release`**
  * Integrated automated surface water masking (`NASADEM.swb(geometry=geometry)` and `mask_non_water_pixels = True`).
  * Transitioned both `water_heat_flux()` and `AquaSEBS()` from returning a single output array/raster to returning comprehensive dictionaries containing all inputs, ETM intermediates (`Tn`, `η`, `S`, `β`, `Te`), `W_Wm2`, radiation terms, and `LE_Wm2`.
  * Tagged `v1.0.1` for initial PyPI publication.

### 2.2 Documentation Release (`v1.0.1` $\rightarrow$ `v1.0.2`)

* **`66e771b` (2025-10-08) — `v1.0.2 documentation`**
  * Bumped package version to `1.0.2` in `pyproject.toml` to publish updated package documentation (`README.md`).

### 2.3 Daylight Upscaling (`v1.0.2` $\rightarrow$ `v1.1.0`)

* **`ec24834` (2025-11-11) — `daylight upscaling`**
  * Added `daylight-evapotranspiration` to `pyproject.toml` dependencies and `UPSCALE_TO_DAYLIGHT = False` to `AquaSEBS/constants.py`.
  * Added `upscale_to_daylight: bool = UPSCALE_TO_DAYLIGHT` to `AquaSEBS()`, calling `daylight_ET_from_instantaneous_LE()` with `TicToc` timing and logging when `time_UTC` is provided, and merging daylight ET products into `results`.

### 2.4 Offline Mode (`v1.1.0` $\rightarrow$ `v1.2.0`)

* **`52be1fe` (2026-01-21) — `v1.2.0 offline mode`**
  * Created `AquaSEBS/exceptions.py` defining `MissingOfflineParameter`.
  * Added `offline_mode: bool = False` to `AquaSEBS()`, raising `MissingOfflineParameter` if `Ta_C` or `SWin_Wm2` are omitted when `offline_mode=True`, and forwarding `offline_mode` to `verma_net_radiation()`.

### 2.5 Constants Update (`v1.2.0` $\rightarrow$ `v1.3.0`)

* **`11cab51` (2026-01-21) — `v1.3.0 GEOS5FP_INPUTS constant`**
  * Added `GEOS5FP_INPUTS = ["Ta_C"]` to `AquaSEBS/constants.py` and bumped version to `1.3.0`.

---

## 3. Current Working Tree (Post-`v1.3.0`)

* **Uncommitted changes in `pyproject.toml`:**
  * Explicitly configured `build-backend = "setuptools.build_meta"`.
  * Corrected package `description` to `"Aquatic Surface Energy Balance System (AquaSEBS) Water Surface Evaporation Model"`.
  * Added missing direct dependencies `pytictoc` and `verma-net-radiation`, removed unused `ECOv002-CMR` and `ECOv002-granules`, deduplicated `pytest` under `dev`, and updated `requires-python = ">=3.11"`.