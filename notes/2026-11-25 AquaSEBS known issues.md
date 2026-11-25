# AquaSEBS Technical Known Issues

**Repository:** `JPL-Evapotranspiration-Algorithms/AquaSEBS` (v1.3.0)

This document tracks software engineering, API consistency, packaging, and test coverage issues identified in the current `AquaSEBS` codebase. (Scientific and methodological limitations are documented separately in `AquaSEBS_science_caveats.md`.)

## 1. Incomplete `offline_mode` Support & Eager Network Initialization

* **`water_heat_flux()` does not support `offline_mode`:**
  * `AquaSEBS()` accepts `offline_mode: bool = False`, validates required inputs (`Ta_C`, `SWin_Wm2`) before querying `GEOS5FP`, forwards `offline_mode` to `verma_net_radiation()`, and raises `MissingOfflineParameter` when required offline parameters are absent.
  * `water_heat_flux()` lacks an `offline_mode` parameter entirely. When `W_Wm2` is not explicitly provided, `water_heat_flux()` will attempt live `NASADEM.swb()` and `GEOS5FP` network queries for `water`, `Ta_C`, `RH`, `windspeed_mps`, `SWin`, and `albedo` rather than raising `MissingOfflineParameter`.
* **Unconditional `NASADEM.swb()` call in `AquaSEBS()`:**
  * In `AquaSEBS/AquaSEBS.py`, when `water is None and geometry is not None`, `NASADEM.swb(geometry=geometry)` is invoked without checking `offline_mode` (or `mask_non_water_pixels`), causing an unexpected external lookup during offline execution.
* **Eager `GEOS5FP()` instantiation:**
  * Both `AquaSEBS()` and `water_heat_flux()` instantiate `GEOS5FP_connection = GEOS5FP()` immediately on entry whenever `GEOS5FP_connection is None`, even if all meteorological and radiative inputs are already supplied or `offline_mode=True`. Instantiation should be deferred until a remote data lookup is actually required.

## 2. Parameter Forwarding & Signature Inconsistencies

* **`SWin_Wm2` and `albedo` are not forwarded to `water_heat_flux()`:**
  * `AquaSEBS()` accepts `SWin_Wm2` and `albedo` (and retrieves `SWin_Wm2` from `GEOS5FP` if needed) for `verma_net_radiation()`, but calls `water_heat_flux()` with `SWnet=SWnet` (defaulting to `None`) without passing `SWin_Wm2` or `albedo`. Consequently, `water_heat_flux()` redundantly re-queries `GEOS5FP` for `SWin` and `ALBEDO` whenever `SWnet` is not explicitly supplied.
* **`RH` is not forwarded to `water_heat_flux()`:**
  * `AquaSEBS()` accepts `RH` and forwards it to `verma_net_radiation()`, but omits `RH` when calling `water_heat_flux()`. If `Td_C` is not provided, `water_heat_flux()` redundantly re-fetches `RH` from `GEOS5FP`.
* **Redundant `None` guards:**
  * `AquaSEBS/AquaSEBS.py` checks `WST_C is not None` inside the `Rn_Wm2 is None` conditional block even though `WST_C is None` already raises `ValueError` earlier in the function.
  * `AquaSEBS/water_heat_flux.py` re-checks `geometry is not None and time_UTC is not None` inside an outer `if Td_C is None and geometry is not None and time_UTC is not None:` block.

## 3. Return Schema & Documentation Mismatches

* **Unicode vs. ASCII return dictionary keys:**
  * `water_heat_flux.py` returns Unicode Greek keys (`"η"`, `"β"`), whereas `README.md` documents ASCII keys (`"eta"`, `"beta"`). Downstream consumers expecting ASCII keys will encounter `KeyError`s unless both keys are returned or the schema and documentation are standardized.
* **Return type docstring mismatch in `water_heat_flux()`:**
  * The docstring in `AquaSEBS/water_heat_flux.py` states `:return: Water heat flux in Watts per square meter.`, implying a scalar or array return, whereas the function signature and implementation return `Dict[str, Union[Raster, np.ndarray, float]]`.
* **Brittle `priestley_taylor` dictionary key lookup:**
  * `AquaSEBS/AquaSEBS.py` directly indexes `priestley_taylor_results["LE_potential_Wm2"]` without a fallback check (e.g., for `"LE_Wm2"`), risking a `KeyError` if upstream `priestley-taylor` output keys vary across versions.
* **Stale `TODO` comment:**
  * `AquaSEBS/AquaSEBS.py` retains `## TODO use NASADEM surface water body extent to mask out land when processing on rasters` despite `NASADEM.swb()` masking already being implemented.

## 4. Packaging & Dependency Hygiene (`pyproject.toml`)

* **Incorrect package description:**
  * `pyproject.toml` currently defines `description = "Priestley-Taylor Jet Propulsion Laboratory Evapotranspiration Model"` (copied from `PTJPL`) instead of describing AquaSEBS.
* **Undeclared direct dependencies vs. unused dependencies:**
  * `verma-net-radiation` and `pytictoc` are imported directly in `AquaSEBS/AquaSEBS.py` but are missing from `pyproject.toml` `dependencies` (relying on transitive installation via `ECOv002-granules`).
  * `ECOv002-CMR` and `ECOv002-granules` are listed in `pyproject.toml` `dependencies` and `tests/test_import_dependencies.py` despite not being imported anywhere inside the `AquaSEBS/` package.
* **Duplicate optional dependency:**
  * `"pytest"` is listed twice under `[project.optional-dependencies] dev` (`"pytest>=6.0"` and `"pytest"`).

## 5. Testing & CI Coverage Gaps

* **Lack of functional and numerical unit tests:**
  * `tests/` currently contains only import smoke tests (`test_import_AquaSEBS.py` and `test_import_dependencies.py`). There are no automated tests verifying:
    * The Equilibrium Temperature Model equations ($T_n, \eta, S, \beta, T_e, W$).
    * `offline_mode=True` execution and `MissingOfflineParameter` exception handling.
    * Water masking behavior (`mask_non_water_pixels=True/False`).
    * Daylight evapotranspiration upscaling (`upscale_to_daylight=True`).
* **Outdated GitHub Actions versions and Python matrix:**
  * `.github/workflows/ci.yml` tests Python `3.10`–`3.13` (missing `3.14`, which is now configured as the default in `makefile`) and uses legacy actions (`actions/checkout@v3`, `actions/setup-python@v4`).
  * `.github/workflows/python-publish.yml` uses `actions/setup-python@v3`.

## Recommended Remediation Priority

1. **Propagate `offline_mode` and lazy `GEOS5FP` initialization** across both `AquaSEBS()` and `water_heat_flux()`, and forward `SWin_Wm2`, `albedo`, and `RH` from `AquaSEBS()` to `water_heat_flux()`.
2. **Add an offline unit test suite** covering ETM math, `MissingOfflineParameter`, water masking, and daylight upscaling so CI validates model behavior without network dependencies.
3. **Clean up `pyproject.toml` and CI workflows** (fix package description, align declared dependencies with actual imports, deduplicate `pytest`, and update GitHub Actions / Python `3.14` matrix).
4. **Reconcile return keys (`eta`/`η`, `beta`/`β`) and docstrings** between `water_heat_flux.py` and `README.md`.