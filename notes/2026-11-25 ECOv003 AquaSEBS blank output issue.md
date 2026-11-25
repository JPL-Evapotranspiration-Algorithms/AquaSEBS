# `AquaSEBS` Operational Integration, Output Consistency Governance, and Edge-Case Analysis for ECOSTRESS Collection 3

- **Package:** [`AquaSEBS`](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS) (Aquatic Surface Energy Balance System)
- **Author & Maintainer:** Gregory H. Halverson ([gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)), NASA Jet Propulsion Laboratory
- **Report Date:** November 25, 2026
- **Context & Upstream Systems:**
  - `ECOv003-L3T-L4T-JET` (`v1.23.0` $\rightarrow$ `v1.24.0`, `BuildID 0700`)
  - `JET3` (`v3.3.0` $\rightarrow$ `v3.4.0`)
  - `GEOS5FP` (`v2.24.0` – `v2.27.0`)
  - Container Runtime: `eco-level-2-3-4-application:v1.14.0` (`linux/amd64`, Python 3.12)
  - Operational DAG Triage Source: `EcostressDAG_Col3.47400.0802.20261115T235815542` & tile `52KBF` reproduction (`orbit=47400`, `scene=030`, `2026-11-14`)

---

## 1. Executive Summary

`AquaSEBS` serves as the operational open-water surface energy balance and evaporation engine within NASA's **ECOSTRESS Collection 3** Level-3/Level-4 evapotranspiration processing chain (`ECOv003-L3T-L4T-JET` $\rightarrow$ `JET3` $\rightarrow$ `AquaSEBS`). While terrestrial pixels across Sentinel-2 MGRS tiles ($1568 \times 1568$ pixels at $70\text{ m}$) are modeled by a four-model land evapotranspiration ensemble (`PT-JPL-SM`, `STIC-JPL`, `PM-JPL`, and `BESS-JPL`), open-water pixels (`water == 1`) are partitioned separately using `AquaSEBS` so inland lakes, reservoirs, estuaries, and coastal waters receive physically consistent latent heat flux ($\text{LE}$) and daylight evaporation ($\text{ET}_{\text{daily}}$) estimates rather than being masked out as missing data.

Operational triage of **762 forward-processing logs** from ECOSTRESS Collection 3 DAG `EcostressDAG_Col3.47400.0802.20261115T235815542` (November 24, 2026) and deterministic reproduction on tile `52KBF` (`orbit=47400`, `scene=030`, `2026-11-14`) exposed a critical operational edge case—alongside two coupled cross-package defects:

1. **Zero-Water-Pixel `BlankOutputError` Crash (`AquaSEBS` / `JET3` Issue `#92`):** When a tile contains zero open-water pixels (`water == 0` everywhere) or all open-water pixels are obscured by clouds (`(water == 1) & np.isfinite(WST_C)` has zero finite elements), `water_heat_flux()` masks non-water pixels to `NaN` (`WST_C = rt.where(water, WST_C, np.nan)`) and subsequently invokes `check_distribution(Tn, "Tn")` on an all-`NaN` array. `check_distribution` raises `BlankOutputError: variable Tn is a blank image`, aborting the entire `JET()` execution and discarding valid terrestrial ET retrievals across 90–94% of the tile.
2. **`GEOS5FPConnection.wind_speed()` API Mismatch (`AquaSEBS/water_heat_flux.py:293`):** When `windspeed_mps` is omitted and `offline_mode=False`, `water_heat_flux()` calls `_get_geos5fp().wind_speed(...)`, which was removed in `GEOS5FP >= 2.24.0` in favor of `.wind_speed_mps(...)`.
3. **Uncaught `check_distribution.BlankOutputError` at the PGE Entrypoint (`ECOv003-L3T-L4T-JET` Issue `#91`):** `ECOv003_L3T_L4T_JET.py` catches `JET3.BlankOutput`, whereas `check_distribution` raises `check_distribution.BlankOutputError`, causing blank-scene aborts to exit with generic status `1` (`run_pge::docker_status= 1`) instead of the designated `ECOv003_exit_codes.BLANK_OUTPUT` skip code.

Because ECOSTRESS Collection 3 Science Data System (SDS) governance requires all maintained software updates under **`BuildID 0700`** to remain **scientifically and numerically consistent** with the deployed late-2025 Collection 3 baseline, remediation must preserve exact floating-point output parity on all valid water and land pixels while eliminating spurious aborts on water-free or cloud-masked water scenes.

---

## 2. Role of `AquaSEBS` in ECOSTRESS Collection 3

### 2.1 Architectural Position in the ECOSTRESS L3T/L4T Pipeline

In the ECOSTRESS Collection 3 processing architecture, `AquaSEBS` is invoked as a specialized surface-type component inside the unified `JET3` science engine, which is orchestrated in production by the `ECOv003-L3T-L4T-JET` Product Generating Executable (PGE):

```text
ECOv003 L2T_LSTE (ST_K, emissivity, water, cloud, elevation_km)
        +
ECOv003 L2T_STARS (NDVI, albedo)
        +
GEOS-5 FP Reanalysis (Ta_C, RH, SM, SWin, wind_speed_mps, AOT, COT, vapor, ozone)
        │
        ▼
┌────────────────────────────────────────────────────────────────────────────┐
│  JET3 Core Science Engine (JET3/JET.py)                                    │
│                                                                            │
│  1. Radiative Transfer & Surface Energy Balance                            │
│     ├─ FLiES-ANN (downwelling shortwave SWin_Wm2, PAR, NIR)                │
│     ├─ BESS-JPL (coupled photosynthesis GPP & latent heat LE_BESS_Wm2)     │
│     └─ Verma Net Radiation (broadband Rn_Wm2)                              │
│                                                                            │
│  2. Terrestrial Evapotranspiration Ensemble (water_mask == 0)              │
│     ├─ PT-JPL-SM  (soil-moisture constrained Priestley–Taylor)             │
│     ├─ STIC-JPL   (surface temperature initiated closure)                  │
│     ├─ PM-JPL     (Penman–Monteith formulation)                            │
│     └─ BESS-JPL   (coupled biophysical conductance)                        │
│        └─► Ensemble Median (ETdaily) & Spread (ETuncertainty)              │
│                                                                            │
│  3. Open-Water Surface Energy Balance (water_mask == 1)                    │
│     └─ AquaSEBS (Equilibrium Temperature Model W_Wm2 + Priestley–Taylor)   │
│        └─► Replaces/populates open-water pixels in L3T_JET                 │
└────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
ECOSTRESS Collection 3 Tiled Products (70 m MGRS, 1568 × 1568 pixels)
  ├─ L3T_JET   (PTJPLSMinst/daily, STICJPLinst/daily, BESSJPLdaily, PMJPLdaily,
  │             ETdaily, ETuncertainty, partitioning fractions, cloud, water)
  ├─ L3T_ETAUX (Ta, RH, SM, Rg, Rn, cloud, water)
  ├─ L4T_ESI   (ESI, PET, cloud, water)
  └─ L4T_WUE   (WUE, GPP, cloud, water)
```

### 2.2 Why `AquaSEBS` Is Required Alongside Terrestrial ET Models

Terrestrial evapotranspiration models (`PT-JPL-SM`, `STIC-JPL`, `PM-JPL`, `BESS-JPL`) rely on vegetation indices (`NDVI`, `fAPAR`, `LAI`), canopy stomatal conductance, and soil moisture (`SM`)—quantities that are undefined or physically invalid over open water (`water == 1`). Conversely, open water bodies exhibit high thermal inertia and vertical convective mixing, storing substantial solar energy during daytime ($W > 0$) and releasing stored heat at night or during cold-air outbreaks ($W < 0$).

- **Collection 1 Behavior:** Open water bodies were masked out entirely as missing data (`NaN`).
- **Collection 2 & Collection 3 Behavior:** Following global validation across 19 open-water sites in Fisher et al. (2023), `AquaSEBS` was integrated into the operational ECOSTRESS processing chain to provide continuous, physically based open-water evaporation over lakes, reservoirs, rivers, and coastal waters identified by the static surface water body mask (`water == 1`).
- **Forcing & Parameter Coupling in `JET3`:** In `JET3/JET.py`, `AquaSEBS` receives pre-computed, tile-aligned rasters from the upstream pipeline:
  - `WST_C = ST_C` (water skin temperature from `L2T_LSTE`)
  - `emissivity` (from `L2T_LSTE`) and `albedo` (from `L2T_STARS`)
  - `Ta_C` and `RH` (downscaled from GEOS-5 FP)
  - `windspeed_mps` (retrieved via `GEOS5FP_connection.wind_speed_mps(...)`)
  - `SWnet` and `SWin_Wm2` (from `FLiES-ANN`)
  - `Rn_Wm2` (from Verma net radiation)
  - `water = water_mask` (from `L2T_LSTE`)
  - `geometry` and `time_UTC`

---

## 3. ECOSTRESS Collection 3 Output Consistency Requirements

### 3.1 Science Data System (SDS) Operational Baseline Governance

Within NASA's ECOSTRESS Science Data System (SDS), Collection 3 forward processing is governed by strict algorithmic and data-product consistency rules across the multi-year mission record:

1. **Static Product `BuildID` (`0700`):** All Collection 3 L3T/L4T output granule filenames (e.g., `ECOv003_L3T_JET_47400_030_52KBF_20261114T060122_0700_01.zip`) carry the fixed Collection 3 **`BuildID 0700`** (`DEFAULT_BUILD = "0700"`), while exact package provenance is recorded in each granule's JSON sidecar under `StandardMetadata.PGEVersion`.
2. **Cross-Epoch Output Consistency:** Collection 3 forward processing began in late 2025 (`eco-level-2-3-4-application` `v1.7.0` / `v1.7.1`) and ran continuously for ~11 months until the November 2, 2026 decommissioning of Suomi NPP VIIRS (`VNP43NRT` $\rightarrow$ JPSS `VCD43NRT`). Under SDS operational requirements, any maintenance release deployed during Collection 3 forward processing **must produce outputs scientifically and numerically consistent with the established Collection 3 baseline** so time-series analyses across the archive do not suffer artificial step-changes or algorithmic discontinuities.
3. **Implications for `AquaSEBS` Changes:**
   - **No uncalibrated physics changes in patch/minor maintenance releases:** Enhancements that alter numerical values over valid open-water pixels—such as salinity scaling (`Issue #11`) or replacing the linear dew-point approximation with Magnus–Tetens (`Issue #12`)—must remain shelved for a future major/collection-level science release rather than being introduced silently into Collection 3 `BuildID 0700`.
   - **Strict numerical invariance on valid water pixels:** Any fix for the zero-water-pixel crash (`Issue #92`) must leave `W_Wm2`, `LE_Wm2`, and daylight-upscaled evaporation **bit-for-bit identical** on every pixel where `(water == 1) & np.isfinite(WST_C)` is non-empty.
   - **Coverage restoration without terrestrial side effects:** Guarding or gracefully handling empty water masks restores valid terrestrial `L3T_JET`, `L3T_ETAUX`, `L4T_ESI`, and `L4T_WUE` granules on land-only and cloud-masked coastal tiles without altering a single terrestrial pixel value.

---

## 4. Detailed Root-Cause Analysis of Observed Edge Cases

### 4.1 Edge Case 1: `BlankOutputError: variable Tn is a blank image` on Zero Valid Water Pixels (`Issue #92`)

#### Operational Manifestation
In forward-processing DAG `EcostressDAG_Col3.47400.0802.20261115T235815542` (triaged November 24, 2026), **19 granule runs across 17 MGRS tiles** (`50JPQ`, `50KRC`, `51KVB`, `51KVV`, `51KWU`, `51KYT`, `51KYU`, `51KZB`, `51KZT`, `51LYC`, `51LZC`, `52KBD`, `52KBF`, `52KBG`, `52KCE`, `52KCF`, `52KCG`) aborted with the exact traceback reproduced locally via `reproduce_aquasebs_blank.py` on tile `52KBF` (`orbit=47400`, `scene=030`, `2026-11-14`):

```text
[2026-11-24 18:23:20 INFO] checking input distributions for AquaSEBS
[2026-11-24 18:23:20 INFO] variable ST_C min: 15.730 mean: 45.982 max: 63.877 nan: 4.09% (nan)
...
[2026-11-24 18:23:21 INFO] running AquaSEBS
[2026-11-24 18:23:21 INFO] variable WST_C min: 15.730 mean: 45.982 max: 63.877 nan: 4.09% (nan)
...
[2026-11-24 18:23:24 INFO] variable Tn (float32) with [nan] unique values
[2026-11-24 18:23:24 INFO] * nan: 2458624
Traceback (most recent call last):
  File ".../JET3/JET.py", line 933, in JET
    AquaSEBS_results = AquaSEBS(
  File ".../AquaSEBS/AquaSEBS.py", line 120, in AquaSEBS
    water_heat_flux_results = water_heat_flux(
  File ".../AquaSEBS/water_heat_flux.py", line 153, in water_heat_flux
    check_distribution(Tn, "Tn")
  File ".../check_distribution/check_distribution.py", line 89, in check_distribution
    raise BlankOutputError(f"variable {variable}{date_message}{target_message} is a blank image")
check_distribution.check_distribution.BlankOutputError: variable Tn is a blank image
```

#### Execution Trace & Mechanism
1. **Upstream terrestrial models succeed completely:** On tile `52KBF` ($1568 \times 1568 = 2,458,624$ pixels), `ST_C` is valid over $95.91\%$ of the tile (`nan: 4.09%`). `FLiES-ANN` (`32.96s`), `BESS-JPL` (`54.49s`), `Verma Rn` (`0.67s`), `STIC-JPL` (`7.84s`), `PT-JPL-SM` (`4.41s`), and `PM-JPL` (`6.69s`) all execute to completion and produce valid terrestrial ET estimates across ~90–94% of the scene.
2. **Unmasked input check passes (`AquaSEBS.py:216` & `water_heat_flux.py:232`):** Both `AquaSEBS()` and `water_heat_flux()` call `check_distribution(WST_C, "WST_C")` **before** applying `water` masking. Because land pixels have valid temperatures (`15.73 °C` to `63.88 °C`), this check succeeds.
3. **Water masking converts `WST_C` to 100% `NaN` (`water_heat_flux.py:249`):**
   ```python
   if mask_non_water_pixels:
       ...
       WST_C = rt.where(water, WST_C, np.nan)
   ```
   When the tile has zero open-water pixels (`np.sum(water == 1) == 0`) or when every open-water pixel falls under cloud cover (`np.sum((water == 1) & np.isfinite(WST_C)) == 0`), `WST_C` becomes an array of $2,458,624$ `NaN` values.
4. **Unconditional `check_distribution` on derived ETM variables raises `BlankOutputError`:**
   - `water_heat_flux.py` computes `Tn = 0.5 * (WST_C - Td_C)` (all-`NaN`) and immediately calls `check_distribution(Tn, "Tn")`, triggering `BlankOutputError`.
   - Even if `Tn` were skipped, subsequent calls to `check_distribution` on `η`, `S`, `β`, `Te`, `W_Wm2` (`water_heat_flux.py:329–366`) and on `W_Wm2`, `epsilon`, `LE_Wm2` (`AquaSEBS.py:285, 297, 314`) would each raise `BlankOutputError` whenever the masked output array is all-`NaN`.

#### Why This Must Be Handled in Both `AquaSEBS` and `JET3`
- **In `AquaSEBS` (`AquaSEBS.py` & `water_heat_flux.py`):** Raising `BlankOutputError` when `mask_non_water_pixels=True` produces zero valid water pixels is a semantic defect in `AquaSEBS` itself. An all-`NaN` output raster is the **mathematically expected result** when a valid scene contains no unmasked open water. Any caller (`JET3`, OpenET, or standalone research scripts) passing a valid land-only or cloud-masked coastal image should receive all-`NaN` flux arrays matching `WST_C`'s shape/geometry rather than an unhandled crash.
- **In `JET3` (`JET3/JET.py:933`):** Guarding `AquaSEBS(...)` with `if np.any((water_mask == 1) & np.isfinite(ST_C)):` provides belt-and-suspenders protection and avoids unnecessary `GEOS5FP` wind-speed queries (`GEOS5FP_connection.wind_speed_mps(...)`) on tiles that contain no valid water pixels.

---

### 4.2 Edge Case 2: Removed `GEOS5FPConnection.wind_speed()` Method (`water_heat_flux.py:293`)

In `AquaSEBS/water_heat_flux.py` (line 293) and `tests/test_offline_mode.py` (line 161), wind speed retrieval still invokes the legacy `GEOS5FP` method name:

```python
if windspeed_mps is None:
    if offline_mode:
        raise MissingOfflineParameter("offline mode is enabled but windspeed_mps is not provided")
    if geometry is not None and time_UTC is not None:
        windspeed_mps = _get_geos5fp().wind_speed(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )
```

In `GEOS5FP >= 2.24.0` (currently pinned at `GEOS5FP~=2.27.0` in `ECOv003-L3T-L4T-JET`), `.wind_speed(...)` was renamed to `.wind_speed_mps(...)` as part of SI unit standardization. In `JET3`, `windspeed_mps` is pre-fetched via `GEOS5FP_connection.wind_speed_mps(...)` and passed explicitly into `AquaSEBS()`, masking this bug in `JET3` runs. However, any standalone `AquaSEBS()` or `water_heat_flux()` call with `offline_mode=False` and `windspeed_mps=None` will raise `AttributeError: 'GEOS5FPConnection' object has no attribute 'wind_speed'` unless updated to `.wind_speed_mps(...)` (with a `getattr` fallback if backward compatibility with older `GEOS5FP` versions is desired).

---

### 4.3 Edge Case 3: Uncaught `check_distribution.BlankOutputError` at the PGE Exit Boundary (`Issue #91`)

In `ECOv003_L3T_L4T_JET/ECOv003_L3T_L4T_JET.py` (line 214), the PGE entrypoint catches `JET3.BlankOutput`, whereas `check_distribution.check_distribution(...)` raises `check_distribution.BlankOutputError`. When a scene is genuinely blank (e.g., 100% cloud-covered `ST_C` or `NDVI` checked in `read_ECOv003_inputs.py`), `BlankOutputError` escapes `L3T_L4T_JET()` uncaught and causes the container to exit with status `1` (`run_pge::docker_status= 1`) instead of `ECOv003_exit_codes.BLANK_OUTPUT`. Catching `(BlankOutput, BlankOutputError)` in `ECOv003_L3T_L4T_JET.py` guarantees proper SDS skip-code reporting.

---

## 5. Recommended Remediation Plan (Preserving Collection 3 Output Consistency)

To resolve these edge cases cleanly without altering a single numerical output on valid water or terrestrial pixels under ECOSTRESS Collection 3 (`BuildID 0700`):

### 5.1 `AquaSEBS` (`AquaSEBS/water_heat_flux.py` & `AquaSEBS/AquaSEBS.py`)

1. **Graceful Zero-Water-Pixel Handling in `water_heat_flux()` (`AquaSEBS/water_heat_flux.py`):**
   - After applying `WST_C = rt.where(water, WST_C, np.nan)` when `mask_non_water_pixels=True`, check `has_valid_water = bool(np.any(np.isfinite(np.asarray(WST_C))))`.
   - If `not has_valid_water`:
     - Log an informational/warning message (`"no valid unmasked water pixels present; returning NaN water heat flux outputs"`).
     - Construct `nan_output = rt.full_like(WST_C, np.nan)` (if `isinstance(WST_C, Raster)`) or `np.full_like(np.asarray(WST_C, dtype=np.float32), np.nan)` so spatial `RasterGeometry` and dtype metadata are preserved.
     - Return the standard dictionary schema (`WST_C`, `Ta_C`, `RH`, `Td_C`, `windspeed_mps`, `SWnet`, `Tn`, `η`, `S`, `β`, `Te`, `W_Wm2`) populated with `nan_output` for all spatial output fields without triggering `BlankOutputError`.
2. **Graceful Zero-Water-Pixel Handling in `AquaSEBS()` (`AquaSEBS/AquaSEBS.py`):**
   - Check `has_valid_water = bool(np.any(np.isfinite(np.asarray(W_Wm2))))` (and `np.any(np.isfinite(np.asarray(WST_C)))`).
   - When `not has_valid_water`, skip `check_distribution` on all-`NaN` water outputs (`W_Wm2`, `LE_Wm2`) and skip `daylight_ET_from_instantaneous_LE` (populating `LE_Wm2`, `LE_potential_Wm2`, `G_Wm2`, and—when `upscale_to_daylight=True`—`ET_daylight_kg` / `LE_daylight_Wm2` with shape/geometry-matched `NaN` arrays).
   - When `has_valid_water` is `True`, execute the **exact existing equations and `check_distribution` checks unchanged**, guaranteeing 100% numerical identity with Collection 3 `BuildID 0700`.
3. **Update `GEOS5FP` Wind Speed Call (`AquaSEBS/water_heat_flux.py:293`):**
   - Replace `_get_geos5fp().wind_speed(...)` with `_get_geos5fp().wind_speed_mps(...)` (and update `FailingGEOS5FP` in `tests/test_offline_mode.py`).
4. **Add Unit Test Coverage (`tests/test_masking_and_upscaling.py`):**
   - Add `test_all_land_or_cloud_masked_water_returns_nan_without_blank_output_error()` verifying both `water = np.zeros(...)` (land-only tile) and `water = np.ones(...)` with all-`NaN` `WST_C` (cloud-masked water tile), with `upscale_to_daylight=True` and `upscale_to_daylight=False`.

### 5.2 Upstream `JET3` & `ECOv003-L3T-L4T-JET` Hardening

| Repository | File | Change |
| :--- | :--- | :--- |
| **`JET3`** (`v3.4.0`) | `JET3/JET.py` (line 933) | Guard `AquaSEBS(...)` with `if np.any((water_mask == 1) & np.isfinite(ST_C)):` so open-water processing and `GEOS5FP` wind-speed queries are bypassed when no valid water pixels exist, leaving terrestrial ET outputs intact. |
| **`ECOv003-L3T-L4T-JET`** (`v1.24.0`) | `ECOv003_L3T_L4T_JET/ECOv003_L3T_L4T_JET.py` | Catch `(BlankOutput, BlankOutputError)` so genuinely blank input granules exit with `ECOv003_exit_codes.BLANK_OUTPUT` (`#91`), and pin `JET3~=3.4.0` (`#92`). |

---

## 6. Summary Matrix of Issues & Architectural Boundaries

| Issue / Defect | Package & Location | Trigger Condition | Impact Without Fix | Resolution & Collection 3 Parity |
| :--- | :--- | :--- | :--- | :--- |
| **Zero-Water `BlankOutputError` (`#92`)** | `AquaSEBS/water_heat_flux.py:153` & `JET3/JET.py:933` | Tile has no water pixels (`water == 0`) or all water pixels are cloud-masked (`WST_C` is `NaN` over `water == 1`) | Aborts entire `JET()` run; loses all valid terrestrial ET granules (19 runs / 17 tiles in `DAG 47400`) | Return all-`NaN` spatial outputs in `AquaSEBS` and guard `AquaSEBS()` call in `JET3`; **zero change** on valid water or land pixels |
| **`GEOS5FP.wind_speed()` Removal** | `AquaSEBS/water_heat_flux.py:293` | Standalone `AquaSEBS()` / `water_heat_flux()` call with `windspeed_mps=None` and `GEOS5FP >= 2.24.0` | `AttributeError: 'GEOS5FPConnection' object has no attribute 'wind_speed'` | Switch to `.wind_speed_mps(...)` |
| **Uncaught `BlankOutputError` (`#91`)** | `ECOv003_L3T_L4T_JET.py:214` | Genuinely blank input scene (e.g., 100% cloud-covered `ST_C` or `NDVI`) | PGE exits with status `1` instead of `BLANK_OUTPUT` | Catch `(BlankOutput, BlankOutputError)` and return `BLANK_OUTPUT` |
| **Salinity & Dew-Point Upgrades (`#11`, `#12`)** | `AquaSEBS` (`v1.5.0+` roadmap) | Marine/hypersaline water bodies or arid $\text{RH} < 50\%$ | Would shift numerical outputs relative to Collection 3 `BuildID 0700` baseline | **Keep shelved** until next Collection/major baseline release |