#!/usr/bin/env fish

# Create pending GitHub issues for AquaSEBS v1.4.0
# Preserves ECOSTRESS Collection 3 (BuildID 0700) numerical consistency

set -l repo "JPL-Evapotranspiration-Algorithms/AquaSEBS"

echo "Publishing pending v1.4.0 issues to $repo..."

# Issue 1: Zero-water-pixel BlankOutputError handling
printf '%s\n' \
    '### Summary' \
    'Target release: **`v1.4.0`**' \
    '' \
    'When `mask_non_water_pixels=True` is applied to a land-only scene (`water == 0` everywhere) or a scene where all open-water pixels are cloud-masked (`(water == 1) & np.isfinite(WST_C)` has zero finite elements), `water_heat_flux()` and `AquaSEBS()` raise `check_distribution.BlankOutputError: variable Tn is a blank image`, aborting callers such as `JET3` (`JET3` Issue #92).' \
    '' \
    '### Details' \
    '- In ECOSTRESS Collection 3 forward-processing DAG `EcostressDAG_Col3.47400.0802.20261115T235815542` (reproduced on tile `52KBF`, `orbit=47400`, `scene=030`), 19 runs across 17 MGRS tiles aborted inside `AquaSEBS` after all terrestrial ET models succeeded, discarding valid terrestrial ET retrievals across 90–94% of each tile.' \
    '- `check_distribution(WST_C, "WST_C")` passes prior to masking because land temperatures are valid, but `WST_C = rt.where(water, WST_C, np.nan)` converts the entire raster to `NaN`, causing subsequent `check_distribution()` calls on `Tn`, `η`, `S`, `β`, `Te`, `W_Wm2`, `epsilon`, and `LE_Wm2` to fail.' \
    '- Under ECOSTRESS Collection 3 `BuildID 0700` governance, remediation must preserve bit-for-bit numerical output parity on all valid water and land pixels.' \
    '' \
    '### Proposed Fix (`v1.4.0`)' \
    '- In `water_heat_flux()` (`AquaSEBS/water_heat_flux.py`), check `bool(np.any(np.isfinite(np.asarray(WST_C))))` after applying the water mask. When `False`, log a warning and return the standard dictionary schema populated with shape- and `RasterGeometry`-matched `NaN` arrays without calling `check_distribution()` on all-`NaN` outputs.' \
    '- In `AquaSEBS()` (`AquaSEBS/AquaSEBS.py`), skip `check_distribution()` on all-`NaN` `W_Wm2` / `LE_Wm2` outputs and skip `daylight_ET_from_instantaneous_LE()` when no valid water pixels exist (populating daylight outputs with `NaN` arrays if `upscale_to_daylight=True`).' \
    '- Add regression tests in `tests/test_masking_and_upscaling.py` covering both land-only (`water == 0`) and cloud-masked water (`water == 1` with all-`NaN` `WST_C`) scenes under `upscale_to_daylight=True` and `False`.' \
| gh issue create --repo "$repo" \
    --title "Handle zero valid water pixels gracefully without raising BlankOutputError (v1.4.0)" \
    --body-file -

# Issue 2: GEOS5FP.wind_speed() -> wind_speed_mps()
printf '%s\n' \
    '### Summary' \
    'Target release: **`v1.4.0`**' \
    '' \
    'In `AquaSEBS/water_heat_flux.py` (line 293), missing `windspeed_mps` values in online mode (`offline_mode=False`) are fetched via `_get_geos5fp().wind_speed(...)`, which was removed in `GEOS5FP >= 2.24.0` in favor of `.wind_speed_mps(...)`.' \
    '' \
    '### Details' \
    '- In `GEOS5FP >= 2.24.0` (pinned at `GEOS5FP~=2.27.0` in `ECOv003-L3T-L4T-JET`), `.wind_speed(...)` was renamed to `.wind_speed_mps(...)` for SI unit standardization.' \
    '- While `JET3` masks this defect by pre-fetching `wind_speed_mps` externally, any standalone `AquaSEBS()` or `water_heat_flux()` call with `windspeed_mps=None` and `offline_mode=False` raises `AttributeError: '\''GEOS5FPConnection'\'' object has no attribute '\''wind_speed'\''`.' \
    '' \
    '### Proposed Fix (`v1.4.0`)' \
    '- Update `_get_geos5fp().wind_speed(...)` to `_get_geos5fp().wind_speed_mps(...)` in `AquaSEBS/water_heat_flux.py`.' \
    '- Update the mock `FailingGEOS5FP` helper in `tests/test_offline_mode.py` to define `wind_speed_mps`.' \
| gh issue create --repo "$repo" \
    --title "Replace removed GEOS5FP.wind_speed() call with GEOS5FP.wind_speed_mps() (v1.4.0)" \
    --body-file -

# Issue 3: Escaped LaTeX formatting cleanup in README.md
printf '%s\n' \
    '### Summary' \
    'Target release: **`v1.4.0`**' \
    '' \
    'Several LaTeX equations, symbols, and a BibTeX author entry in `README.md` contain literal control characters or missing backslashes that break Markdown math rendering.' \
    '' \
    '### Details' \
    '- Affected tokens include `$ lpha = 1.26$`, `$ eta$`, `$	ext{W}$`, `$
ho_w$`, `$rac{...}$`, ` pprox`, `
ight)`, and `Zolt{'"'"'a}n` in the BibTeX citation block across Sections 1.1, 1.4, 2.2, 2.3, 3.1, 3.2, 4.2, 4.3, 5, and 6.' \
    '' \
    '### Proposed Fix (`v1.4.0`)' \
    '- Restore properly escaped LaTeX macros (`$\alpha = 1.26$`, `$\beta$`, `\text{...}`, `\rho_w`, `\frac{...}`, `\approx`, `\right)`, and `Zolt{\\'"'"'a}n`) across `README.md`.' \
| gh issue create --repo "$repo" \
    --title "Fix escaped LaTeX control characters in README.md (v1.4.0)" \
    --body-file -

echo "Done. Current open issues:"
make show-issues