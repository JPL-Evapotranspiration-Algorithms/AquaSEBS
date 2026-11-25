# AquaSEBS Comprehensive Development Plan & Release Roadmap

**Package:** `AquaSEBS` (Aquatic Surface Energy Balance System)  
**Repository:** [JPL-Evapotranspiration-Algorithms/AquaSEBS](https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS)  
**Author & Maintainer:** Gregory H. Halverson ([gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)), NASA Jet Propulsion Laboratory  
**Current Target Release:** `v1.4.0`  
**Operational Baseline Compatibility:** NASA ECOSTRESS Collection 3 (`BuildID 0700`)  
**Default Development Runtime:** Python `3.14` (CI support matrix: Python `3.11` – `3.14`)  

---

## 1. Executive Summary

`AquaSEBS` implements the open-water surface energy balance and evaporation model combining the Equilibrium Temperature Model (ETM) of Edinger et al. (1968) and Abdelrady et al. (2016) with the energy-limited Priestley–Taylor (1972) evaporation parameterization validated globally across 19 open-water sites by Fisher et al. (2023) for NASA's **ECOSTRESS** (`L3_ET_PT-JPL` Collection 2+ / `L3T_JET`) and **OpenET** systems.

All post-`v1.3.0` engineering improvements and open GitHub issues (`#6`–`#16`) are consolidated into **`v1.4.0`** while guaranteeing that **default science outputs remain bit-for-bit identical to ECOSTRESS Collection 3 (`BuildID 0700`)**:
- **Operational Production Fixes (`#14`, `#15`)**: Graceful all-`NaN` handling when scenes contain zero valid water pixels (eliminating `BlankOutputError` aborts in `JET3` Issue `#92`) and updating `GEOS5FP.wind_speed()` to `GEOS5FP.wind_speed_mps()`.
- **Opt-In / Non-Breaking Scientific Enhancements (`#11`, `#12`, `#13`)**: Native Turk (1970) salinity correction (`salinity_gL=None` by default), opt-in Magnus–Tetens dew-point calculation (`dewpoint_method="linear"` by default), and non-intrusive high-wind diagnostic flagging (`"high_wind_flag"`).
- **Documentation & Release Finalization (`#16`)**: Escaped LaTeX cleanup in `README.md`, full `pytest` and doctest verification, closure of `#11`–`#16`, and PyPI/GitHub release of `v1.4.0`.
