# HBEF Watershed 3 — TOPMODEL Explorer

A Shiny web application for running TOPMODEL, a rainfall-runoff model that predicts streamflow based on watershed topography, soil properties, and climate inputs. Built for the **Hubbard Brook Experimental Forest (HBEF) Watershed 3** as part of the EDS Capstone project.

**Authors:** Samuel Handel, Kaelyn Harvey, Max Hughes

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Before Running the App: Two-Step Precompute](#️-before-running-the-app-two-step-precompute)
  - [Step 1 — precompute_twi.R](#step-1--precompute_twir)
  - [Step 2 — precompute_all.R](#step-2--precompute_allr)
- [Running the App](#running-the-app)
- [Using the App](#using-the-app)
- [Data Requirements](#data-requirements)

---

## Overview

This app simulates streamflow for HBEF Watershed 3 using TOPMODEL, a topography-driven hydrological model. It allows users to:

- Manually tune model parameters via interactive sliders
- Automatically calibrate parameters using the Nelder-Mead optimizer over a user-defined calibration window
- Load a custom catchment by uploading your own time series and TOPIDX files
- Visualize simulated vs. observed discharge, flow duration curves, and input forcing data
- Explore an animated **saturation map** showing which parts of the watershed are predicted to be saturated on any given day — click directly on the discharge plot to jump the map to that date
- View the Topographic Wetness Index (TWI) map and histogram with summary statistics
- Export all model output to a CSV file

---

## Project Structure

```
project/
├── app.R                    # Main Shiny application
├── precompute_twi.R         # Step 1: compute TWI from raw DEM via WhiteboxTools
├── precompute_all.R         # Step 2: convert spatial outputs to .rds files
├── Data/
│   ├── dem_10m.tif                          # 10-meter resolution DEM (input)
│   ├── dem_ws_clipped.tif                   # DEM clipped to watershed (produced by Step 1)
│   ├── Watershed3HB.shp                     # Watershed boundary shapefile (+ .dbf, .prj, .shx)
│   ├── DailyWatershed.csv                   # Daily precipitation (W3)
│   ├── HBEF_DailyStreamflow_1956-2024.csv   # Observed daily streamflow
│   ├── HBEF_air_temp_daily.csv              # Daily air temperature
│   ├── twi_precomputed.tif                  # Produced by Step 1
│   ├── twi_matrix.rds                       # Produced by Step 2
│   ├── twi_extent.rds                       # Produced by Step 2
│   ├── dem_matrix.rds                       # Produced by Step 2
│   ├── dem_extent.rds                       # Produced by Step 2
│   ├── ws_boundary.rds                      # Produced by Step 2
│   ├── topidx.rds                           # Produced by Step 2
│   └── twi_lam.rds                          # Produced by Step 2
└── README.md
```

---

## Prerequisites

Install the following R packages before running anything:

```r
install.packages(c(
  "shiny", "bslib", "plotly", "dplyr", "lubridate",  # app dependencies
  "whitebox", "terra", "sf",                           # Step 1 (precompute_twi.R)
  "raster"                                             # Step 2 (precompute_all.R)
))
```

WhiteboxTools (WBT) also requires a one-time binary installation after the R package is installed:

```r
whitebox::install_whitebox()
```

---

## Before Running the App: Two-Step Precompute

**Both steps are required and must be completed in order before launching the app for the first time.**

The app loads spatial data from lightweight `.rds` files at startup rather than reading rasters or shapefiles directly. This means the deployed app on shinyapps.io requires **no spatial packages** (`terra`, `sf`, `raster`, etc.) at runtime, which avoids binary dependency issues and reduces startup time.

The two precompute scripts produce all seven `.rds` files the app needs. Run them once, then you can launch or re-deploy the app as many times as you like without re-running them.

---

### Step 1 — `precompute_twi.R`

**What it does:** Computes the Topographic Wetness Index (TWI) from the raw 10-meter DEM using WhiteboxTools (WBT) and saves the result as `Data/twi_precomputed.tif`.

**Run it:**

```r
source("precompute_twi.R")
```

You will see progress messages as it works through each step:

```
Loading DEM and watershed ...
Breaching depressions ...
FD8 flow accumulation ...
Computing slope ...
Computing TWI ...
Done! Saved to: Data/twi_precomputed.tif
You can now run precompute_all.R next.
```

**What it produces:** `Data/twi_precomputed.tif`

**Processing steps inside the script:**
1. Loads the 10m DEM and clips it to the Watershed 3 boundary
2. Breaches depressions in the DEM (fills flow-blocking sinks) using WBT
3. Computes FD8 flow accumulation (specific contributing area)
4. Computes slope in degrees
5. Calculates TWI as `log(SCA / tan(slope))` and masks it to the watershed

---

### Step 2 — `precompute_all.R`

**What it does:** Reads the TWI raster, clipped DEM, and watershed shapefile and converts them to plain R matrices and vectors saved as `.rds` files. The app reads these files at startup without needing any spatial packages.

**Run it** (after Step 1 is complete):

```r
source("precompute_all.R")
```

You will see:

```
Loading TWI raster ...
  TWI raster: 245 rows x 189 cols, extent [...] x [...]
Loading clipped DEM ...
  DEM raster: 245 rows x 189 cols, elevation range 488.2 – 774.6 m
Loading watershed shapefile ...
  Watershed boundary: 312 vertices
Computing TWI classes (TOPIDX) ...
  16 TWI classes | lambda (mean TWI) = 7.431
Saving .rds files to Data/ ...

════════════════════════════════════════════
  DONE! Created 7 .rds files in Data/
  Next step: deploy app.R to shinyapps.io
  The app loads these .rds files at startup
  and does NOT require raster, sf, or terra.
════════════════════════════════════════════
```

**What it produces** (7 files in `Data/`):

| File | Contents |
|---|---|
| `twi_matrix.rds` | TWI raster as a 2D numeric matrix |
| `twi_extent.rds` | Named vector `c(xmin, xmax, ymin, ymax)` |
| `dem_matrix.rds` | Clipped DEM as a 2D numeric matrix |
| `dem_extent.rds` | Named vector `c(xmin, xmax, ymin, ymax)` |
| `ws_boundary.rds` | Watershed polygon as an (X, Y) coordinate matrix |
| `topidx.rds` | 16-class TWI histogram table (midpoints + area fractions) |
| `twi_lam.rds` | Area-weighted mean TWI scalar (lambda, λ) |

---

## Running the App

Once all seven `.rds` files exist in `Data/`, launch the app:

```r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**.

---

## Using the App

### Sidebar Controls

**Catchment selector** — choose between the default HBEF Watershed 3 dataset or upload a custom catchment (see the **Custom Data Help** tab in the app for the required file formats).

**Viewing Period** — year range slider that controls the time window displayed on all plots. The model runs on the full selected period.

**▶ Run Model** — runs TOPMODEL with the current slider parameters and updates all plots.

**⬇ Export Output (CSV)** — downloads a CSV of all daily model outputs (discharge, ET, SWE, storage deficit, saturated fraction, and more).

**⚡ Optimize (Nelder-Mead)** — automatically calibrates parameters over the calibration window using the Nelder-Mead algorithm. Tries 3 starting points and keeps the best result. Updates the sliders and re-runs the model. Note: this is a local optimizer — results depend on starting values.

**Spin-Up & Calibration** — set the number of spin-up days excluded from NSE scoring, and the calibration period used by the optimizer.

### Parameter Groups

| Group | Parameters | What They Control |
|---|---|---|
| Subsurface Flow | `ln(Te)`, `m`, `td` | Baseflow magnitude, drainage speed, and time delay |
| Root Zone & ET | `Srmax` | Maximum soil moisture storage and evapotranspiration |
| Snow | `Tmelt`, `Melt Factor` | Snow/rain threshold and degree-day melt rate |

### NSE Score

The plot title reports NSE for both the **calibration period** and the current **viewing window** separately. NSE interpretation: ≥ 0.5 is generally considered acceptable, ≥ 0.7 is good, and 1.0 is a perfect fit.

### Plot Tabs

| Tab | Description |
|---|---|
| **Model Output** | 7-panel stacked time series: discharge, overland flow, baseflow, ET, SWE, storage deficit, and saturated area. X-axis is restricted to the viewing window. |
| **Input Data** | 3-panel plot of the raw forcing data: precipitation, air temperature, and observed discharge. |
| **Observed vs Simulated** | Scatter plot of observed vs. simulated discharge with a 1:1 reference line. |
| **Flow Duration** | Log-scale exceedance curve for observed and simulated flow. |
| **TWI Map** | Static map of the Topographic Wetness Index with watershed boundary and elevation contours. |
| **TWI Histogram** | Distribution of TWI values with adjustable class count and summary statistics. |
| **Saturation Map** | Animated daily map of saturated vs. unsaturated cells. Use the slider or click on the discharge mini-plot to jump to any date. |
| **Custom Data Help** | Instructions and R code for preparing custom catchment files. |

---

## Data Requirements

All input files should be placed in the `Data/` directory (note: capital D to match the repository structure):

| File | Description |
|---|---|
| `dem_10m.tif` | 10-meter resolution Digital Elevation Model for HBEF W3 |
| `Watershed3HB.shp` | Watershed 3 boundary polygon (+ `.dbf`, `.prj`, `.shx`) |
| `DailyWatershed.csv` | Daily precipitation — must include `DATE` and `watershed` columns; filtered to `W3` |
| `HBEF_DailyStreamflow_1956-2024.csv` | Daily streamflow — must include `DATE`, `WS`, and `Streamflow` columns; filtered to `WS == 3` |
| `HBEF_air_temp_daily.csv` | Daily air temperature — must include `date`, `STA`, and `AVE` columns; stations `STA14` and `STA17` are used |
| `twi_precomputed.tif` | **Produced by `precompute_twi.R`** |
| `twi_matrix.rds` | **Produced by `precompute_all.R`** |
| `twi_extent.rds` | **Produced by `precompute_all.R`** |
| `dem_matrix.rds` | **Produced by `precompute_all.R`** |
| `dem_extent.rds` | **Produced by `precompute_all.R`** |
| `ws_boundary.rds` | **Produced by `precompute_all.R`** |
| `topidx.rds` | **Produced by `precompute_all.R`** |
| `twi_lam.rds` | **Produced by `precompute_all.R`** |