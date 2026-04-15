---
editor_options: 
  markdown: 
    wrap: 72
---

# HBEF Watershed 3 — TOPMODEL Explorer

**EDS Capstone Project \| Virginia Tech** Samuel Handel · Kaelyn Harvey
· Max Hughes

------------------------------------------------------------------------

## Table of Contents

1.  [Introduction & Motivation](#1-introduction--motivation)
2.  [Application Overview](#2-application-overview)
3.  [Summary of Functions](#3-summary-of-functions)
4.  [How to Use the App](#4-how-to-use-the-app)
    -   [Setup & Precompute Workflow](#setup--precompute-workflow)
    -   [Running the App](#running-the-app)
    -   [Sidebar Controls](#sidebar-controls)
    -   [Plot Tabs](#plot-tabs)
    -   [Custom Catchment Upload](#custom-catchment-upload)
5.  [Libraries & Package Versions](#5-libraries--package-versions)
6.  [Data Used & Metadata](#6-data-used--metadata)
    -   [Data Sources](#data-sources)
    -   [Data Structure](#data-structure)
    -   [Preprocessing](#preprocessing)
7.  [Project Structure](#7-project-structure)
8.  [Known Issues](#8-known-issues)
9.  [Future Features](#9-future-features)
10. [References & Resources](#10-references--resources)

------------------------------------------------------------------------

## 1. Introduction & Motivation

Hydrological models are fundamental tools in environmental science,
enabling researchers and land managers to understand how precipitation
moves through a watershed — as streamflow, groundwater recharge,
evapotranspiration, and soil storage. However, most operational
hydrological models are implemented in static scripts or proprietary
software, making them inaccessible to non-specialists and difficult to
use interactively for exploration and education.

This project addresses that gap by wrapping **TOPMODEL** (Beven &
Kirkby, 1979) — one of the most widely taught and used conceptual
rainfall-runoff models — in an interactive Shiny web application. The
app is built around **Watershed 3 (W3) of the Hubbard Brook Experimental
Forest (HBEF)** in the White Mountains of New Hampshire, one of the
longest-running and most thoroughly documented small watershed research
sites in North America. HBEF has been continuously monitored since 1955,
providing over 60 years of daily precipitation, streamflow, and air
temperature records that make it ideal for hydrological model
calibration and validation.

**Why TOPMODEL?** TOPMODEL is a semi-distributed, physically-based model
that predicts streamflow by linking the watershed's topographic
structure — specifically the Topographic Wetness Index (TWI) — to
spatial patterns of soil moisture and saturation. Unlike lumped models,
TOPMODEL produces spatially explicit predictions of which areas of the
landscape are saturated on any given day, making it well-suited for
visualizing hydrological processes rather than just fitting discharge
curves.

**Why an interactive app?** Interactive parameter exploration allows
users to develop an intuitive understanding of how each model parameter
controls the shape of the hydrograph — something that is difficult to
convey through static figures or tables. The optimizer further
accelerates learning by automatically finding parameter combinations
that maximize model skill (Nash-Sutcliffe Efficiency), freeing users to
focus on interpreting results rather than manual trial-and-error
calibration.

The app supports both the default HBEF dataset and user-uploaded custom
catchments, making it a generalizable teaching and research tool beyond
the specific HBEF context.

------------------------------------------------------------------------

## 2. Application Overview

The app is built with **R Shiny** and consists of three scripts:

| Script | Purpose | When to run |
|------------------------|------------------------|------------------------|
| `precompute_twi.R` | Computes TWI from the raw DEM using WhiteboxTools | Once, before first use |
| `precompute_all.R` | Converts spatial outputs to `.rds` files for the app | Once, after `precompute_twi.R` |
| `app.R` | The main Shiny application | Every session |

The precompute scripts are run locally once and produce a set of plain
`.rds` files in `Data/`. The deployed app reads only these `.rds` files
at startup — it does not require any spatial packages (`terra`, `sf`,
`raster`) at runtime, which keeps the shinyapps.io deployment lean and
avoids binary dependency issues.

The model itself is implemented entirely in base R inside `app.R` with
no external TOPMODEL package dependency. This gives full transparency
over the equations and makes the code easier to modify and extend.

------------------------------------------------------------------------

## 3. Summary of Functions

All core functions are defined in `app.R`. The precompute scripts
contain standalone processing code without functions (except where
noted).

------------------------------------------------------------------------

### `load_hbef_data()`

**File:** `app.R` \| **Section 1**

Reads and merges the three HBEF input time series: precipitation,
streamflow, and air temperature. Filters precipitation to Watershed 3,
filters streamflow to `WS == 3`, restricts temperature to stations STA14
and STA17, and inner-joins all three on date. Missing precipitation
values are filled with zero; missing temperature values are linearly
interpolated between valid points.

**Returns:** A data frame with columns `date`, `precip_mm`, `flow_mm`,
`tavg`.

------------------------------------------------------------------------

### `hamon_pet(tavg, doy, lat_deg)`

**File:** `app.R` \| **Section 2**

Estimates daily potential evapotranspiration (PET) using the Hamon
(1961) equation. PET is computed from the daily mean air temperature and
the astronomical day length at the watershed latitude. Returns zero on
days with mean temperature at or below 0°C.

**Arguments:** - `tavg` — daily mean air temperature (°C) - `doy` — day
of year (integer 1–365/366) - `lat_deg` — watershed latitude in decimal
degrees

**Returns:** A numeric vector of daily PET values (mm/day).

**Equation:** \> PET = 0.1651 × Daylight × e_sat / (T + 273.3) × 29.8

where daylight is computed from the sunset hour angle at the given
latitude.

------------------------------------------------------------------------

### `run_topmodel(pars, data, topidx, lat_deg, spinup_days)`

**File:** `app.R` \| **Section 3**

Full daily TOPMODEL simulation. Implements the six core process
components in a daily time-step loop:

1.  **Snow** — degree-day accumulation and melt
2.  **Saturated area** — binary search on the sorted TWI threshold
    lookup table
3.  **Overland flow** — saturation-excess runoff on the wet fraction
4.  **Root zone** — water balance and actual ET via Hamon PET
5.  **Unsaturated zone** — time-delayed drainage to the water table
6.  **Baseflow** — exponential decay from saturated zone transmissivity

The initial saturation deficit is inferred from the first non-zero
observed flow record rather than being a free parameter, reducing the
number of parameters and the risk of a poor starting state.

**Arguments:** - `pars` — named numeric vector of model parameters (see
[Parameter Groups](#sidebar-controls)) - `data` — daily data frame from
`load_hbef_data()` or a custom upload - `topidx` — 2-column matrix of
TWI class midpoints and area fractions - `lat_deg` — watershed latitude
(decimal degrees) for Hamon PET - `spinup_days` — number of warm-up days
excluded from NSE scoring (default 365)

**Returns:** A named list containing daily time series for all model
outputs (`Qsim`, `Qobs`, `Qb`, `Qof`, `ET`, `SWE`, `Sd`, `SF`), model
dates, NSE score, warm-up length, and the mean TWI lambda.

------------------------------------------------------------------------

### `topmodel_nse(pv, P, tavg, obs, pet, topidx, n, warmup, lo, up)`

**File:** `app.R` \| **Section 3b**

Lightweight version of `run_topmodel()` used exclusively by the
optimizer. Implements identical TOPMODEL physics but skips storing daily
output arrays, accumulating only the NSE numerator and denominator
on-the-fly. Returns negative NSE (for minimization by `optim()`).
Returns `1e6` for out-of-bounds or non-finite parameter vectors.

------------------------------------------------------------------------

### `optimize_params(data, topidx, lat_deg, cal_start, cal_end, spinup_days, progress_fn)`

**File:** `app.R` \| **Section 4**

Runs Nelder-Mead optimization (`optim()`) from three different parameter
starting points over a user-defined calibration window. Keeps the result
with the highest NSE across all three starts. The calibration window
subsets the data before calling `topmodel_nse()`, which keeps run time
manageable. Parameter bounds are enforced both inside the objective
function and by clamping the final result.

**Optimized parameters:** `ln(Te)`, `m`, `Srmax`, `td`, `snow_t`,
`snow_mf`

**Starting points (three):**

| Start | ln(Te) | m   | Srmax | td  | Tmelt | Melt |
|-------|--------|-----|-------|-----|-------|------|
| 1     | 2.0    | 23  | 68    | 1   | −0.7  | 2.9  |
| 2     | 1.0    | 30  | 100   | 10  | 0.0   | 3.0  |
| 3     | −2.0   | 40  | 150   | 20  | −0.5  | 3.5  |

------------------------------------------------------------------------

### `compute_nse(res, start_date, end_date)`

**File:** `app.R` \| **Section 4 (helper)**

Computes Nash-Sutcliffe Efficiency over any arbitrary date window from
an existing model result. Used to display both calibration-period NSE
and viewing-window NSE in the plot title.

------------------------------------------------------------------------

### `snap(val, mn, mx, step)`

**File:** `app.R` \| **Section 4 (helper)**

Rounds a continuous value (e.g., from `optim()`) to the nearest slider
step increment so that `updateSliderInput()` places the slider on a
valid tick mark.

------------------------------------------------------------------------

### `make_topidx_classes(twi_raster, n_classes)` *(precompute only)*

**File:** `precompute_twi.R` \| **Section 4**

Converts a TWI raster into a TOPIDX lookup table by binning cell values
into `n_classes` histogram classes and computing the area fraction in
each class. Returns a 2-column matrix with columns `twi` (class
midpoints) and `frac` (area fractions summing to 1).

------------------------------------------------------------------------

## 4. How to Use the App

### Setup & Precompute Workflow

Both precompute scripts must be run **once** before launching the app.
They require the raw input files in `Data/` and produce the `.rds` files
the app depends on.

#### Prerequisites

Install all required R packages:

``` r
install.packages(c(
  # App dependencies
  "shiny", "bslib", "plotly", "dplyr", "lubridate",
  # Precompute dependencies (not needed by the deployed app)
  "terra", "sf", "whitebox", "raster"
))
```

Initialize WhiteboxTools (one-time binary download):

``` r
whitebox::install_whitebox()
```

#### Step 1 — Run `precompute_twi.R`

Computes TWI from the raw 10m DEM using WhiteboxTools.

``` r
source("precompute_twi.R")
```

Expected output:

```         
Loading DEM and watershed ...
Breaching depressions ...
FD8 flow accumulation ...
Computing slope ...
Computing TWI ...
Done! Saved to: Data/twi_precomputed.tif
```

#### Step 2 — Run `precompute_all.R`

Converts the TWI raster, clipped DEM, and watershed boundary into `.rds`
files.

``` r
source("precompute_all.R")
```

Expected output:

```         
Loading TWI raster ...
Loading and clipping DEM to watershed ...
Loading watershed boundary ...
Computing TWI classes (TOPIDX) ...
  16 TWI classes | lambda (mean TWI) = 7.431
Saving .rds files to Data/ ...
Done! Saved to: Data/
========================================
  DONE! Created 7 .rds files in Data/
========================================
```

------------------------------------------------------------------------

### Running the App {#running-the-app}

``` r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**. The app auto-runs the
model on startup so all plots are populated immediately.

------------------------------------------------------------------------

### Sidebar Controls {#sidebar-controls}

#### Catchment

Select **HBEF Watershed 3 (default)** to use the pre-loaded HBEF data,
or **Custom upload** to run the model on a different catchment (see
[Custom Catchment Upload](#custom-catchment-upload)).

#### Viewing Period

Year range slider that controls the time window displayed on all plots.
The x-axis of every plot is clipped to this window; the model itself
runs on the full selected period.

#### Run Model

Runs TOPMODEL with the current slider values and refreshes all plots.
The NSE score for both the calibration period and viewing window is
displayed in the Model Output plot title.

#### Export Output (CSV)

Downloads a daily CSV of all model outputs. Columns: `date`, `P_mm`,
`T_C`, `Qobs_mm`, `Qsim_mm`, `Qbaseflow_mm`, `Qsat_excess_mm`, `ET_mm`,
`SWE_mm`, `Storage_Deficit_mm`, `Saturated_Fraction`.

#### Optimize (Nelder-Mead)

Calibrates model parameters over the selected calibration window. Tries
three starting points; keeps the best result. Updates the sliders with
the optimized values and re-runs the model automatically.

> **Note:** Nelder-Mead is a local optimizer — results can vary
> depending on starting values and the length of the calibration period.
> Longer calibration periods (3–5 years) generally produce more stable
> results.

#### Parameter Groups

| Group | Parameter | Full Name | Units | Default | Range | Effect |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| Subsurface Flow | `ln(Te)` | Log Soil Transmissivity | — | 1.9 | −7 to 5 | Controls baseflow magnitude; higher = more baseflow |
| Subsurface Flow | `m` | Transmissivity Decay with Depth | mm | 19 | 5–100 | Small = flashy response; large = slow recession |
| Subsurface Flow | `td` | Unsaturated Zone Time Delay | days | 3 | 1–60 | Delays drainage to the water table |
| Root Zone & ET | `Srmax` | Maximum Root Zone Storage | mm | 65 | 5–300 | Controls ET; larger = more summer ET |
| Snow | `Tmelt` | Snow/Rain Temperature Threshold | °C | −0.7 | −2 to 2 | Below this, precipitation falls as snow |
| Snow | `Melt Factor` | Degree-Day Melt Coefficient | mm/°C/day | 3.4 | 1–6 | Snowmelt per degree above Tmelt per day |

#### Spin-Up & Calibration

-   **Spin-up (days):** Days excluded from NSE computation while the
    model state stabilizes. Default 365 days. A one-year warm-up is
    recommended for most applications.
-   **Calibration Years:** The date window used by the Optimize button.
    Should overlap with a period of good data quality and represent the
    full range of hydrological conditions.

------------------------------------------------------------------------

### Plot Tabs {#plot-tabs}

#### Model Output

Seven stacked time series panels sharing an x-axis, restricted to the
viewing window:

| Panel | Variable                         | Color      |
|-------|----------------------------------|------------|
| 1     | Observed and simulated discharge | Blue / Red |
| 2     | Saturation-excess overland flow  | Orange     |
| 3     | Baseflow (subsurface)            | Purple     |
| 4     | Actual evapotranspiration        | Green      |
| 5     | Snow water equivalent            | Blue       |
| 6     | Watershed storage deficit        | Red        |
| 7     | Saturated catchment fraction     | Dark blue  |

The plot title reports NSE for both the calibration period and the
viewing window.

#### Input Data

Three panels showing the raw forcing data (precipitation, air
temperature, observed discharge) for the viewing window — useful for
diagnosing model performance relative to individual storm events or
seasons.

#### Observed vs Simulated

Scatter plot of daily observed vs. simulated discharge for all
post-warm-up days in the viewing window. A dashed 1:1 line indicates a
perfect model. Points clustered close to this line indicate good model
skill across all flow magnitudes.

#### Flow Duration Curve

Log-scale exceedance probability plot. Shows how often each discharge
magnitude is exceeded. A good model fit across both the high-exceedance
(low flows) and low-exceedance (high flows) portions of the curve
indicates skill at both baseflow recession and peak storm response.

#### TWI Map

Static map of the Topographic Wetness Index raster with watershed
boundary and DEM elevation contours overlaid. The color scale (yellow →
blue) shows where the landscape is topographically predisposed to
accumulate water (high TWI = valley bottoms and concave slopes).

#### TWI Histogram

Distribution of TWI cell values across the watershed. The adjustable
class slider (4–40 classes) allows inspection of the TWI distribution
shape. Summary statistics (min, max, mean, median, SD, N, lambda) are
displayed alongside the histogram.

#### Saturation Map

Animated daily binary map showing saturated (blue) vs. unsaturated (tan)
cells. For each day, cells with TWI ≥ (λ + Sd/m) are predicted to have
their water table at or above the surface. Use the **Day slider** to
step through time, or click **▶** to animate. **Click any point on the
discharge mini-plot** to jump the map directly to that date.

------------------------------------------------------------------------

### Custom Catchment Upload {#custom-catchment-upload}

To run the model on a catchment other than HBEF W3, select **Custom
upload** in the sidebar and provide:

**Required:** - **Time series CSV** — columns: `date` (YYYY-MM-DD),
`precip_mm`, `tavg`, `flow_mm` - **TOPIDX CSV** — columns: `twi` (class
midpoints), `area_frac` (area fractions summing to 1) - **Latitude** —
decimal degrees

**Optional (for spatial maps):** - **TWI matrix (.rds)** — 2D numeric
matrix with `attr(., "extent") = c(xmin, xmax, ymin, ymax)`; max 250,000
cells - **DEM matrix (.rds)** — same format, elevation values in
metres - **Watershed boundary (.rds)** — 2-column matrix of polygon (x,
y) coordinates

R code to prepare these files from a local raster:

``` r
library(terra)
twi <- rast("twi.tif")
twi_mat <- as.matrix(twi, wide = TRUE)
attr(twi_mat, "extent") <- c(xmin(twi), xmax(twi), ymin(twi), ymax(twi))
saveRDS(twi_mat, "twi_matrix.rds")

vals <- na.omit(values(twi))
h <- hist(vals, breaks = 30, plot = FALSE)
topidx <- data.frame(twi = h$mids, area_frac = h$counts / sum(h$counts))
write.csv(topidx, "topidx.csv", row.names = FALSE)
```

------------------------------------------------------------------------

## 5. Libraries & Package Versions

### App Dependencies (`app.R`)

| Package     | Version (tested) | Purpose                                  |
|-------------|------------------|------------------------------------------|
| `shiny`     | ≥ 1.8.0          | Web application framework                |
| `bslib`     | ≥ 0.6.0          | Bootstrap 5 UI themes (`bs_theme`)       |
| `plotly`    | ≥ 4.10.0         | Interactive plots and subplots           |
| `dplyr`     | ≥ 1.1.0          | Data wrangling (`group_by`, `summarise`) |
| `lubridate` | ≥ 1.9.0          | Date handling (`yday`)                   |

### Precompute Dependencies (not required at runtime)

| Package | Version (tested) | Purpose |
|------------------------|------------------------|------------------------|
| `terra` | ≥ 1.7.0 | Raster/vector operations (`rast`, `mask`, `crop`, `ext`) |
| `sf` | ≥ 1.0.0 | Shapefile reading and CRS transformation |
| `whitebox` | ≥ 2.3.0 | WhiteboxTools interface for DEM processing |
| `raster` | ≥ 3.6.0 | Legacy raster support (used only in `precompute_all.R`) |

To check your installed versions:

``` r
packageVersion("shiny")
packageVersion("terra")
# etc.
```

> **Note on spatial packages:** `terra`, `sf`, and `raster` are **not**
> listed in the app's `DESCRIPTION` or `rsconnect` dependencies. They
> are only needed to run the two precompute scripts locally. The
> deployed app reads only `.rds` files and requires none of these
> packages.

------------------------------------------------------------------------

## 6. Data Used & Metadata

### Data Sources {#data-sources}

All HBEF data are publicly available from the Hubbard Brook Ecosystem
Study data repository and the Environmental Data Initiative (EDI).

| Dataset | File | Source | Coverage |
|------------------|------------------|------------------|------------------|
| Daily watershed precipitation | `DailyWatershed.csv` | HBEF Data Archive | 1956–2024 |
| Daily streamflow | `HBEF_DailyStreamflow_1956-2024.csv` | HBEF Data Archive / EDI | 1956–2024 |
| Daily air temperature | `HBEF_air_temp_daily.csv` | HBEF Data Archive | \~1956–2024 |
| Digital Elevation Model | `dem_10m.tif` | USGS National Elevation Dataset (NED) | HBEF W3 extent |
| Watershed boundary | `Watershed3HB.shp` | HBEF GIS Data Archive | W3 polygon |

**Hubbard Brook Experimental Forest** is located in the White Mountain
National Forest, New Hampshire (43.94°N, 71.70°W). Watershed 3 is a
42.4-hectare forested headwater catchment with a mean elevation of
approximately 620 m and a mean annual precipitation of \~1400 mm
(roughly 30% as snow).

------------------------------------------------------------------------

### Data Structure {#data-structure}

#### `DailyWatershed.csv`

Daily precipitation totals for all HBEF watersheds.

| Column      | Type      | Units      | Description               |
|-------------|-----------|------------|---------------------------|
| `DATE`      | Date      | YYYY-MM-DD | Measurement date          |
| `watershed` | Character | —          | Watershed ID (e.g., `W3`) |
| `Precip`    | Numeric   | mm         | Total daily precipitation |

The app filters to `watershed == "W3"`.

#### `HBEF_DailyStreamflow_1956-2024.csv`

Daily streamflow for all HBEF watersheds.

| Column | Type | Units | Description |
|------------------|------------------|------------------|------------------|
| `DATE` | Date | YYYY-MM-DD | Measurement date |
| `WS` | Integer | — | Watershed number (app uses `WS == 3`) |
| `Streamflow` | Numeric | mm | Daily mean streamflow, expressed as mm over watershed area |

#### `HBEF_air_temp_daily.csv`

Daily air temperature from multiple HBEF meteorological stations.

| Column | Type      | Units      | Description               |
|--------|-----------|------------|---------------------------|
| `date` | Date      | YYYY-MM-DD | Measurement date          |
| `STA`  | Character | —          | Station ID                |
| `AVE`  | Numeric   | °C         | Daily average temperature |

The app uses stations `STA14` and `STA17` only and averages across them
for each day.

#### `dem_10m.tif`

10-meter resolution Digital Elevation Model in GeoTIFF format, projected
in a UTM coordinate reference system. Covers the full extent of HBEF
Watershed 3 with a boundary buffer. Elevation values are in metres above
sea level.

#### `Watershed3HB.shp`

Polygon shapefile of the Watershed 3 boundary. Used to clip the DEM and
TWI rasters and to draw the watershed outline on spatial maps.

------------------------------------------------------------------------

### Preprocessing {#preprocessing}

#### Precipitation

No gap-filling is applied. Remaining `NA` values are replaced with `0`
(conservative assumption: no data = no precipitation).

#### Streamflow

Missing daily values (`NA`) are retained and excluded from NSE
computation. The model is run on all days regardless of whether observed
flow is available.

#### Air Temperature

Missing values are linearly interpolated between valid observations
using `approx()` with `rule = 2` (flat extrapolation at the edges). This
is necessary because Hamon PET requires a complete temperature record.

#### TWI Computation (`precompute_twi.R`)

1.  The raw DEM is clipped and masked to the Watershed 3 boundary
2.  Depressions in the DEM are breached using
    `wbt_breach_depressions_least_cost()` with a maximum breach distance
    of 10 cells — this removes artificial sinks while preserving true
    depressions
3.  FD8 (D-infinity multi-flow) accumulation is computed as specific
    contributing area (m²/m)
4.  Slope is computed in degrees
5.  TWI = ln(SCA / tan(β)), where SCA is specific contributing area and
    β is slope. A minimum tan(β) of 0.001 is enforced to prevent
    infinite TWI values on flat cells
6.  The TWI raster is masked back to the watershed boundary

#### TOPIDX Construction (`precompute_all.R`)

The TWI raster is binned into 16 histogram classes. The area fraction in
each class is computed as `cell count / total non-NA cells`. The
16-class resolution matches the classic TOPMODEL parameterization (Beven
& Kirkby, 1979) and provides sufficient resolution for the relatively
small HBEF W3 watershed.

------------------------------------------------------------------------

## 7. Project Structure

```         
project/
├── app.R                    # Main Shiny application (model + UI + server)
├── precompute_twi.R         # Step 1: compute TWI via WhiteboxTools (run once)
├── precompute_all.R         # Step 2: convert spatial data to .rds (run once)
├── README.md                # This file
└── Data/
    │
    ├── ── Raw inputs (checked into repository) ──────────────────────────────
    ├── dem_10m.tif                          # 10m DEM (USGS NED)
    ├── Watershed3HB.shp                     # Watershed boundary
    ├── Watershed3HB.dbf                     # Shapefile attribute table
    ├── Watershed3HB.prj                     # Shapefile projection
    ├── Watershed3HB.shx                     # Shapefile index
    ├── DailyWatershed.csv                   # Daily precipitation (all watersheds)
    ├── HBEF_DailyStreamflow_1956-2024.csv   # Daily streamflow (all watersheds)
    ├── HBEF_air_temp_daily.csv              # Daily air temperature
    │
    └── ── Generated by precompute scripts (not committed to repository) ─────
        ├── twi_precomputed.tif              # TWI raster (from precompute_twi.R)
        ├── twi_matrix.rds                   # TWI as 2D matrix
        ├── twi_extent.rds                   # TWI spatial extent vector
        ├── dem_matrix.rds                   # DEM as 2D matrix (clipped)
        ├── dem_extent.rds                   # DEM spatial extent vector
        ├── ws_boundary.rds                  # Watershed polygon vertices
        ├── topidx.rds                       # TOPIDX class table (16 x 2)
        └── twi_lam.rds                      # Mean TWI scalar (lambda)
```

> **Git note:** The generated `.rds` and `.tif` files should be added to
> `.gitignore`. Collaborators clone the repo and run the two precompute
> scripts to regenerate them locally before launching the app.

------------------------------------------------------------------------

## 8. Known Issues

**Optimizer sensitivity to calibration period length.** The Nelder-Mead
optimizer is a local search method and can converge to different local
optima depending on the selected calibration window. Short windows (\< 2
years) are particularly prone to overfitting. Using a 3–5 year
calibration period that spans wet and dry years is recommended. The app
tries three starting points as a partial mitigation, but a global search
(e.g., Monte Carlo sampling) would be more robust.

**Snow model limitations.** The degree-day snow model is a
simplification. It uses a single melt factor and temperature threshold
applied uniformly across the watershed, ignoring elevation gradients and
aspect effects. In a watershed with significant topographic relief,
actual snowmelt timing can differ substantially from the model's
predictions, particularly in spring.

**Single temperature station.** The app averages across two HBEF
meteorological stations (STA14, STA17). Spatial variability in
temperature across the watershed — important for snow accumulation and
melt — is not captured.

**No routing.** TOPMODEL assumes instantaneous routing of generated
runoff to the watershed outlet. In larger catchments or those with long
channel networks, travel time effects could introduce lag between
modeled and observed peak flows that TOPMODEL cannot reproduce.

**ET is limited by Hamon PET.** The Hamon equation is a simple
temperature-based PET method. It does not account for vapor pressure
deficit, wind speed, or net radiation, which are the primary drivers of
actual evapotranspiration at daily timescales. This limits model
accuracy during periods where energy supply (not water supply) limits
ET.

**Custom catchments: no spatial maps if extent attribute is missing.**
If a user uploads a TWI matrix `.rds` file without an
`attr(., "extent")` attribute, the app skips spatial maps and shows a
placeholder message. The time series model still runs correctly.

**ShinyApps.io memory limits.** Large custom catchments (\> 250,000 TWI
cells) are blocked by the app to prevent memory issues on the free
shinyapps.io tier. This limit can be adjusted in `app.R` by changing the
`CELL_LIMIT` constant.

------------------------------------------------------------------------

## 9. Future Features

### Additional Catchments

-   **Multi-watershed support:** Add a dropdown to switch between all
    nine HBEF watersheds (W1–W9) using pre-computed `.rds` files for
    each, enabling paired watershed comparisons (e.g., logged vs.
    reference).
-   **Larger catchment support:** Implement catchment aggregation or
    raster resampling so that larger DEM extents can be used without
    exceeding memory limits.
-   **Automated TOPIDX generation:** Add a built-in workflow inside the
    app (or a companion script) that computes TWI and TOPIDX from a
    user-uploaded raw DEM, removing the need for users to run
    WhiteboxTools separately.

### Model Parameters & Physics

-   **Calibrate `qs0` (initial subsurface flow):** The current version
    derives the initial saturation deficit from the first observed flow
    value. Exposing `qs0` as an optimizable parameter would give more
    control over the warm-up state.
-   **Spatial ET heterogeneity:** Allow `Srmax` to vary by land cover
    class (e.g., using a simple forest/non-forest mask), which would
    improve ET simulation in partially disturbed catchments.
-   **Improved snow model:** Implement a simple energy-balance snow
    model or an elevation-banded degree-day model to better capture the
    spatial pattern of snowmelt across watersheds with significant
    topographic relief.
-   **Uncertainty quantification:** Add Monte Carlo or GLUE (Generalized
    Likelihood Uncertainty Estimation) sampling to produce ensemble
    predictions and visualize parameter uncertainty rather than a single
    best-fit solution.
-   **Multi-objective calibration:** Currently NSE is the only
    calibration metric. Adding log-NSE (which weights low flows more
    heavily) or a composite objective function would allow users to
    calibrate for different flow regimes.

------------------------------------------------------------------------

## 10. References & Resources

### TOPMODEL

Beven, K.J., & Kirkby, M.J. (1979). A physically based, variable
contributing area model of basin hydrology. *Hydrological Sciences
Bulletin*, 24(1), 43–69. <https://doi.org/10.1080/02626667909491834>

Beven, K.J. (2012). *Rainfall-Runoff Modelling: The Primer* (2nd ed.).
Wiley-Blackwell. — Comprehensive textbook covering TOPMODEL theory and
implementation.

Beven, K.J., Lamb, R., Quinn, P., Romanowicz, R., & Freer, J. (1995).
TOPMODEL. In V.P. Singh (Ed.), *Computer Models of Watershed Hydrology*
(pp. 627–668). Water Resources Publications. — Original TOPMODEL
implementation reference.

### Hubbard Brook Experimental Forest

Likens, G.E. (Ed.) (2013). *Biogeochemistry of a Forested Ecosystem*
(3rd ed.). Springer. — Foundational reference for HBEF ecosystem
science.

Hubbard Brook Ecosystem Study. (2024). *Hubbard Brook Experimental
Forest: Long-Term Ecological Research*. <https://hubbardbrook.org> —
Data portal and site description.

Campbell, J.L., Driscoll, C.T., Pourmokhtarian, A., & Hayhoe, K. (2011).
Streamflow responses to past and projected future changes in climate at
the Hubbard Brook Experimental Forest, New Hampshire, United States.
*Water Resources Research*, 47, W02514.

### Topographic Wetness Index & WhiteboxTools

Quinn, P., Beven, K., Chevallier, P., & Planchon, O. (1991). The
prediction of hillslope flow paths for distributed hydrological
modelling using digital terrain models. *Hydrological Processes*, 5(1),
59–79.

Lindsay, J.B. (2016). Whitebox GAT: A case study in geomorphometric
analysis. *Computers & Geosciences*, 95, 75–84.
<https://doi.org/10.1016/j.cageo.2016.07.003>

Wu, Q. (2020). *whitebox: WhiteboxTools R Frontend*. R package version
2.x. <https://cran.r-project.org/package=whitebox>

### Evapotranspiration

Hamon, W.R. (1961). Estimating potential evapotranspiration. *Journal of
the Hydraulics Division*, 87(3), 107–120.

### Model Evaluation

Nash, J.E., & Sutcliffe, J.V. (1970). River flow forecasting through
conceptual models part I — A discussion of principles. *Journal of
Hydrology*, 10(3), 282–290.
[https://doi.org/10.1016/0022-1694(70)90255-6](https://doi.org/10.1016/0022-1694(70)90255-6){.uri}
— Original paper defining the Nash-Sutcliffe Efficiency (NSE) criterion
used throughout this app.

Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel,
R.D., & Veith, T.L. (2007). Model evaluation guidelines for systematic
quantification of accuracy in watershed simulations. *Transactions of
the ASABE*, 50(3), 885–900. — Standard reference for interpreting NSE
thresholds (≥ 0.5 satisfactory, ≥ 0.65 good, ≥ 0.75 very good).

### R Packages & Tools

Chang, W., Cheng, J., Allaire, J.J., Sievert, C., Schloerke, B., Xie,
Y., Allen, J., McPherson, J., Dipert, A., & Borges, B. (2024). *shiny:
Web Application Framework for R*. R package version 1.8.x.
<https://cran.r-project.org/package=shiny>

Sievert, C. (2020). *Interactive Web-Based Data Visualization with R,
plotly, and shiny*. Chapman and Hall/CRC. <https://plotly-r.com>

Hijmans, R.J. (2024). *terra: Spatial Data Analysis*. R package version
1.7.x. <https://cran.r-project.org/package=terra>

Pebesma, E. (2018). Simple features for R: Standardized support for
spatial vector data. *The R Journal*, 10(1), 439–446.
<https://doi.org/10.32614/RJ-2018-009>

### Additional Reading

Beven, K.J. (2001). How far can we go in distributed hydrological
modelling? *Hydrology and Earth System Sciences*, 5(1), 1–12. —
Discusses the limits of spatial complexity in watershed models.

Seibert, J., & McGlynn, B.L. (2007). A new triangular multiple flow
direction algorithm for computing upslope areas from gridded digital
elevation models. *Water Resources Research*, 43, W04501. — Background
on multi-directional flow accumulation methods (FD8) used in TWI
computation.
