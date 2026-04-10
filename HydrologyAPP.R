# ═══════════════════════════════════════════════════════════════════════════════
# app.R — HBEF Watershed 3 TOPMODEL Explorer
#
# A Shiny app for running TOPMODEL, a rainfall-runoff model that predicts
# streamflow based on watershed topography, soil properties, and climate inputs.
# Built for the Hubbard Brook Experimental Forest (HBEF) Watershed 3.
#
# Supports both the default HBEF dataset and user-uploaded custom catchments.
#
# Samuel Handel, Kaelyn Harvey, Max Hughes
# ═══════════════════════════════════════════════════════════════════════════════

# ── Core libraries ─────────────────────────────────────────────────────────────
library(shiny)      # web app framework
library(bslib)      # bootstrap-based UI themes
library(plotly)     # interactive plots
library(dplyr)      # data wrangling
library(lubridate)  # data handling

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  FILE PATHS
# ═══════════════════════════════════════════════════════════════════════════════
DATA_DIR <- "Data"  # all input data files go here

# Maximum number of raster cells allowed for custom spatial uploads
CELL_LIMIT <- 250000  # 500 x 500 cells

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD HBEF DATA
#
# Reads in & merges the 3 HBEF input time series:
#   - Precipitation:   DailyWatershed.csv
#   - Streamflow:      HBEF_DailyStreamflow_1956-2024.csv
#   - Air Temperature: HBEF_air_temp_daily.csv
# ═══════════════════════════════════════════════════════════════════════════════
load_hbef_data <- function() {
  # ── Precipitation ────────────────────────────────────────────────────────────
  precip <- read.csv(file.path(DATA_DIR, "DailyWatershed.csv"), stringsAsFactors = FALSE)
  precip$DATE <- as.Date(trimws(precip$DATE))
  precip <- precip[trimws(precip$watershed) == "W3", c("DATE", "Precip")]
  names(precip) <- c("date", "precip_mm")
  
  # ── Streamflow ───────────────────────────────────────────────────────────────
  flow <- read.csv(file.path(DATA_DIR, "HBEF_DailyStreamflow_1956-2024.csv"), stringsAsFactors = FALSE)
  flow$DATE <- as.Date(trimws(flow$DATE))
  flow <- flow[flow$WS == 3, c("DATE", "Streamflow")]
  names(flow) <- c("date", "flow_mm")
  
  # ── Air temperature ───────────────────────────────────────────────────────────
  temp <- read.csv(file.path(DATA_DIR, "HBEF_air_temp_daily.csv"), stringsAsFactors = FALSE)
  temp$date <- as.Date(trimws(temp$date))
  
  # Restrict to the 2 primary HBEF stations
  temp <- temp[temp$STA %in% c("STA14", "STA17"), ]
  temp_avg <- temp |> 
    group_by(date) |> 
    summarize(tavg = mean(AVE, na.rm = TRUE), .groups = "drop")
  
  # Remove any days where the temperature average is non-finite
  temp_avg <- temp_avg[!is.na(temp_avg$tavg) & is.finite(temp_avg$tavg), ]
  
  # ── Merge & clean ─────────────────────────────────────────────────────────────
  # Inner-join where it keeps only dates that appear in all 3 datasets
  df <- merge(precip, flow, by = "date", all = FALSE)
  df <- merge(df, temp_avg, by = "date", all = FALSE)
  df <- df[order(df$date), ]
  
  # Replace any remaining NA precipitation with 0
  df$precip_mm[is.na(df$precip_mm)] <- 0
  df
  }

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  HAMON PET (mm/day)
#
# Estimates potential evapotranspiration (PET) from temperature & day length.
# Returns zero on days where tavg ≤ 0°C.
# ═══════════════════════════════════════════════════════════════════════════════
hamon_pet <- function(tavg, doy, lat_deg) {
  lat_rad  <- lat_deg * pi / 180
  dec      <- 0.4093 * sin(2 * pi * (284 + doy) / 365)           # solar declination (radians)
  ws       <- acos(pmax(pmin(-tan(lat_rad) * tan(dec), 1), -1))  # sunset hour angle
  daylight <- 24 * ws / pi                                       # day length (hours)
  esat     <- 0.611 * exp(17.27 * tavg / (tavg + 237.3))        # saturation vapor pressure (kPa)
  pet      <- ifelse(tavg > 0,
                     0.1651 * daylight * esat / (tavg + 273.3) * 29.8, 0)
  pmax(pet, 0)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  TOPMODEL — full simulation with all outputs.
#
# Implements the TOPMODEL formulation:
#    - Snow accumulation & degree-day melt
#    - Saturated area fraction via TWI-based threshold (binary search)
#    - Saturation-excess overland flow on the wet fraction
#    - Root zone water balance and Hamon actual ET
#    - Unsaturated zone drainage delayed by td
#    - Exponential baseflow decay from saturated zone transmissivity
#
# The initial saturation deficit (Sd) is inferred from the first non-zero
# observed flow record, rather than being a free parameter, reducing the risk
# of a bad starting state.
#
# Returns a named list with daily time series for all fluxes and states, plus
# the full-record NSE and spin-up length (warmup).
# ═══════════════════════════════════════════════════════════════════════════════
run_topmodel <- function(pars, data, topidx, lat_deg, spinup_days = 365) {
  n    <- nrow(data)
  P    <- data$precip_mm
  tavg <- data$tavg
  obs  <- data$flow_mm
  
  # Pre-compute daily PET for the entire record
  pet  <- hamon_pet(tavg, yday(data$date), lat_deg)
  
  # ── Unpack parameters ────────────────────────────────────────────────────────
  lnTe    <- pars["lnTe"]          # log transmissivity at saturation
  m       <- pars["m"]             # transmissivity decay with depth (mm)
  Srmax   <- pars["Srmax"]         # maximum root zone storage capacity (mm)
  td      <- max(pars["td"], 0.1)  # unsaturated zone time delay (days)
  snow_t  <- pars["snow_t"]        # snow/rain temperature threshold (°C)
  snow_mf <- pars["snow_mf"]       # degree-day melt factor (mm/°C/day)
  
  # ── Build sorted TWI look-up table ────────────────────────────────────────────
  tv  <- topidx[, 1]; tf  <- topidx[, 2]
  lam <- sum(tv * tf)                   # area-weighted mean TWI (lambda)
  tt  <- m * (tv - lam)                 # local saturation thresholds
  si  <- order(tt, decreasing = TRUE)
  tts <- tt[si]; tfc <- cumsum(tf[si])  # sorted thresholds & cumulative fractions
  nc  <- length(tts)
  
  # ── Smart initial saturation deficit from first observed flow ─────────────────
  # Using observed flow avoids a hard-coded Sd0 parameter. If no valid flow
  # record exists, fall back to q0 = 1 mm/d.
  q0 <- obs[which(!is.na(obs) & obs > 0)[1]]
  if (is.na(q0) || !is.finite(q0)) q0 <- 1
  Sd <- -m * (log(max(q0, 0.01)) - lnTe)
  if (Sd > 2000) Sd <- 2000; if (Sd < 0) Sd <- 0
  
  # ── Initialize state variables ────────────────────────────────────────────────
  Srz  <- Srmax / 2  # root zone storage deficit
  Suz  <- 0          # unsaturated zone storage (mm)
  SWE  <- 0          # snow water equivalent (mm)
  Te   <- exp(lnTe)  # transmissivity at saturation
  td_m <- td / m     # scaled time delay
  
  # ── Pre-allocate output arrays ────────────────────────────────────────────────
  Qsim    <- numeric(n); Qb_out  <- numeric(n); Qof_out <- numeric(n)
  SWE_out <- numeric(n); ET_out  <- numeric(n); Sd_out  <- numeric(n)
  SF_out  <- numeric(n)
  
  # ════════════════════════════════════════════════════════════════════════════
  # Daily time-step loop
  # ════════════════════════════════════════════════════════════════════════════
  for (t in seq_len(n)) {
    tv2 <- tavg[t]
    
    # ── Snow accumulation & melt (degree-day model) ───────────────────────────
    if (tv2 <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0        # all precip stored as snow
    } else {
      ma  <- snow_mf * (tv2 - snow_t)     # potential snowmelt
      if (ma > SWE) ma <- SWE             # can't melt more than exists
      SWE <- SWE - ma; rain <- P[t] + ma  # liquid water = rain + melt
    }
    
    SWE_out[t] <- SWE
    
    # ── Saturated area fraction (binary search on pre-sorted thresholds) ──────
    # sf = fraction of watershed where TWI ≥ (lambda + Sd/m)
    if (Sd >= tts[1L]) { sf <- 0} 
    else if (Sd <  tts[nc]) { sf <- tfc[nc]} 
    else {
      lo2 <- 1L; hi <- nc
      while (hi - lo2 > 1L) {
        mid <- (lo2 + hi) %/% 2L
        if (tts[mid] > Sd) lo2 <- mid else hi <- mid
        }
      sf <- tfc[lo2]
      }
    
    SF_out[t] <- sf
    
    # ── Saturation-excess overland flow ───────────────────────────────────────
    qof   <- sf * rain   # rain falling on saturated cells runs off immediately
    infil <- rain - qof  # remainder infiltrates
    Qof_out[t] <- qof
    
    # ── Root zone water balance ───────────────────────────────────────────────
    Srz <- Srz - infil                           # infiltration reduces the root zone deficit
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }  # excess drains to unsaturated zone
    
    # Actual ET is proportional to how full the root zone is
    if (Srz < Srmax) { ae <- pet[t] * (1 - Srz / Srmax) } else { ae <- 0 }
    Srz <- Srz + ae; if (Srz > Srmax) Srz <- Srmax
    ET_out[t] <- ae
    
    # ── Unsaturated zone drainage (time-delayed by td) ────────────────────────
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd)
      if (quz > Suz) quz <- Suz
      Suz <- Suz - quz
      } else { quz <- 0 }
    
    # ── Saturated zone: update deficit then compute exponential baseflow ───────
    Sd <- Sd - quz; if (Sd < 0) Sd <- 0
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + 50) qb <- Sd + 50; if (qb < 0) qb <- 0
    Sd <- Sd + qb; if (Sd < 0) Sd <- 0  # outflow increases the deficit
    
    Qb_out[t]  <- qb
    Qsim[t]    <- qb + qof  # total discharge = baseflow + overland flow
    Sd_out[t]  <- Sd
  }
  
  # Remove any non-finite simulated values (numerical edge cases)
  Qsim[!is.finite(Qsim)] <- 0
  Qsim <- pmax(Qsim, 0)
  
  # ── Nash-Sutcliffe Efficiency (full record, post spin-up) ─────────────────
  warmup <- min(spinup_days, n - 1)
  idx    <- (warmup + 1):n
  o <- obs[idx]; s <- Qsim[idx]
  ok <- !is.na(o) & is.finite(s)
  if (sum(ok) < 30) {
    nse <- NA_real_
    } 
  else {
    ss_r <- sum((o[ok] - s[ok])^2)
    ss_t <- sum((o[ok] - mean(o[ok]))^2)
    nse  <- ifelse(ss_t > 0, 1 - ss_r / ss_t, NA_real_)
  }
  
  list(Qsim = Qsim, Qobs = obs, Qb = Qb_out, Qof = Qof_out,
       SWE = SWE_out, ET = ET_out, Sd = Sd_out, SF = SF_out,
       dates = data$date, nse = nse, warmup = warmup, m = m, lam = lam)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3b.  NSE-ONLY TOPMODEL (optimizer objective function)
#
# Returns negative NSE (optim() minimizes, so −NSE → maximizes NSE).
# Accepts pre-computed lower/upper bounds and returns 1e6 for invalid inputs.
# ═══════════════════════════════════════════════════════════════════════════════
topmodel_nse <- function(pv, P, tavg, obs, pet, topidx, n, warmup, lo, up) {
  
  # Rejects out-of-bounds or non-finite parameter vectors immediately
  if (any(!is.finite(pv)) || any(pv < lo) || any(pv > up)) return(1e6)
  
  lnTe <- pv[1]; m <- pv[2]; Srmax <- pv[3]
  td   <- pv[4]; snow_t <- pv[5]; snow_mf <- pv[6]
  if (td < 0.1) td <- 0.1
  
  # Build sorted TWI lookup table
  tv  <- topidx[, 1];    tf <- topidx[, 2];     lam <- sum(tv * tf)
  tt  <- m * (tv - lam); si <- order(tt, decreasing = TRUE)
  tts <- tt[si];         tfc <- cumsum(tf[si]); nc <- length(tts)
  
  # Smart Sd0 from first valid observed flow
  q0 <- obs[which(!is.na(obs) & obs > 0)[1]]
  if (is.na(q0) || !is.finite(q0)) q0 <- 1
  Sd <- -m * (log(max(q0, 0.01)) - lnTe)
  if (Sd > 2000) Sd <- 2000; if (Sd < 0) Sd <- 0
  
  Srz <- Srmax / 2; Suz <- 0; SWE <- 0
  Te  <- exp(lnTe); td_m <- td / m
  
  # Pre-check that there are enough post-warmup observations to score
  oa  <- obs[(warmup + 1):n]; ok2 <- oa[!is.na(oa)]
  if (length(ok2) < 30) return(1e6)
  om  <- mean(ok2)  # observed mean for NSE denominator
  
  # Accumulate residuals on-the-fly without storing arrays
  ss_r <- 0; ss_o <- 0; nk <- 0L
  
  for (t in seq_len(n)) {
    tv2 <- tavg[t]
    if (tv2 <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0
    } 
    else {
      ma <- snow_mf * (tv2 - snow_t)
      if (ma > SWE) ma <- SWE
      SWE <- SWE - ma; rain <- P[t] + ma
    }
    
    if (Sd >= tts[1L]) { sf <- 0
    } 
    else if (Sd <  tts[nc]) { sf <- tfc[nc]
    } 
    else {
      lo2 <- 1L; hi <- nc
      while (hi - lo2 > 1L) {
        mid <- (lo2 + hi) %/% 2L
        if (tts[mid] > Sd) lo2 <- mid else hi <- mid
      }
      sf <- tfc[lo2]
    }
    
    infil <- (1 - sf) * rain
    Srz   <- Srz - infil
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }
    if (Srz < Srmax) { ae <- pet[t] * (1 - Srz / Srmax) } else { ae <- 0 }
    Srz <- Srz + ae; if (Srz > Srmax) Srz <- Srmax
    
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd)
      if (quz > Suz) quz <- Suz
      Suz <- Suz - quz
    } 
    else { quz <- 0 }
    
    Sd <- Sd - quz; if (Sd < 0) Sd <- 0
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + 50) qb <- Sd + 50; if (qb < 0) qb <- 0
    Sd <- Sd + qb; if (Sd < 0) Sd <- 0
    qsim <- qb + sf * rain
    
    # Only accumulate residuals for post-warmup days with valid observations
    if (t > warmup) {
      ot <- obs[t]
      if (!is.na(ot) && is.finite(qsim)) {
        ss_r <- ss_r + (ot - qsim)^2
        ss_o <- ss_o + (ot - om)^2
        nk   <- nk + 1L
      }
    }
  }
  
  if (nk < 30 || ss_o <= 0) return(1e6)
  -(1 - ss_r / ss_o)  # negative NSE for minimization
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  OPTIMIZER — Nelder-Mead on a user-defined calibration window
#
# Runs optim() with method = "Nelder-Mead" from three different starting
# points to reduce the chance of converging to a local optimum, then returns
# the best (highest NSE) parameter set clamped to valid bounds.
# ═══════════════════════════════════════════════════════════════════════════════
optimize_params <- function(data, topidx, lat_deg, cal_start, cal_end,
                            spinup_days = 365, progress_fn = NULL) {
  
  # Parameter names, lower & upper bounds
  pn        <- c("lnTe", "m", "Srmax", "td", "snow_t", "snow_mf")
  lo.       <- c(-7,  5,   5,   1, -2, 1)
  up        <- c( 5, 100, 300, 60,  2, 6)
  names(lo) <- pn; names(up) <- pn
  
  # Subset data to calibration window; fall back to full record if too short
  cal <- data[data$date >= cal_start & data$date <= cal_end, ]
  if (nrow(cal) < spinup_days + 30) cal <- data
  
  cn <- nrow(cal)
  cP <- cal$precip_mm; ct <- cal$tavg; co <- cal$flow_mm
  cp <- hamon_pet(ct, yday(cal$date), lat_deg)
  cw <- min(spinup_days, cn - 1)
  
  # Objective: negative NSE
  obj <- function(pv) {
    tryCatch(topmodel_nse(pv, cP, ct, co, cp, topidx, cn, cw, lo, up),
             error = function(e) 1e6)
  }
  
  # 3 starting points covering different parts of parameter space
  starts <- list(
    c(lnTe =  2,  m = 23, Srmax =  68, td =  1, snow_t = -0.7, snow_mf = 2.9),
    c(lnTe =  1,  m = 30, Srmax = 100, td = 10, snow_t =  0.0, snow_mf = 3.0),
    c(lnTe = -2,  m = 40, Srmax = 150, td = 20, snow_t = -0.5, snow_mf = 3.5)
  )
  
  bv <- Inf; bp <- starts[[1]]
  
  for (i in seq_along(starts)) {
    
    # Notify the UI of progress if a callback is provided
    if (!is.null(progress_fn))
      progress_fn(i, length(starts), sprintf("Start %d of %d \u2026", i, length(starts)))
    opt <- tryCatch(
      optim(starts[[i]], obj, method = "Nelder-Mead",
            control = list(maxit = 400, reltol = 1e-6)),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value) && opt$value < bv) {
      bv <- opt$value; bp <- opt$par
    }
  }
  
  names(bp) <- pn
  pmax(pmin(bp, up), lo)  # clamp final result to valid bounds
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: compute NSE over an arbitrary date window
#
# Used to display calibration-period and viewing-window NSE separately in the
# plot title, independently of the full-record NSE stored in the result list.
# ═══════════════════════════════════════════════════════════════════════════════
compute_nse <- function(res, start_date, end_date) {
  
  # Only include post-warmup days that fall within the requested window
  vi <- which(res$dates >= start_date & res$dates <= end_date &
                seq_along(res$Qsim) > res$warmup)
  if (length(vi) < 30) return(NA_real_)
  o  <- res$Qobs[vi]; s <- res$Qsim[vi]
  ok <- !is.na(o) & is.finite(s)
  if (sum(ok) < 30) return(NA_real_)
  ss_r <- sum((o[ok] - s[ok])^2)
  ss_t <- sum((o[ok] - mean(o[ok]))^2)
  if (ss_t <= 0) return(NA_real_)
  1 - ss_r / ss_t
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: snap a slider value to the nearest valid step
#
# This helper rounds an optimized value to the nearest slider increment so
# updateSliderInput() doesn't produce a value between two ticks.
# ═══════════════════════════════════════════════════════════════════════════════
snap <- function(val, mn, mx, step) {
  s <- round((val - mn) / step) * step + mn
  max(min(s, mx), mn)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  LOAD HBEF DEFAULTS AT STARTUP
#
# These objects are computed once when the app launches and shared across all
# user sessions (global scope in Shiny = evaluated before server() runs).
# ═══════════════════════════════════════════════════════════════════════════════
message("Loading HBEF defaults ...")
HBEF_DATA    <- load_hbef_data()
HBEF_TOPIDX  <- readRDS(file.path(DATA_DIR, "topidx.rds"))       # TWI class table (16 rows)
HBEF_TWI_MAT <- readRDS(file.path(DATA_DIR, "twi_matrix.rds"))   # TWI raster as matrix
HBEF_TWI_EXT <- readRDS(file.path(DATA_DIR, "twi_extent.rds"))   # c(xmin, xmax, ymin, ymax)
HBEF_DEM_MAT <- readRDS(file.path(DATA_DIR, "dem_matrix.rds"))   # DEM raster as matrix
HBEF_DEM_EXT <- readRDS(file.path(DATA_DIR, "dem_extent.rds"))   # c(xmin, xmax, ymin, ymax)
HBEF_WS_POLY <- readRDS(file.path(DATA_DIR, "ws_boundary.rds"))  # watershed boundary polygon (x, y)
HBEF_TWI_LAM <- readRDS(file.path(DATA_DIR, "twi_lam.rds"))      # area-weighted mean TWI (lambda)
HBEF_LAT     <- 43.56                                            # latitude of HBEF W3 centroid

message(sprintf("  HBEF: %d records: %s to %s",
                nrow(HBEF_DATA), min(HBEF_DATA$date), max(HBEF_DATA$date)))
message("Ready.")

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  SHINY UI
# ═══════════════════════════════════════════════════════════════════════════════
# Custom CSS:
#    - .pg: parameter group card (light grey background)
#    - .action-btn: full-width sidebar buttons
#    - .opt-params: monospaced readout of optimized parameter values
#    - .catchment-status: blue info banner showing the active dataset summary
# ═══════════════════════════════════════════════════════════════════════════════
css <- "
  body { font-family: 'Segoe UI', Tahoma, sans-serif; }
  .pg { background: #f8f9fa; border-radius: 8px;
        padding: 12px 14px 4px 14px; margin-bottom: 10px; }
  .pg h6 { margin-top: 0; color: #495057; font-weight: 600; }
  .action-btn { width: 100%; margin-top: 6px; margin-bottom: 6px; }
  .opt-params { font-size: 12px; color: #495057; margin-top: 4px;
                padding: 6px 8px; background: #e9ecef; border-radius: 6px;
                font-family: monospace; line-height: 1.6; }
  .pg .help-block { font-size: 11px; color: #888;
                    margin-top: 2px; margin-bottom: 2px; line-height: 1.3; }
  .catchment-status { font-size: 11px; padding: 6px 8px; border-radius: 6px;
                      background: #d1ecf1; color: #0c5460; margin-top: 4px; }
"

ui <- page_sidebar(
  title = "TOPMODEL Explorer",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  tags$head(tags$style(HTML(css))),
  
  # ── Left sidebar ─────────────────────────────────────────────────────────────
  sidebar = sidebar(width = 350,
                    # ── Catchment selector ─────────────────────────────────────────────────────
                    div(class = "pg", h6("\U0001F30D Catchment"),
                        # Switch between HBEF defaults and a user-uploaded custom catchment
                        radioButtons("mode", "Data source",
                                     choices  = c("HBEF Watershed 3 (default)" = "hbef",
                                                  "Custom upload"              = "custom"),
                                     selected = "hbef"),
                        
                        # Upload panel (only shown in custom mode)
                        conditionalPanel("input.mode == 'custom'",
                                         helpText("Upload pre-computed catchment data. See format guide on the 'Custom Data Help' tab."),
                                         fileInput("up_ts",     "Time series CSV (required)", accept = ".csv"),
                                         fileInput("up_topidx", "TOPIDX CSV (required)",      accept = ".csv"),
                                         numericInput("up_lat", "Latitude (degrees)", value = 43.56,
                                                      min = -90, max = 90, step = 0.01),
                                         helpText("Optional spatial files for maps:"),
                                         fileInput("up_twi",  "TWI matrix (.rds)",        accept = ".rds"),
                                         fileInput("up_dem",  "DEM matrix (.rds)",        accept = ".rds"),
                                         fileInput("up_poly", "Watershed boundary (.rds)", accept = ".rds"),
                                         actionButton("load_custom", "\u2705 Load Catchment",
                                                      class = "btn-info action-btn")
                        ),
                        
                        # Displays the name, date range, TWI class count, and latitude of the currently active dataset
                        uiOutput("catchment_status")
                        ),
                    
                    # Year slider is built dynamically once the dataset is known
                    uiOutput("yr_slider_ui"),
                    helpText("Adjust the time window shown on the plots."),
                    
                    actionButton("run", "\u25B6  Run Model",        class = "btn-success action-btn"),
                    downloadButton("download_csv", "\u2B07  Export Output (CSV)",
                                   class = "btn-secondary action-btn"),
                    
                    # Displays the exact (un-snapped) optimized parameter values after calibration
                    uiOutput("opt_disp"),
                    
                    # ── Subsurface flow parameters ─────────────────────────────────────────────
                    div(class = "pg", h6("\u2193 Subsurface Flow"),
                        helpText("Log of soil transmissivity at saturation. Higher = more baseflow."),
                        sliderInput("lnTe", "ln(Te) \u2013 Log Soil Transmissivity", -7, 5, 1.9, 0.1),
                        
                        helpText("Decay rate of transmissivity with depth. Small m = flashy response."),
                        sliderInput("m", "m \u2013 Transmissivity Decay (mm)", 5, 100, 19, 1),
                        
                        helpText("Time delay for water to drain through the unsaturated zone to the water table."),
                        sliderInput("td", "td \u2013 Unsaturated Zone Time Delay (days)", 1, 60, 3, 1)
                    ),
                    
                    # ── Root zone & ET parameters ──────────────────────────────────────────────
                    div(class = "pg", h6("\U0001F331 Root Zone & ET"),
                        helpText("Maximum water the root zone can hold before draining. Controls ET — larger = more ET in summer."),
                        sliderInput("Srmax", "Srmax \u2013 Max Root Zone Storage (mm)", 5, 300, 65, 5)
                    ),
                    
                    # ── Snow parameters ────────────────────────────────────────────────────────
                    div(class = "pg", h6("\u2744 Snow"),
                        helpText("Temperature threshold: below this, precipitation falls as snow."),
                        sliderInput("snow_t",  "Snow/Rain Temperature Threshold \u2013 Tmelt (\u00B0C)", -2, 2, -0.7, 0.1),
                        
                        helpText("Degree-day melt factor: mm of snowmelt per degree above Tmelt per day."),
                        sliderInput("snow_mf", "Degree-Day Melt Factor (mm/\u00B0C/day)", 1, 6, 3.4, 0.1)
                    ),
                    
                    # ── Spin-up & calibration window ───────────────────────────────────────────
                    div(class = "pg", h6("\u23F1 Spin-Up & Calibration"),
                        helpText("Spin-up: model warms up for this many days before NSE is scored."),
                        sliderInput("spinup", "Spin-up (days)", 90, 1095, 365, 30),
                        
                        helpText("Calibration period used by the Optimize button."),
                        uiOutput("cal_yr_ui")
                    ),
                    
                    actionButton("opt", "\u26A1 Optimize (Nelder-Mead)", class = "btn-primary action-btn"),
                    helpText("NOTE: Nelder-Mead is a local optimizer — results depend on starting values.
             This optimizer tries 3 starting points and keeps the best NSE result.")
  ),
  
  # ── Main panel: tabbed plots & help ──────────────────────────────────────────
  navset_card_tab(
    nav_panel("Model Output",             plotlyOutput("model_out",    height = "950px")),
    nav_panel("Input Data",               plotlyOutput("input_plot",   height = "640px")),
    nav_panel("Observed vs Simulated",    plotlyOutput("scatter",      height = "520px")),
    nav_panel("Flow Duration",            plotlyOutput("fdc",          height = "520px")),
    nav_panel("TWI Map",                  plotOutput("twi_map",        height = "560px")),
    
    nav_panel("TWI Histogram",
              fluidRow(
                column(4,
                       sliderInput("twi_classes", "Number of TWI classes", 4, 40, 16, 1),
                       helpText("More classes = finer detail."),
                       verbatimTextOutput("twi_stats")),
                column(8, plotOutput("twi_hist", height = "500px"))
              )),
    
    nav_panel("Saturation Map",
              fluidRow(
                column(8,
                       sliderInput("sat_day", "Day", 1, 730, 1, step = 1, width = "100%",
                                   animate = animationOptions(interval = 1000, loop = FALSE)),
                       plotOutput("sat_map", height = "460px")),
                column(4,
                       plotlyOutput("sat_q",   height = "240px"),
                       plotlyOutput("sat_pct", height = "240px"),
                       helpText("Tip: click any point on the discharge plot to jump the map to that day."))
              )),
    
    nav_panel("Custom Data Help",
              div(style = "padding:20px;max-width:800px",
                  h4("Loading a Custom Catchment"),
                  p("To run TOPMODEL on a different catchment, prepare these files locally and upload them in the sidebar."),
                  h5("Required Files"),
                  tags$ol(
                    tags$li(tags$b("Time series CSV"),
                            tags$ul(
                              tags$li("Columns: ", tags$code("date"), ", ", tags$code("precip_mm"),
                                      ", ", tags$code("tavg"), ", ", tags$code("flow_mm")),
                              tags$li("Date format: YYYY-MM-DD"),
                              tags$li("Daily values; missing flow values OK (NA); missing precip filled with 0"),
                              tags$li("Example row: ", tags$code("2010-01-01,5.2,-3.1,1.8"))
                            )),
                    tags$li(tags$b("TOPIDX CSV"),
                            tags$ul(
                              tags$li("Columns: ", tags$code("twi"), ", ", tags$code("area_frac")),
                              tags$li("Each row is one TWI class; area fractions must sum to ~1.0"),
                              tags$li("Typically 20–50 classes is sufficient"),
                              tags$li("Compute locally using WhiteBox + your DEM, then bin into classes")
                            )),
                    tags$li(tags$b("Latitude"),
                            tags$ul(tags$li("Decimal degrees (positive = North). Used by Hamon PET.")))
                  ),
                  h5("Optional Files (for maps)"),
                  tags$ul(
                    tags$li(tags$b("TWI matrix (.rds)"),
                            " — a 2D matrix of TWI values with attribute ",
                            tags$code('attr(., "extent")'), " = c(xmin, xmax, ymin, ymax). Max ",
                            format(CELL_LIMIT, big.mark = ","), " cells."),
                    tags$li(tags$b("DEM matrix (.rds)"),
                            " — same format as TWI matrix; holds elevations in metres."),
                    tags$li(tags$b("Watershed boundary (.rds)"),
                            " — 2-column matrix or data.frame of polygon (x, y) coordinates.")
                  ),
                  h5("R helper to prepare files"),
                  tags$pre(style = "background:#f5f5f5;padding:10px;border-radius:6px;font-size:11px",
                           paste(
                             "# After computing TWI raster with WhiteBox:",
                             "library(terra)",
                             "twi <- rast('twi.tif')",
                             "twi_mat <- as.matrix(twi, wide = TRUE)",
                             "attr(twi_mat, 'extent') <- c(xmin(twi), xmax(twi), ymin(twi), ymax(twi))",
                             "saveRDS(twi_mat, 'twi_matrix.rds')",
                             "",
                             "# TOPIDX classes:",
                             "vals <- na.omit(values(twi))",
                             "h <- hist(vals, breaks = 30, plot = FALSE)",
                             "topidx <- data.frame(twi = h$mids, area_frac = h$counts / sum(h$counts))",
                             "write.csv(topidx, 'topidx.csv', row.names = FALSE)",
                             sep = "\n"
                           ))
              ))
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  SHINY SERVER
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  rv  <- reactiveVal(NULL)  # stores the most recent run_topmodel() result
  otx <- reactiveVal(NULL)  # stores formatted optimized parameter string
  
  # ── Active dataset (reactive) ─────────────────────────────────────────────
  # Returns a named list describing the current catchment:
  #   $data, $topidx, $twi_mat, $twi_ext, $dem_mat, $dem_ext,
  #   $ws_poly, $twi_lam, $lat, $has_spatial, $name
  custom_data <- reactiveVal(NULL)
  
  ds <- reactive({
    if (input$mode == "hbef" || is.null(custom_data())) {
      # Default: return the pre-loaded HBEF objects
      list(name = "HBEF Watershed 3",
           data     = HBEF_DATA,    topidx   = HBEF_TOPIDX,
           twi_mat  = HBEF_TWI_MAT, twi_ext  = HBEF_TWI_EXT,
           dem_mat  = HBEF_DEM_MAT, dem_ext  = HBEF_DEM_EXT,
           ws_poly  = HBEF_WS_POLY, twi_lam  = HBEF_TWI_LAM,
           lat      = HBEF_LAT,     has_spatial = TRUE)
    } 
    else {
      custom_data()
    }
  })
  
  # ── Load custom catchment when the button is clicked ──────────────────────
  observeEvent(input$load_custom, {
    req(input$up_ts, input$up_topidx)
    
    # 1. Read and validate the time series CSV
    df <- tryCatch(read.csv(input$up_ts$datapath, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || !all(c("date", "precip_mm", "tavg", "flow_mm") %in% names(df))) {
      showNotification("Time series CSV must have columns: date, precip_mm, tavg, flow_mm",
                       type = "error", duration = 8)
      return()
    }
    df$date <- as.Date(df$date)
    df <- df[order(df$date), ]
    df$precip_mm[is.na(df$precip_mm)] <- 0
    
    # Linearly interpolate missing temperature values
    if (any(is.na(df$tavg))) {
      idx <- which(!is.na(df$tavg))
      df$tavg <- approx(idx, df$tavg[idx], seq_len(nrow(df)), rule = 2)$y
    }
    
    # 2. Read and validate the TOPIDX CSV; normalize area fractions to sum to 1
    ti <- tryCatch(read.csv(input$up_topidx$datapath, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(ti) || !all(c("twi", "area_frac") %in% names(ti))) {
      showNotification("TOPIDX CSV must have columns: twi, area_frac",
                       type = "error", duration = 8)
      return()
    }
    af_sum <- sum(ti$area_frac)
    if (af_sum < 0.95 || af_sum > 1.05)
      showNotification(sprintf("TOPIDX area_frac sums to %.3f (should be ~1.0)", af_sum),
                       type = "warning", duration = 6)
    ti$area_frac  <- ti$area_frac / af_sum   # normalize
    topidx_mat    <- as.matrix(ti[, c("twi", "area_frac")])
    
    # 3. Optional spatial files — fail gracefully if malformed or oversized
    has_spatial <- FALSE
    twi_mat <- NULL; twi_ext <- NULL
    dem_mat <- NULL; dem_ext <- NULL
    ws_poly <- NULL
    twi_lam <- sum(topidx_mat[, 1] * topidx_mat[, 2])  # lambda from TOPIDX classes
    
    if (!is.null(input$up_twi)) {
      twi_mat <- tryCatch(readRDS(input$up_twi$datapath), error = function(e) NULL)
      if (!is.null(twi_mat)) {
        if (length(twi_mat) > CELL_LIMIT) {
          showNotification(sprintf("TWI matrix has %s cells (limit %s). Skipping maps.",
                                   format(length(twi_mat), big.mark = ","),
                                   format(CELL_LIMIT,      big.mark = ",")),
                           type = "warning", duration = 8)
          twi_mat <- NULL
        } 
        else {
          twi_ext <- attr(twi_mat, "extent")
          if (is.null(twi_ext)) {
            showNotification("TWI matrix missing 'extent' attribute. Skipping maps.",
                             type = "warning", duration = 6)
            twi_mat <- NULL
          } 
          else {
            names(twi_ext) <- c("xmin", "xmax", "ymin", "ymax")
          }
        }
      }
    }
    
    # Only attempt to read DEM if the TWI matrix loaded successfully
    if (!is.null(input$up_dem) && !is.null(twi_mat)) {
      dem_mat <- tryCatch(readRDS(input$up_dem$datapath), error = function(e) NULL)
      if (!is.null(dem_mat)) {
        dem_ext <- attr(dem_mat, "extent")
        if (!is.null(dem_ext)) names(dem_ext) <- c("xmin", "xmax", "ymin", "ymax")
      }
    }
    if (!is.null(input$up_poly))
      ws_poly <- tryCatch(readRDS(input$up_poly$datapath), error = function(e) NULL)
    
    has_spatial <- !is.null(twi_mat)
    
    custom_data(list(
      name        = sprintf("Custom (%s)", basename(input$up_ts$name)),
      data        = df,          topidx      = topidx_mat,
      twi_mat     = twi_mat,     twi_ext     = twi_ext,
      dem_mat     = dem_mat,     dem_ext     = dem_ext,
      ws_poly     = ws_poly,     twi_lam     = twi_lam,
      lat         = input$up_lat, has_spatial = has_spatial
    ))
    
    showNotification(sprintf("Loaded: %d days, %d TWI classes%s",
                             nrow(df), nrow(topidx_mat),
                             if (has_spatial) ", spatial maps available" else ", no spatial maps"),
                     type = "message", duration = 5)
    rv(NULL)   # clear old results so stale plots don't show
  })
  
  # ── Catchment status banner ────────────────────────────────────────────────
  output$catchment_status <- renderUI({
    d <- ds()
    div(class = "catchment-status",
        tags$b("Active: "), d$name, tags$br(),
        sprintf("%d days  |  %d TWI classes  |  lat %.2f\u00B0",
                nrow(d$data), nrow(d$topidx), d$lat))
  })
  
  # ── Year sliders (driven by the active dataset's date range) ──────────────
  yr_range <- reactive({
    d <- ds()$data
    c(as.integer(format(min(d$date), "%Y")),
      as.integer(format(max(d$date), "%Y")))
  })
  
  output$yr_slider_ui <- renderUI({
    yr <- yr_range()
    sliderInput("yr", "Viewing Period",
                min = yr[1], max = yr[2], value = yr,
                step = 1, sep = "", ticks = TRUE)
  })
  
  output$cal_yr_ui <- renderUI({
    yr <- yr_range()
    sliderInput("cal_yr", "Calibration Years",
                min = yr[1], max = yr[2], value = yr,
                step = 1, sep = "", ticks = TRUE)
  })
  
  # Convenience reactive: start/end Date objects for the viewing window
  view_range <- reactive({
    req(input$yr)
    list(start = as.Date(paste0(input$yr[1], "-01-01")),
         end   = as.Date(paste0(input$yr[2], "-12-31")))
  })
  
  # Bundle current slider values into a named parameter vector
  cpars <- reactive(c(lnTe    = input$lnTe,    m       = input$m,
                      Srmax   = input$Srmax,   td      = input$td,
                      snow_t  = input$snow_t,  snow_mf = input$snow_mf))
  
  # ── Run Model button ───────────────────────────────────────────────────────
  observeEvent(input$run, {
    d <- ds()
    res <- tryCatch(
      run_topmodel(cpars(), d$data, d$topidx, d$lat, spinup_days = input$spinup),
      error = function(e) { showNotification(e$message, type = "error"); NULL }
    )
    if (!is.null(res)) {
      rv(res)
      # Compute NSE over the calibration window for the notification message
      cs  <- as.Date(paste0(input$cal_yr[1], "-01-01"))
      ce  <- as.Date(paste0(input$cal_yr[2], "-12-31"))
      nc  <- compute_nse(res, cs, ce)
      nse_str <- if (is.finite(nc)) sprintf("%.3f", nc) else "N/A"
      showNotification(sprintf("NSE (cal) = %s", nse_str), type = "message", duration = 3)
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)
  
  # ── Auto-run on first load so plots are not blank at startup ──────────────
  once <- reactiveVal(FALSE)
  observe({
    req(!once())
    d <- ds()
    res <- tryCatch(
      run_topmodel(isolate(cpars()), d$data, d$topidx, d$lat,
                   spinup_days = isolate(input$spinup) %||% 365),
      error = function(e) NULL
    )
    if (!is.null(res)) rv(res)
    once(TRUE)
  })
  
  # ── Re-run automatically when the user switches datasets ──────────────────
  observeEvent(ds()$name, {
    d <- ds()
    res <- tryCatch(
      run_topmodel(isolate(cpars()), d$data, d$topidx, d$lat,
                   spinup_days = isolate(input$spinup) %||% 365),
      error = function(e) NULL
    )
    if (!is.null(res)) rv(res)
  }, ignoreInit = TRUE)
  
  # ── Optimize button ────────────────────────────────────────────────────────
  observeEvent(input$opt, {
    d   <- ds()
    cs  <- as.Date(paste0(input$cal_yr[1], "-01-01"))
    ce  <- as.Date(paste0(input$cal_yr[2], "-12-31"))
    
    showNotification("Optimizing \u2026 Start 1 of 3",
                     type = "message", duration = NULL, id = "op")
    
    best <- optimize_params(d$data, d$topidx, d$lat, cs, ce,
                            spinup_days = input$spinup,
                            progress_fn = function(i, n, msg)
                              showNotification(sprintf("Optimizing \u2026 %s", msg),
                                               type = "message", duration = NULL, id = "op"))
    removeNotification("op")
    
    # Push optimized values back to the sliders (snapped to step increments)
    updateSliderInput(session, "lnTe",    value = snap(best["lnTe"],    -7,  5,   0.1))
    updateSliderInput(session, "m",       value = snap(best["m"],        5,  100, 1))
    updateSliderInput(session, "td",      value = snap(best["td"],       1,  60,  1))
    updateSliderInput(session, "Srmax",   value = snap(best["Srmax"],    5,  300, 5))
    updateSliderInput(session, "snow_t",  value = snap(best["snow_t"],  -2,  2,   0.1))
    updateSliderInput(session, "snow_mf", value = snap(best["snow_mf"],  1,  6,   0.1))
    
    # Store the exact (un-snapped) optimized values for display
    otx(sprintf("ln(Te)=%.2f  m=%.1f  Srmax=%.1f\nUnsat. Delay td=%.1f  Tmelt=%.2f  Melt Coeff=%.2f",
                best["lnTe"], best["m"], best["Srmax"],
                best["td"],   best["snow_t"], best["snow_mf"]))
    
    # Re-run the model with the exact optimized parameters
    res <- tryCatch(
      run_topmodel(best, d$data, d$topidx, d$lat, spinup_days = input$spinup),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      rv(res)
      nc  <- compute_nse(res, cs, ce)
      nse_str <- if (is.finite(nc)) sprintf("%.3f", nc) else "N/A"
      showNotification(sprintf("Optimized! NSE (cal) = %s", nse_str),
                       type = "message", duration = 5)
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # ── Optimized parameter readout (below Run button) ─────────────────────────
  output$opt_disp <- renderUI({
    txt <- otx()
    if (is.null(txt)) return(NULL)
    div(class = "opt-params",
        tags$strong("Optimized values:"), tags$br(),
        tags$pre(style = "margin:2px 0 0 0; white-space:pre-wrap;", txt))
  })
  
  # ── CSV export ────────────────────────────────────────────────────────────
  output$download_csv <- downloadHandler(
    filename = function() sprintf("topmodel_output_%s.csv", format(Sys.Date(), "%Y%m%d")),
    content  = function(file) {
      res <- rv(); req(res)
      d   <- ds()$data
      out <- data.frame(
        date                 = res$dates,
        P_mm                 = d$precip_mm,
        T_C                  = d$tavg,
        Qobs_mm              = res$Qobs,
        Qsim_mm              = res$Qsim,
        Qbaseflow_mm         = res$Qb,
        Qsat_excess_mm       = res$Qof,
        ET_mm                = res$ET,
        SWE_mm               = res$SWE,
        Storage_Deficit_mm   = res$Sd,
        Saturated_Fraction   = res$SF
      )
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  # ════════════════════════════════════════════════════════════════════════════
  # MODEL OUTPUT — 7-panel stacked time series
  #
  # Panels (top to bottom):
  #   1. Observed & simulated discharge
  #   2. Saturation-excess overland flow (surface runoff)
  #   3. Baseflow (subsurface)
  #   4. Actual evapotranspiration (ET)
  #   5. Snow water equivalent (SWE)
  #   6. Watershed storage deficit (Sd)
  #   7. Saturated catchment fraction (%)
  #
  # The x-axis is restricted to the user's viewing window; NSE is reported
  # for both the calibration period and the viewing window in the plot title.
  # ════════════════════════════════════════════════════════════════════════════
  output$model_out <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    
    d <- data.frame(date = res$dates, Qobs = res$Qobs, Qsim = res$Qsim,
                    Qof  = res$Qof,   Qb   = res$Qb,   ET   = res$ET,
                    SWE  = res$SWE,   Sd   = res$Sd,   SF   = res$SF * 100)
    
    vr  <- view_range()
    cs  <- as.Date(paste0(input$cal_yr[1], "-01-01"))
    ce  <- as.Date(paste0(input$cal_yr[2], "-12-31"))
    nse_cal  <- compute_nse(res, cs, ce)
    nse_view <- compute_nse(res, vr$start, vr$end)
    nse_cal_str  <- if (is.finite(nse_cal))  sprintf("%.3f", nse_cal)  else "N/A"
    nse_view_str <- if (is.finite(nse_view)) sprintf("%.3f", nse_view) else "N/A"
    
    # Shared y-axis style: compact with auto-margins for label room
    yax <- list(title = "", nticks = 4, automargin = TRUE,
                titlefont = list(size = 13), tickfont = list(size = 12))
    
    # Panel 1 — Discharge
    p1 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qobs, name = "Observed Discharge",
                line = list(color = "#2980b9", width = 1.5)) |> 
      add_lines(y = ~Qsim, name = "Simulated Discharge",
                line = list(color = "#e74c3c", width = 1.5)) |> 
      layout(yaxis = yax)
    
    # Panel 2 — Saturation-excess overland flow
    p2 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qof, name = "Saturation-Excess Overland Flow",
                line = list(color = "#e67e22", width = 1.2)) |> 
      layout(yaxis = yax)
    
    # Panel 3 — Baseflow
    p3 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qb, name = "Baseflow (Subsurface)",
                line = list(color = "#8e44ad", width = 1.2)) |> 
      layout(yaxis = yax)
    
    # Panel 4 — Actual ET
    p4 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~ET, name = "Actual Evapotranspiration",
                line = list(color = "#27ae60", width = 1.2)) |> 
      layout(yaxis = yax)
    
    # Panel 5 — Snow Water Equivalent
    p5 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~SWE, name = "Snow Water Equivalent",
                line = list(color = "#3498db", width = 1.2)) |> 
      layout(yaxis = yax)
    
    # Panel 6 — Watershed storage deficit
    p6 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~Sd, name = "Watershed Storage Deficit",
                line = list(color = "#c0392b", width = 1.2)) |>
      layout(yaxis = yax)
    
    # Panel 7 — Saturated catchment fraction
    p7 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~SF, name = "Saturated Catchment Area",
                line = list(color = "#2166ac", width = 1.2)) |> 
      layout(yaxis = yax, xaxis = list(title = ""))
    
    # Panel height proportions and y-axis label annotations
    ht <- c(.22, .11, .11, .11, .15, .15, .15)
    lb <- c("Discharge (mm/d)",
            "Saturation-Excess Overland Flow (mm/d)",
            "Baseflow (mm/d)",
            "Actual ET (mm/d)",
            "Snow Water Equivalent (mm)",
            "Watershed Storage Deficit (mm)",
            "Saturated Catchment Area (%)")
    ch   <- cumsum(ht)
    tops <- 1 - c(0, ch[-length(ch)])
    
    # Place a bold label at the top-left of each panel
    an <- lapply(seq_along(lb), function(i)
      list(text     = paste0("<b>", lb[i], "</b>"),
           x = 0.01, y = tops[i] - 0.005,
           xref = "paper", yref = "paper",
           xanchor = "left", yanchor = "top",
           showarrow = FALSE,
           font      = list(size = 12, color = "#333"),
           bgcolor   = "rgba(255,255,255,0.8)", borderpad = 2))
    
    subplot(p1, p2, p3, p4, p5, p6, p7,
            nrows = 7, shareX = TRUE, titleY = FALSE, heights = ht) |> 
      layout(
        title      = list(text = sprintf("Model Output  |  NSE (cal) = %s  |  NSE (view) = %s",
                                         nse_cal_str, nse_view_str),
                          font = list(size = 16)),
        showlegend = FALSE,
        hovermode  = "x unified",
        annotations = an,
        
        # Restrict x-axis to the viewing window selected in the sidebar
        xaxis      = list(range = as.character(c(vr$start, vr$end))),
        margin     = list(t = 55, b = 30, l = 70)
      )
  })
  
  # ════════════════════════════════════════════════════════════════════════════
  # INPUT DATA — precipitation, temperature, and observed streamflow
  #
  # Three stacked panels showing the raw forcing data.  Useful for sanity-
  # checking a newly uploaded custom catchment before running the model.
  # ════════════════════════════════════════════════════════════════════════════
  output$input_plot <- renderPlotly({
    d  <- ds()$data
    vr <- view_range()
    
    p1 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~precip_mm, line = list(color = "#2980b9", width = 1)) |>
      layout(yaxis = list(title = "Precip (mm/d)", nticks = 4, automargin = TRUE,
                          titlefont = list(size = 13), tickfont = list(size = 12)))
    
    p2 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~tavg, line = list(color = "#e67e22", width = 1)) |>
      layout(yaxis = list(title = "Air Temp (\u00B0C)", nticks = 4, automargin = TRUE,
                          titlefont = list(size = 13), tickfont = list(size = 12)))
    
    p3 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~flow_mm, line = list(color = "#27ae60", width = 1)) |>
      layout(yaxis = list(title = "Obs Q (mm/d)", nticks = 4, automargin = TRUE,
                          titlefont = list(size = 13), tickfont = list(size = 12)))
    
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE) |>
      layout(title      = list(text = "Input Data: Precipitation, Temperature, Observed Discharge",
                               font = list(size = 15)),
             showlegend = FALSE,
             hovermode  = "x unified",
             xaxis      = list(range = as.character(c(vr$start, vr$end))),
             margin     = list(t = 50, l = 90))
  })
  
  # ════════════════════════════════════════════════════════════════════════════
  # OBSERVED vs SIMULATED — scatter plot
  #
  # Only post-warmup days within the current viewing window are plotted.
  # The dashed 1:1 line marks where a perfect model would fall.
  # ════════════════════════════════════════════════════════════════════════════
  output$scatter <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    vr  <- view_range()
    vi  <- which(res$dates >= vr$start & res$dates <= vr$end &
                   seq_along(res$Qsim) > res$warmup)
    d   <- data.frame(obs = res$Qobs[vi], sim = res$Qsim[vi])
    d   <- d[!is.na(d$obs) & is.finite(d$sim), ]
    validate(need(nrow(d) > 0, "No valid data in the selected viewing window."))
    rng <- range(c(d$obs, d$sim), na.rm = TRUE)
    
    plot_ly(d, x = ~obs, y = ~sim, type = "scattergl", mode = "markers",
            marker = list(size = 3, color = "#3498db", opacity = 0.4)) |>
      add_lines(x = rng, y = rng, name = "1:1 Line",
                line = list(color = "black", dash = "dash", width = 1.5)) |>
      layout(title  = list(text = "Observed vs Simulated Discharge (viewing window)",
                           font = list(size = 15)),
             xaxis  = list(title = "Observed Discharge (mm/d)",
                           titlefont = list(size = 13), tickfont = list(size = 12)),
             yaxis  = list(title = "Simulated Discharge (mm/d)",
                           titlefont = list(size = 13), tickfont = list(size = 12)),
             showlegend = FALSE,
             margin = list(t = 55, l = 70, b = 60))
  })
  
  # ════════════════════════════════════════════════════════════════════════════
  # FLOW DURATION CURVE — log-scale exceedance plot
  #
  # Shows the frequency distribution of daily flows on a log y-axis.
  # Good fit across the full range (high, mid, and low flows) indicates a
  # well-calibrated model.
  # ════════════════════════════════════════════════════════════════════════════
  output$fdc <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    vr  <- view_range()
    vi  <- which(res$dates >= vr$start & res$dates <= vr$end &
                   seq_along(res$Qsim) > res$warmup)
    ov  <- sort(res$Qobs[vi][!is.na(res$Qobs[vi])],   decreasing = TRUE)
    sv  <- sort(res$Qsim[vi][is.finite(res$Qsim[vi])], decreasing = TRUE)
    eo  <- seq_along(ov) / length(ov) * 100
    es  <- seq_along(sv) / length(sv) * 100
    
    plot_ly() |> 
      add_lines(x = eo, y = ov, name = "Observed",
                line = list(color = "#2980b9", width = 2)) |> 
      add_lines(x = es, y = sv, name = "Simulated",
                line = list(color = "#e74c3c", width = 2)) |> 
      layout(title  = list(text = "Flow Duration Curve (viewing window)",
                           font = list(size = 15)),
             xaxis  = list(title = "Exceedance Probability (%)",
                           titlefont = list(size = 13), tickfont = list(size = 12)),
             yaxis  = list(title = "Discharge (mm/d)", type = "log",
                           titlefont = list(size = 13), tickfont = list(size = 12)),
             legend = list(orientation = "h", y = -0.15, font = list(size = 13)),
             margin = list(t = 55, l = 70, b = 65))
  })
  
  # ════════════════════════════════════════════════════════════════════════════
  # TWI MAP — static spatial wetness map
  #
  # Displays the Topographic Wetness Index raster as a colour image with the
  # watershed boundary polygon and elevation contours overlaid.  The map is
  # identical regardless of which year or parameters are selected — it reflects
  # only watershed topography, not model state.
  # ════════════════════════════════════════════════════════════════════════════
  output$twi_map <- renderPlot({
    d <- ds()
    if (!d$has_spatial) {
      plot.new()
      title(main = "Spatial maps not available for this catchment")
      text(0.5, 0.5, "Upload a TWI matrix .rds file with an extent attribute to enable maps.",
           cex = 0.9, col = "#555")
      return()
    }
    # Build x/y coordinate vectors from the stored extent attribute
    twi_x   <- seq(d$twi_ext["xmin"], d$twi_ext["xmax"], length.out = ncol(d$twi_mat))
    twi_y   <- seq(d$twi_ext["ymin"], d$twi_ext["ymax"], length.out = nrow(d$twi_mat))
    # Flip the matrix vertically so north is up
    twi_img <- t(d$twi_mat[nrow(d$twi_mat):1, ])
    
    par(mar = c(2, 2, 3, 4), bg = "white")
    image(twi_x, twi_y, twi_img,
          col   = hcl.colors(50, "YlGnBu", rev = TRUE),
          asp   = 1, axes = FALSE, xlab = "", ylab = "",
          main  = sprintf("Topographic Wetness Index  (\u03BB = %.1f)", d$twi_lam),
          cex.main = 1.3)
    
    # Watershed boundary polygon
    if (!is.null(d$ws_poly))
      polygon(d$ws_poly[, 1], d$ws_poly[, 2], border = "#222", lwd = 2.5, col = NA)
    
    # Elevation contour lines (DEM, if available)
    if (!is.null(d$dem_mat) && !is.null(d$dem_ext)) {
      dem_x   <- seq(d$dem_ext["xmin"], d$dem_ext["xmax"], length.out = ncol(d$dem_mat))
      dem_y   <- seq(d$dem_ext["ymin"], d$dem_ext["ymax"], length.out = nrow(d$dem_mat))
      dem_img <- t(d$dem_mat[nrow(d$dem_mat):1, ])
      contour(dem_x, dem_y, dem_img, add = TRUE, nlevels = 10,
              col = adjustcolor("#333", alpha.f = 0.5),
              lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
      mtext("Contour elevations in meters", side = 1, line = 0.5,
            cex = 0.8, col = "#666")
    }
  }, res = 120)
  
  # ════════════════════════════════════════════════════════════════════════════
  # TWI HISTOGRAM — distribution of TWI values
  #
  # If spatial data are available, the histogram is computed from all raster
  # cells; otherwise the TOPIDX class midpoints are used (weighted sampling).
  # ════════════════════════════════════════════════════════════════════════════
  output$twi_hist <- renderPlot({
    d <- ds()
    vals <- if (d$has_spatial) {
      v <- as.vector(d$twi_mat); v[!is.na(v)]
    } else {
      # Reconstruct approximate distribution from TOPIDX class weights
      rep(d$topidx[, 1], times = round(d$topidx[, 2] * 10000))
    }
    par(mar = c(4, 4, 3, 1), bg = "white")
    hist(vals, breaks = input$twi_classes,
         col = "#3498db", border = "white",
         main = sprintf("Distribution of TWI values (%d classes)", input$twi_classes),
         xlab = "TWI",
         ylab = if (d$has_spatial) "Cell count" else "Weighted count",
         cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.1)
    abline(v = mean(vals),   col = "#e74c3c", lwd = 2, lty = 2)
    abline(v = median(vals), col = "#27ae60", lwd = 2, lty = 2)
    legend("topright",
           legend = c(sprintf("Mean = %.2f", mean(vals)),
                      sprintf("Median = %.2f", median(vals))),
           col = c("#e74c3c", "#27ae60"), lwd = 2, lty = 2, bty = "n")
  }, res = 110)
  
  # Summary statistics panel beside the histogram
  output$twi_stats <- renderText({
    d <- ds()
    vals <- if (d$has_spatial) {
      v <- as.vector(d$twi_mat); v[!is.na(v)]
    } else {
      rep(d$topidx[, 1], times = round(d$topidx[, 2] * 10000))
    }
    sprintf("TWI Summary\n  Min:    %.2f\n  Max:    %.2f\n  Mean:   %.2f\n  Median: %.2f\n  SD:     %.2f\n  N:      %d\n  Lambda: %.2f",
            min(vals), max(vals), mean(vals), median(vals), sd(vals),
            length(vals), d$twi_lam)
  })
  
  # ════════════════════════════════════════════════════════════════════════════
  # SATURATION MAP — animated daily map of saturated/unsaturated cells
  #
  # For each selected day the watershed-average deficit (Sd) is converted to a
  # TWI threshold: cells with TWI ≥ (lambda + Sd/m) are predicted to have their
  # water table at or above the surface (saturated).  The slider can be animated
  # automatically or the user can click points on the discharge plot to jump
  # directly to any date of interest.
  # ════════════════════════════════════════════════════════════════════════════
  
  # Indices of model output days that fall within the current viewing window
  view_idx <- reactive({
    res <- rv(); req(!is.null(res))
    vr  <- view_range()
    which(res$dates >= vr$start & res$dates <= vr$end)
  })
  
  # Keep the day slider in sync with the viewing window length
  observe({
    vi <- view_idx(); req(length(vi) > 0)
    updateSliderInput(session, "sat_day", min = 1, max = length(vi), value = 1)
  })
  
  # Jump the map slider when the user clicks on the discharge mini-plot
  observeEvent(event_data("plotly_click", source = "sat_q_src"), {
    cd <- event_data("plotly_click", source = "sat_q_src")
    if (is.null(cd)) return()
    clicked_date <- as.Date(cd$x)
    vi  <- view_idx(); req(length(vi) > 0)
    res <- rv()
    # Find the nearest date in the viewing window
    pos <- which.min(abs(as.numeric(res$dates[vi] - clicked_date)))
    updateSliderInput(session, "sat_day", value = pos)
  })
  
  # Render the binary saturated/unsaturated raster for the selected day
  output$sat_map <- renderPlot({
    res <- rv()
    validate(need(!is.null(res), "Run the model first."))
    d   <- ds()
    if (!d$has_spatial) {
      plot.new(); title(main = "Saturation map not available")
      text(0.5, 0.5, "Upload a TWI matrix .rds file to enable saturation maps.",
           cex = 0.9, col = "#555")
      return()
    }
    vi  <- view_idx(); req(length(vi) > 0)
    si  <- min(input$sat_day %||% 1, length(vi))
    di  <- vi[si]   # index into the full model result arrays
    
    # Compute the TWI threshold for this day's saturation deficit
    thresh <- res$lam + res$Sd[di] / res$m
    label  <- format(res$dates[di], "%Y-%m-%d")
    
    # Classify every raster cell as saturated (1) or unsaturated (0)
    twi_vals  <- as.vector(d$twi_mat)
    twi_valid <- !is.na(twi_vals)
    sat_v     <- rep(NA_real_, length(twi_vals))
    sat_v[twi_valid] <- ifelse(twi_vals[twi_valid] >= thresh, 1, 0)
    sat_mat   <- matrix(sat_v, nrow = nrow(d$twi_mat), ncol = ncol(d$twi_mat))
    sat_img   <- t(sat_mat[nrow(sat_mat):1, ])
    
    twi_x <- seq(d$twi_ext["xmin"], d$twi_ext["xmax"], length.out = ncol(d$twi_mat))
    twi_y <- seq(d$twi_ext["ymin"], d$twi_ext["ymax"], length.out = nrow(d$twi_mat))
    pct   <- round(100 * sum(sat_v == 1, na.rm = TRUE) / sum(twi_valid), 1)
    
    par(mar = c(2, 2, 3, 1), bg = "white")
    image(twi_x, twi_y, sat_img,
          col    = c("#f5e6c8", "#2166ac"),
          breaks = c(-0.5, 0.5, 1.5),
          asp = 1, axes = FALSE, xlab = "", ylab = "",
          main = sprintf("%s \u2014 %.1f%% saturated", label, pct),
          cex.main = 1.2)
    
    if (!is.null(d$ws_poly))
      polygon(d$ws_poly[, 1], d$ws_poly[, 2], border = "#222", lwd = 2.5, col = NA)
    
    # Elevation contours for terrain context
    if (!is.null(d$dem_mat) && !is.null(d$dem_ext)) {
      dem_x   <- seq(d$dem_ext["xmin"], d$dem_ext["xmax"], length.out = ncol(d$dem_mat))
      dem_y   <- seq(d$dem_ext["ymin"], d$dem_ext["ymax"], length.out = nrow(d$dem_mat))
      dem_img <- t(d$dem_mat[nrow(d$dem_mat):1, ])
      contour(dem_x, dem_y, dem_img, add = TRUE, nlevels = 10,
              col = adjustcolor("#555", alpha.f = 0.5),
              lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
    }
    legend("bottomright",
           legend = c("Saturated", "Unsaturated"),
           fill   = c("#2166ac", "#f5e6c8"),
           border = c("#2166ac", "#f5e6c8"),
           bty = "n", cex = 0.9)
  }, res = 110)
  
  # ── Saturation panel: discharge time series (click to jump map) ───────────
  output$sat_q <- renderPlotly({
    res <- rv(); validate(need(!is.null(res), ""))
    vi  <- view_idx(); req(length(vi) > 0)
    si  <- min(input$sat_day %||% 1, length(vi))
    di  <- vi[si]
    d   <- data.frame(date = res$dates, Qobs = res$Qobs, Qsim = res$Qsim)
    
    p <- plot_ly(d, x = ~date, source = "sat_q_src") |> 
      add_lines(y = ~Qobs, name = "Observed",
                line = list(color = "#2980b9", width = 1.5),
                hovertemplate = "%{x|%Y-%m-%d}<br>Obs: %{y:.2f} mm/d<extra></extra>") |> 
      add_lines(y = ~Qsim, name = "Simulated",
                line = list(color = "#e74c3c", width = 1.5),
                hovertemplate = "%{x|%Y-%m-%d}<br>Sim: %{y:.2f} mm/d<extra></extra>") |> 
      add_markers(x = res$dates[di], y = res$Qsim[di], name = "Current Day",
                  marker = list(color = "black", size = 10),
                  hovertemplate = paste0("<b>", format(res$dates[di], "%Y-%m-%d"),
                                         "</b><extra></extra>")) |> 
      layout(
        title      = list(text = "Discharge (mm/d) \u2014 click to jump map",
                          font = list(size = 13)),
        yaxis      = list(title = "", titlefont = list(size = 12), tickfont = list(size = 12)),
        xaxis      = list(title = "", tickfont = list(size = 12),
                          range = as.character(c(view_range()$start, view_range()$end))),
        showlegend = FALSE,
        margin     = list(t = 40, b = 20, l = 45, r = 10),
        hovermode  = "closest"
      )
    event_register(p, "plotly_click")
  })
  
  # ── Saturation panel: saturated fraction time series ──────────────────────
  output$sat_pct <- renderPlotly({
    res <- rv(); validate(need(!is.null(res), ""))
    vi  <- view_idx(); req(length(vi) > 0)
    si  <- min(input$sat_day %||% 1, length(vi))
    di  <- vi[si]
    d   <- data.frame(date = res$dates, SF = res$SF * 100)
    
    plot_ly(d, x = ~date) |> 
      add_lines(y = ~SF, name = "Saturated %",
                line = list(color = "#2166ac", width = 1.5),
                hovertemplate = "%{x|%Y-%m-%d}<br>Sat: %{y:.1f}%<extra></extra>") |> 
      add_markers(x = res$dates[di], y = d$SF[di], name = "Current Day",
                  marker = list(color = "black", size = 10),
                  hovertemplate = paste0("<b>", format(res$dates[di], "%Y-%m-%d"),
                                         "</b><extra></extra>")) |> 
      layout(
        title      = list(text = "Saturated Catchment (%)",
                          font = list(size = 13)),
        yaxis      = list(title = "", titlefont = list(size = 12), tickfont = list(size = 12)),
        xaxis      = list(title = "", tickfont = list(size = 12),
                          range = as.character(c(view_range()$start, view_range()$end))),
        showlegend = FALSE,
        margin     = list(t = 40, b = 20, l = 45, r = 10),
        hovermode  = "closest"
      )
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  LAUNCH
# ═══════════════════════════════════════════════════════════════════════════════
shinyApp(ui, server)
