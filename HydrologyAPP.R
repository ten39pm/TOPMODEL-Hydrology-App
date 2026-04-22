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
# Constants
DATA_DIR <- "Data"  # all input data files go here
HBEF_LAT  <- 43.56    # latitude of Hubbard Brook (used by Hamon PET)
SD_CAP    <- 2000     # mm, upper limit on storage deficit so it can't blow up
QB_BUFFER <- 50       # mm/d, how much baseflow can exceed Sd before we cap it
MIN_NSE_N <- 30       # need at least 30 valid obs days before computing NSE

PARAM_CONFIG <- list(
  lnTe    = list(lo = -7, up =   5, default =  1.9, step = 0.1,
                 label = "Log of Soil Transmissivity",
                 help  = "Logarithm of soil transmissivity at saturation. Higher = more baseflow."),
  m       = list(lo =  5, up = 100, default = 19,   step = 1,
                 label = "Transmissivity decay with Depth (mm)",
                 help  = "Controls how fast transmissivity declines with depth. Small = flashy, large = slow."),
  Srmax   = list(lo =  5, up = 300, default = 65,   step = 5,
                 label = "Maximum root zone storage capacity (mm)",
                 help  = "Controls how much water plants can access for ET."),
  td      = list(lo =  1, up =  60, default =  3,   step = 1,
                 label = "Unsaturated Zone Time Delay (days)",
                 help  = "Time for water to drain through the unsaturated zone to the water table."),
  snow_t  = list(lo = -2, up =   2, default = -0.7, step = 0.1,
                 label = "Snow/Rain Temperature Threshold (°C)",
                 help  = "Below this, precipitation falls as snow instead of rain."),
  snow_mf = list(lo =  1, up =   6, default =  3.4, step = 0.1,
                 label = "Melt Factor (mm/°C/d)",
                 help  = "Degree-day melt factor: mm of snowmelt per degree above Tmelt per day.")
)

# Pull the bounds out as named vectors so optim() can use them directly.
PARAM_NAMES <- names(PARAM_CONFIG)
PARAM_LO    <- setNames(sapply(PARAM_CONFIG, `[[`, "lo"),      PARAM_NAMES)
PARAM_UP    <- setNames(sapply(PARAM_CONFIG, `[[`, "up"),      PARAM_NAMES)
PARAM_DEF   <- setNames(sapply(PARAM_CONFIG, `[[`, "default"), PARAM_NAMES)
PARAM_STEP  <- setNames(sapply(PARAM_CONFIG, `[[`, "step"),    PARAM_NAMES)
# ═══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD HBEF DATA
#
# Reads in & merges the 3 HBEF input time series:
#   - Precipitation:   DailyWatershed.csv
#   - Streamflow:      HBEF_DailyStreamflow_1956-2024.csv
#   - Air Temperature: HBEF_air_temp_daily.csv
# ═══════════════════════════════════════════════════════════════════════════════
load_hbef_data <- function() {
  # ── Daily Watershed Precipitation ────────────────────────────────────────────
  precip <- read.csv(file.path(DATA_DIR, "DailyWatershed.csv"),
                     stringsAsFactors = FALSE)
  precip$DATE <- as.Date(trimws(precip$DATE))
  precip <- precip[trimws(precip$watershed) == "W3", c("DATE", "Precip")]
  names(precip) <- c("date", "precip_mm")
  
  # ── Daily Streamflow ─────────────────────────────────────────────────────────
  flow <- read.csv(file.path(DATA_DIR, "HBEF_DailyStreamflow_1956-2024.csv"),
                   stringsAsFactors = FALSE)
  flow$DATE <- as.Date(trimws(flow$DATE))
  flow <- flow[flow$WS == 3, c("DATE", "Streamflow")]
  names(flow) <- c("date", "flow_mm")
  
  # ── Daily Air temperature ────────────────────────────────────────────────────
  temp <- read.csv(file.path(DATA_DIR, "HBEF_air_temp_daily.csv"),
                   stringsAsFactors = FALSE)
  temp$date <- as.Date(trimws(temp$date))
  
  # Restrict to the 2 primary HBEF stations
  temp <- temp[temp$STA %in% c("STA14", "STA17"), ]
  
  temp_avg <- temp |>
    group_by(date) |>
    summarise(tavg = mean(AVE, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(tavg) & is.finite(tavg))
  
  # ── Merge & clean ─────────────────────────────────────────────────────────────
  # Inner-join where it keeps only dates that appear in all 3 datasets
  df <- merge(precip, flow,     by = "date", all = FALSE)
  df <- merge(df,     temp_avg, by = "date", all = FALSE)
  df <- df[order(df$date), ]
  df$precip_mm[is.na(df$precip_mm)] <- 0   # treat any leftover P gaps as zero
  df
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  HAMON PET (mm/day)
#
# Estimates potential evapotranspiration (PET) from temperature & day length.
# Returns zero on days where tavg ≤ 0°C.
# ═══════════════════════════════════════════════════════════════════════════════
hamon_pet <- function(tavg, doy, lat_deg = HBEF_LAT) {
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
run_topmodel <- function(pars, data, topidx, spinup_days = 365) {
  
  n    <- nrow(data)
  P    <- data$precip_mm
  tavg <- data$tavg
  obs  <- data$flow_mm
  
  # Pre-compute daily PET for the entire record
  pet  <- hamon_pet(tavg, yday(data$date))
  
  # ── Unpack parameters ────────────────────────────────────────────────────────
  lnTe    <- pars["lnTe"]          # log transmissivity at saturation
  m       <- pars["m"]             # transmissivity decay with depth (mm)
  Srmax   <- pars["Srmax"]         # maximum root zone storage capacity (mm)
  td      <- max(pars["td"], 0.1)  # unsaturated zone time delay (days)
  snow_t  <- pars["snow_t"]        # snow/rain temperature threshold (°C)
  snow_mf <- pars["snow_mf"]       # degree-day melt factor (mm/°C/day)
  
  # ── Build sorted TWI look-up table ────────────────────────────────────────────
  # A cell is saturated when m*(TWI - lambda) >= Sd. Sort thresholds high -> low
  # so we can look up the saturated fraction with a binary search each day.
  twi_vals   <- topidx[, 1]               # TWI class values
  twi_frac   <- topidx[, 2]               # fraction of catchment in each class
  lam        <- sum(twi_vals * twi_frac)  # catchment-average TWI (lambda)
  thr_sorted <- m * (twi_vals - lam)      # deficit threshold for each class
  ord        <- order(thr_sorted, decreasing = TRUE)
  thr_desc   <- thr_sorted[ord]           # thresholds, sorted high to low
  cum_area   <- cumsum(twi_frac[ord])     # running sum of saturated fraction
  n_classes  <- length(thr_desc)
  
  # ── Smart initial saturation deficit from first observed flow ─────────────────
  # Using observed flow avoids a hard-coded Sd0 parameter. If no valid flow
  # record exists, fall back to q0 = 1 mm/d.
  q0 <- obs[which(!is.na(obs) & obs > 0)[1]]
  if (is.na(q0) || !is.finite(q0)) q0 <- 1
  Sd <- -m * (log(max(q0, 0.01)) - lnTe)
  if (Sd > SD_CAP) Sd <- SD_CAP
  if (Sd < 0)      Sd <- 0
  
  # ── Initialize state variables ────────────────────────────────────────────────
  Srz  <- Srmax / 2  # root zone storage deficit
  Suz  <- 0          # unsaturated zone storage (mm)
  SWE  <- 0          # snow water equivalent (mm)
  Te   <- exp(lnTe)  # transmissivity at saturation
  td_m <- td / m     # scaled time delay
  
  # ── Pre-allocate output arrays ────────────────────────────────────────────────
  Qsim    <- numeric(n)
  Qb_out  <- numeric(n)
  Qof_out <- numeric(n)
  SWE_out <- numeric(n)
  ET_out  <- numeric(n)
  Sd_out  <- numeric(n)
  SF_out  <- numeric(n)
  
  # ════════════════════════════════════════════════════════════════════════════
  # Daily time-step loop
  # ════════════════════════════════════════════════════════════════════════════
  for (t in seq_len(n)) {
    
    Tt <- tavg[t]
    
    # ── Snow accumulation & melt (degree-day model) ───────────────────────────
    if (Tt <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0             # all precip stored as snow
    } else {
      melt <- snow_mf * (Tt - snow_t)          # potential snowmelt
      if (melt > SWE) melt <- SWE              # can't melt more than exists
      SWE  <- SWE - melt; rain <- P[t] + melt  # liquid water = rain + melt
    }
    
    SWE_out[t] <- SWE
    
    # ── Saturated area fraction (binary search on pre-sorted thresholds) ──────
    # sf = fraction of watershed where TWI ≥ (lambda + Sd/m)
    if (Sd >= thr_desc[1L]) {
      sat_frac <- 0                    # Sd bigger than every threshold -> nothing saturated
    } else if (Sd < thr_desc[n_classes]) {
      sat_frac <- cum_area[n_classes]  # Sd below every threshold -> fully saturated
    } else {
      lo <- 1L; hi <- n_classes
      while (hi - lo > 1L) {
        mid <- (lo + hi) %/% 2L
        if (thr_desc[mid] > Sd) lo <- mid else hi <- mid
      }
      sat_frac <- cum_area[lo]
    }
    SF_out[t] <- sat_frac
    
    # ── Saturation-excess overland flow ───────────────────────────────────────
    qof   <- sat_frac * rain   # rain falling on saturated cells runs off immediately
    infil <- rain - qof  # remainder infiltrates
    Qof_out[t] <- qof
    
    # ── Root zone water balance ───────────────────────────────────────────────
    Srz <- Srz - infil              # infiltration reduces the root zone deficit
    if (Srz < 0) {
      Suz <- Suz - Srz; Srz <- 0 }  # excess drains to unsaturated zone
    
    # Actual ET is proportional to how full the root zone is
    ae <- if (Srz < Srmax) pet[t] * (1 - Srz / Srmax) else 0
    Srz <- Srz + ae
    if (Srz > Srmax) Srz <- Srmax
    ET_out[t] <- ae
    
    # ── Unsaturated zone drainage (time-delayed by td) ────────────────────────
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd)
      if (quz > Suz) quz <- Suz
      Suz <- Suz - quz
    } else {
      quz <- 0 }
    
    # ── Saturated zone: update deficit then compute exponential baseflow ───────
    Sd <- Sd - quz; if (Sd < 0) Sd <- 0
    
    # Baseflow — falls off exponentially as the deficit grows
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + QB_BUFFER) qb <- Sd + QB_BUFFER
    if (qb < 0)              qb <- 0
    Sd <- Sd + qb
    if (Sd < 0) Sd <- 0
    
    # Total simulated streamflow this day
    Qb_out[t]  <- qb
    Qsim[t]    <- qb + qof  # total discharge = baseflow + overland flow
    Sd_out[t]  <- Sd
  }
  
  # Remove any non-finite simulated values (numerical edge cases)
  Qsim[!is.finite(Qsim)] <- 0
  Qsim <- pmax(Qsim, 0)
  
  # ── Nash-Sutcliffe Efficiency (full record, post spin-up) ─────────────────
  # Compute NSE after the spin-up period so the initial-state guesses don't
  # pollute the score.
  warmup <- min(spinup_days, n - 1)
  idx    <- (warmup + 1):n
  obs_w  <- obs[idx]
  sim_w  <- Qsim[idx]
  valid  <- !is.na(obs_w) & is.finite(sim_w)
  
  if (sum(valid) < MIN_NSE_N) {
    nse <- NA_real_
  } else {
    ss_res <- sum((obs_w[valid] - sim_w[valid])^2)
    ss_tot <- sum((obs_w[valid] - mean(obs_w[valid]))^2)
    nse    <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  }
  
  list(
    Qsim = Qsim, Qobs = obs, Qb = Qb_out, Qof = Qof_out,
    SWE = SWE_out, ET = ET_out, Sd = Sd_out, SF = SF_out,
    dates = data$date, nse = nse, warmup = warmup, m = m, lam = lam
  )
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
  
  lnTe    <- pv[1]
  m       <- pv[2]
  Srmax   <- pv[3]
  td      <- pv[4]
  snow_t  <- pv[5]
  snow_mf <- pv[6]
  if (td < 0.1) td <- 0.1
  
  # TWI lookup (same as run_topmodel).
  twi_vals   <- topidx[, 1]
  twi_frac   <- topidx[, 2]
  lam        <- sum(twi_vals * twi_frac)
  thr_sorted <- m * (twi_vals - lam)
  ord        <- order(thr_sorted, decreasing = TRUE)
  thr_desc   <- thr_sorted[ord]
  cum_area   <- cumsum(twi_frac[ord])
  n_classes  <- length(thr_desc)
  
  # Smart Sd0 from first valid observed flow
  q0 <- obs[which(!is.na(obs) & obs > 0)[1]]
  if (is.na(q0) || !is.finite(q0)) q0 <- 1
  Sd <- -m * (log(max(q0, 0.01)) - lnTe)
  if (Sd > SD_CAP) Sd <- SD_CAP
  if (Sd < 0)      Sd <- 0
  
  Srz <- Srmax / 2
  Suz <- 0
  SWE <- 0
  Te  <- exp(lnTe)
  td_m <- td / m
  
  # Pre-check that there are enough post-warmup observations to score
  obs_after <- obs[(warmup + 1):n]
  obs_valid <- obs_after[!is.na(obs_after)]
  if (length(obs_valid) < MIN_NSE_N) return(1e6)
  obs_mean <- mean(obs_valid)  # observed mean for NSE denominator
  
  # Accumulate residuals on-the-fly without storing arrays
  # Running totals for NSE = 1 - ss_res/ss_tot.
  ss_res <- 0
  ss_tot <- 0
  n_used <- 0L
  
  for (t in seq_len(n)) {
    
    Tt <- tavg[t]
    
    # Snow
    if (Tt <= snow_t) {
      SWE  <- SWE + P[t]
      rain <- 0
    } else {
      melt <- snow_mf * (Tt - snow_t)
      if (melt > SWE) melt <- SWE
      SWE  <- SWE - melt
      rain <- P[t] + melt
    }
    
    # Saturated fraction (binary search)
    if (Sd >= thr_desc[1L]) {
      sat_frac <- 0
    } else if (Sd < thr_desc[n_classes]) {
      sat_frac <- cum_area[n_classes]
    } else {
      lo2 <- 1L; hi <- n_classes
      while (hi - lo2 > 1L) {
        mid <- (lo2 + hi) %/% 2L
        if (thr_desc[mid] > Sd) lo2 <- mid else hi <- mid
      }
      sat_frac <- cum_area[lo2]
    }
    
    # Infiltration, root zone, ET
    infil <- (1 - sat_frac) * rain
    Srz   <- Srz - infil
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }
    
    ae <- if (Srz < Srmax) pet[t] * (1 - Srz / Srmax) else 0
    Srz <- Srz + ae
    if (Srz > Srmax) Srz <- Srmax
    
    # Recharge, baseflow
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd)
      if (quz > Suz) quz <- Suz
      Suz <- Suz - quz
    } else {
      quz <- 0
    }
    Sd <- Sd - quz
    if (Sd < 0) Sd <- 0
    
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + QB_BUFFER) qb <- Sd + QB_BUFFER
    if (qb < 0)              qb <- 0
    Sd <- Sd + qb
    if (Sd < 0) Sd <- 0
    
    qsim <- qb + sat_frac * rain
    
    # Only score days past the spin-up.
    if (t > warmup) {
      ot <- obs[t]
      if (!is.na(ot) && is.finite(qsim)) {
        ss_res <- ss_res + (ot - qsim)^2
        ss_tot <- ss_tot + (ot - obs_mean)^2
        n_used <- n_used + 1L
      }
    }
  }
  
  if (n_used < MIN_NSE_N || ss_tot <= 0) return(1e6)
  -(1 - ss_res / ss_tot)  # minimize -NSE
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  OPTIMIZER — Nelder-Mead on a user-defined calibration window
#
# Runs optim() with method = "Nelder-Mead" from three different starting
# points to reduce the chance of converging to a local optimum, then returns
# the best (highest NSE) parameter set clamped to valid bounds.
# ═══════════════════════════════════════════════════════════════════════════════
optimize_params <- function(data, topidx, cal_start, cal_end,
                            spinup_days = 365, progress_fn = NULL) {
  
  # Bounds come from PARAM_CONFIG so they match the sliders.
  lo <- PARAM_LO
  up <- PARAM_UP
  
  # Subset data to calibration window; fall back to full record if too short
  cal <- data[data$date >= cal_start & data$date <= cal_end, ]
  if (nrow(cal) < spinup_days + MIN_NSE_N) cal <- data
  
  # Pre-compute the forcings once so each optim call can reuse them.
  cal_n      <- nrow(cal)
  cal_P      <- cal$precip_mm
  cal_tavg   <- cal$tavg
  cal_obs    <- cal$flow_mm
  cal_pet    <- hamon_pet(cal_tavg, yday(cal$date))
  cal_warmup <- min(spinup_days, cal_n - 1)
  
  # Wrap the NSE function so optim gets back a plain scalar and can't crash.
  obj <- function(pv) {
    tryCatch(
      topmodel_nse(pv, cal_P, cal_tavg, cal_obs, cal_pet,
                   topidx, cal_n, cal_warmup, lo, up),
      error = function(e) 1e6
    )
  }
  
  # 3 starting points covering different parts of parameter space
  starts <- list(
    c(lnTe =  2, m = 23, Srmax =  68, td =  1, snow_t = -0.7, snow_mf = 2.9),
    c(lnTe =  1, m = 30, Srmax = 100, td = 10, snow_t =  0.0, snow_mf = 3.0),
    c(lnTe = -2, m = 40, Srmax = 150, td = 20, snow_t = -0.5, snow_mf = 3.5)
  )
  
  best_val <- Inf
  best_par <- starts[[1]]
  
  # Try each start, keep the best final value.
  for (i in seq_along(starts)) {
    
    # Notify the UI of progress if a callback is provided
    if (!is.null(progress_fn)) {
      progress_fn(i, length(starts),
                  sprintf("Start %d of %d \u2026", i, length(starts)))
    }
    opt <- tryCatch(
      optim(starts[[i]], obj,
            method  = "Nelder-Mead",
            control = list(maxit = 400, reltol = 1e-6)),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value) && opt$value < best_val) {
      best_val <- opt$value
      best_par <- opt$par
    }
  }
  
  # Clamp to bounds as a final safety net.
  names(best_par) <- PARAM_NAMES
  pmax(pmin(best_par, up), lo)
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: compute NSE over an arbitrary date window
#
# Used to display calibration-period and viewing-window NSE separately in the
# plot title, independently of the full-record NSE stored in the result list.
# ═══════════════════════════════════════════════════════════════════════════════
compute_nse <- function(res, start_date, end_date) {
  
  # Only include post-warmup days that fall within the requested window
  vi <- which(res$dates >= start_date &
                res$dates <= end_date  &
                seq_along(res$Qsim) > res$warmup)
  if (length(vi) < MIN_NSE_N) return(NA_real_)
  
  o  <- res$Qobs[vi]
  s  <- res$Qsim[vi]
  ok <- !is.na(o) & is.finite(s)
  if (sum(ok) < MIN_NSE_N) return(NA_real_)
  
  ss_res <- sum((o[ok] - s[ok])^2)
  ss_tot <- sum((o[ok] - mean(o[ok]))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - ss_res / ss_tot
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
# Load precomputed spatial data 
# TWI, DEM, and watershed boundary are built offline by a WhiteBox script
# and saved as .rds. That way the deployed Shiny app doesn't need WhiteBox.
message("Loading data ...")

ALL_DATA <- load_hbef_data()
message(sprintf("  %d records: %s to %s",
                nrow(ALL_DATA), min(ALL_DATA$date), max(ALL_DATA$date)))

TWI_MAT <- readRDS(file.path(DATA_DIR, "twi_matrix.rds"))    # TWI raster as matrix
TWI_EXT <- readRDS(file.path(DATA_DIR, "twi_extent.rds"))    # xmin/xmax/ymin/ymax
DEM_MAT <- readRDS(file.path(DATA_DIR, "dem_matrix.rds"))    # DEM raster as matrix
DEM_EXT <- readRDS(file.path(DATA_DIR, "dem_extent.rds"))
WS_POLY <- readRDS(file.path(DATA_DIR, "ws_boundary.rds"))   # watershed outline (x,y)
TOPIDX  <- readRDS(file.path(DATA_DIR, "topidx.rds"))        # TWI class/fraction table
TWI_LAM <- readRDS(file.path(DATA_DIR, "twi_lam.rds"))       # catchment mean TWI

# Build axis vectors, and flip+transpose the matrices into the orientation
# image() expects (x along cols, y along rows, going bottom-to-top).
twi_x   <- seq(TWI_EXT["xmin"], TWI_EXT["xmax"], length.out = ncol(TWI_MAT))
twi_y   <- seq(TWI_EXT["ymin"], TWI_EXT["ymax"], length.out = nrow(TWI_MAT))
TWI_IMG <- t(TWI_MAT[nrow(TWI_MAT):1, ])

dem_x   <- seq(DEM_EXT["xmin"], DEM_EXT["xmax"], length.out = ncol(DEM_MAT))
dem_y   <- seq(DEM_EXT["ymin"], DEM_EXT["ymax"], length.out = nrow(DEM_MAT))
DEM_IMG <- t(DEM_MAT[nrow(DEM_MAT):1, ])

# Flatten TWI for the histogram tab.
TWI_VALS  <- as.vector(TWI_MAT)
TWI_VALID <- !is.na(TWI_VALS)

# Year range used by the slider widgets.
YEAR_MIN <- as.integer(format(min(ALL_DATA$date), "%Y"))
YEAR_MAX <- as.integer(format(max(ALL_DATA$date), "%Y"))

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
  title = "HBEF Watershed 3 \u2014 TOPMODEL Explorer",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  tags$head(tags$style(HTML(css))),
  
  # ---- Sidebar: controls ----
  sidebar = sidebar(
    width = 340,
    
    # Viewing window (just changes the plot x-axis, not the model run).
    sliderInput("yr", "Viewing Period",
                min = YEAR_MIN, max = YEAR_MAX,
                value = c(YEAR_MIN, YEAR_MAX),
                step = 1, sep = "", ticks = TRUE),
    helpText("Adjust the time window shown on the plots. The model always runs on the full dataset."),
    
    # Action buttons.
    actionButton("run", "\u25B6  Run Model",
                 class = "btn-success action-btn"),
    downloadButton("download_csv", "\u2B07  Export Output (CSV)",
                   class = "btn-secondary action-btn"),
    uiOutput("opt_disp"),    # text box showing the latest optimized values
    
    # Subsurface flow parameter group
    div(class = "pg",
        h6("\u2193 Subsurface Flow Parameters"),
        helpText(PARAM_CONFIG$lnTe$help),
        sliderInput("lnTe", PARAM_CONFIG$lnTe$label,
                    PARAM_CONFIG$lnTe$lo, PARAM_CONFIG$lnTe$up,
                    PARAM_CONFIG$lnTe$default, PARAM_CONFIG$lnTe$step),
        helpText(PARAM_CONFIG$m$help),
        sliderInput("m", PARAM_CONFIG$m$label,
                    PARAM_CONFIG$m$lo, PARAM_CONFIG$m$up,
                    PARAM_CONFIG$m$default, PARAM_CONFIG$m$step),
        helpText(PARAM_CONFIG$td$help),
        sliderInput("td", PARAM_CONFIG$td$label,
                    PARAM_CONFIG$td$lo, PARAM_CONFIG$td$up,
                    PARAM_CONFIG$td$default, PARAM_CONFIG$td$step)
    ),
    
    # Root zone & ET group
    div(class = "pg",
        h6("\U0001F331 Root Zone & Evapotranspiration"),
        helpText(PARAM_CONFIG$Srmax$help),
        sliderInput("Srmax", PARAM_CONFIG$Srmax$label,
                    PARAM_CONFIG$Srmax$lo, PARAM_CONFIG$Srmax$up,
                    PARAM_CONFIG$Srmax$default, PARAM_CONFIG$Srmax$step)
    ),
    
    # Snow group
    div(class = "pg",
        h6("\u2744 Snow Accumulation & Melt"),
        helpText(PARAM_CONFIG$snow_t$help),
        sliderInput("snow_t", PARAM_CONFIG$snow_t$label,
                    PARAM_CONFIG$snow_t$lo, PARAM_CONFIG$snow_t$up,
                    PARAM_CONFIG$snow_t$default, PARAM_CONFIG$snow_t$step),
        helpText(PARAM_CONFIG$snow_mf$help),
        sliderInput("snow_mf", PARAM_CONFIG$snow_mf$label,
                    PARAM_CONFIG$snow_mf$lo, PARAM_CONFIG$snow_mf$up,
                    PARAM_CONFIG$snow_mf$default, PARAM_CONFIG$snow_mf$step)
    ),
    
    # Spin-up + calibration window controls
    div(class = "pg",
        h6("\u23F1 Spin-Up & Calibration"),
        helpText("Spin-up days: model runs but NSE is not computed during this period (lets storages stabilize from initial guesses)."),
        sliderInput("spinup", "Spin-up (days)", 90, 1095, 365, 30),
        helpText("Calibration period (used only by the Optimize button)."),
        sliderInput("cal_yr", "Calibration Years",
                    min = YEAR_MIN, max = YEAR_MAX,
                    value = c(YEAR_MIN, YEAR_MAX),
                    step = 1, sep = "", ticks = TRUE)
    ),
    
    actionButton("opt", "\u26A1 Optimize (Nelder-Mead)",
                 class = "btn-primary action-btn"),
    helpText("NOTE: Nelder-Mead is a local optimizer — results depend on starting parameter values. Different seeds may converge to different local optima. This optimizer tries 3 different starting points and keeps the best one.")
  ),
  
  #Main panel: one tab per view
  navset_card_tab(
    nav_panel("Model Output",                    plotlyOutput("model_out",  height = "950px")),
    nav_panel("Input Data",                      plotlyOutput("input_plot", height = "640px")),
    nav_panel("Observed vs Simulated Discharge", plotlyOutput("scatter",    height = "520px")),
    nav_panel("Flow Duration",                   plotlyOutput("fdc",        height = "520px")),
    nav_panel("TWI Map",                         plotOutput ("twi_map",     height = "560px")),
    
    nav_panel("TWI Histogram",
              fluidRow(
                column(4,
                       sliderInput("twi_classes", "Number of TWI classes", 4, 40, 16, 1),
                       helpText("TWI classes discretize the topography for TOPMODEL. More classes = finer detail."),
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
              ))
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  SHINY SERVER
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  rv  <- reactiveVal(NULL)  # stores the most recent run_topmodel() result
  otx <- reactiveVal(NULL)  # stores formatted optimized parameter string
  
  # Current viewing window as Date objects.
  view_range <- reactive({
    req(input$yr)
    list(start = as.Date(paste0(input$yr[1], "-01-01")),
         end   = as.Date(paste0(input$yr[2], "-12-31")))
  })
  
  # Current calibration window as Date objects.
  cal_range <- reactive({
    list(start = as.Date(paste0(input$cal_yr[1], "-01-01")),
         end   = as.Date(paste0(input$cal_yr[2], "-12-31")))
  })
  
  # Collect slider values into a named parameter vector.
  cpars <- reactive({
    c(lnTe    = input$lnTe,
      m       = input$m,
      Srmax   = input$Srmax,
      td      = input$td,
      snow_t  = input$snow_t,
      snow_mf = input$snow_mf)
  })
  
  # Run the model once on startup so the tabs aren't empty.
  once <- reactiveVal(FALSE)
  observe({
    req(!once())
    res <- tryCatch(
      run_topmodel(isolate(cpars()), ALL_DATA, TOPIDX,
                   spinup_days = isolate(input$spinup)),
      error = function(e) NULL
    )
    if (!is.null(res)) rv(res)
    once(TRUE)
  })
  
  # Manual "Run Model" button: rerun with the current slider values.
  observeEvent(input$run, {
    res <- tryCatch(
      run_topmodel(cpars(), ALL_DATA, TOPIDX, spinup_days = input$spinup),
      error = function(e) {
        showNotification(e$message, type = "error")
        NULL
      }
    )
    if (!is.null(res)) {
      rv(res)
      cr <- cal_range()
      nc <- compute_nse(res, cr$start, cr$end)
      showNotification(
        sprintf("NSE (cal) = %s",
                if (is.finite(nc)) sprintf("%.3f", nc) else "N/A"),
        type = "message", duration = 3
      )
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)
  
  # "Optimize" button: runs the 3-start Nelder-Mead and updates the sliders.
  observeEvent(input$opt, {
    showNotification("Optimizing \u2026 Start 1 of 3",
                     type = "message", duration = NULL, id = "op")
    
    cr <- cal_range()
    best <- optimize_params(
      ALL_DATA, TOPIDX, cr$start, cr$end,
      spinup_days = input$spinup,
      progress_fn = function(i, n, msg) {
        # Update the status notification as each start finishes.
        showNotification(sprintf("Optimizing \u2026 %s", msg),
                         type = "message", duration = NULL, id = "op")
      }
    )
    removeNotification("op")
    
    # Push the optimized values back onto the sliders (snapped to the grid).
    for (p in PARAM_NAMES) {
      updateSliderInput(session, p,
                        value = snap(best[p],
                                     PARAM_CONFIG[[p]]$lo,
                                     PARAM_CONFIG[[p]]$up,
                                     PARAM_CONFIG[[p]]$step))
    }
    
    # Show the optimized values as text in the sidebar.
    otx(sprintf("lnTe=%.2f  m=%.1f  Srmax=%.1f\ntd=%.1f  Tmelt=%.2f  Melt=%.2f",
                best["lnTe"], best["m"], best["Srmax"],
                best["td"],   best["snow_t"], best["snow_mf"]))
    
    # Rerun the model with the optimized values so the plots update.
    res <- tryCatch(
      run_topmodel(best, ALL_DATA, TOPIDX, spinup_days = input$spinup),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      rv(res)
      nc <- compute_nse(res, cr$start, cr$end)
      showNotification(
        sprintf("Optimized! NSE (cal) = %s",
                if (is.finite(nc)) sprintf("%.3f", nc) else "N/A"),
        type = "message", duration = 5
      )
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Renders the optimized-values box in the sidebar (blank until Optimize runs).
  output$opt_disp <- renderUI({
    txt <- otx()
    if (is.null(txt)) return(NULL)
    div(class = "opt-params",
        tags$strong("Optimized values:"), tags$br(),
        tags$pre(style = "margin:2px 0 0 0;white-space:pre-wrap;", txt))
  })
  
  # CSV download: assembles the full output table and writes it.
  output$download_csv <- downloadHandler(
    filename = function() sprintf("topmodel_output_%s.csv",
                                  format(Sys.Date(), "%Y%m%d")),
    content  = function(file) {
      res <- rv(); req(res)
      d <- data.frame(
        date               = res$dates,
        P_mm               = ALL_DATA$precip_mm,
        T_C                = ALL_DATA$tavg,
        Qobs_mm            = res$Qobs,
        Qsim_mm            = res$Qsim,
        Qbaseflow_mm       = res$Qb,
        Qsat_excess_mm     = res$Qof,
        ET_mm              = res$ET,
        SWE_mm             = res$SWE,
        Storage_Deficit_mm = res$Sd,
        Saturated_Fraction = res$SF
      )
      write.csv(d, file, row.names = FALSE)
    }
  )
  
  #Model Output tab: 7 stacked time-series panels
  output$model_out <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    
    # Wide-format table for plotly.
    d <- data.frame(
      date = res$dates,
      Qobs = res$Qobs, Qsim = res$Qsim,
      Qof  = res$Qof,  Qb   = res$Qb,
      ET   = res$ET,   SWE  = res$SWE,
      Sd   = res$Sd,   SF   = res$SF * 100
    )
    
    # Compute two NSE values: one over the calibration window, one over
    # whatever the user is currently looking at.
    vr <- view_range()
    cr <- cal_range()
    nse_cal_val  <- compute_nse(res, cr$start, cr$end)
    nse_view_val <- compute_nse(res, vr$start, vr$end)
    nse_cal  <- if (is.finite(nse_cal_val))  sprintf("%.3f", nse_cal_val)  else "N/A"
    nse_view <- if (is.finite(nse_view_val)) sprintf("%.3f", nse_view_val) else "N/A"
    
    yax <- list(title = "", nticks = 4, automargin = TRUE)
    
    # Build each of the 7 subplots.
    p1 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~Qobs, name = "Observed",  line = list(color = "#2980b9", width = 1.2)) |> 
      add_lines(y = ~Qsim, name = "Simulated", line = list(color = "#e74c3c", width = 1.2)) |> 
      layout(yaxis = yax)
    
    p2 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qof, name = "SatExc", line = list(color = "#e67e22", width = 1)) |>
      layout(yaxis = yax)
    
    p3 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qb, name = "Baseflow", line = list(color = "#8e44ad", width = 1)) |>
      layout(yaxis = yax)
    
    p4 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~ET, name = "ET", line = list(color = "#27ae60", width = 1)) |>
      layout(yaxis = yax)
    
    p5 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~SWE, name = "SWE", line = list(color = "#3498db", width = 1)) |>
      layout(yaxis = yax)
    
    p6 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~Sd, name = "Deficit", line = list(color = "#c0392b", width = 1)) |>
      layout(yaxis = yax)
    
    p7 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~SF, name = "Sat%", line = list(color = "#2166ac", width = 1)) |>
      layout(yaxis = yax, xaxis = list(title = ""))
    
    # Relative heights of each panel (must sum to 1).
    panel_h <- c(.22, .11, .11, .11, .15, .15, .15)
    labels  <- c(
      "Discharge (mm/d)",
      "Surface Runoff as Saturation-Excess Overland Flow (mm/d)",
      "Baseflow (mm/d)",
      "ET (mm/d)",
      "SWE (mm)",
      "Storage Deficit (mm)",
      "Saturated Area (%)"
    )
    # Figure out where the top of each panel sits (paper coords) so we can
    # drop a label annotation at that y position.
    cum_h <- cumsum(panel_h)
    tops  <- 1 - c(0, cum_h[-length(cum_h)])
    
    ann <- lapply(seq_along(labels), function(i) {
      list(text      = paste0("<b>", labels[i], "</b>"),
           x = 0.01, y = tops[i] - 0.005,
           xref = "paper", yref = "paper",
           xanchor = "left", yanchor = "top",
           showarrow = FALSE,
           font      = list(size = 11, color = "#333"),
           bgcolor   = "rgba(255,255,255,0.8)",
           borderpad = 2)
    })
    
    # Stack all 7 with shared x-axis and zoom into the viewing window.
    subplot(p1, p2, p3, p4, p5, p6, p7,
            nrows = 7, shareX = TRUE, titleY = FALSE, heights = panel_h) |>
      layout(title = list(
        text = sprintf("Model Output  |  NSE (cal) = %s  |  NSE (view) = %s",
                       nse_cal, nse_view),
        font = list(size = 16)),
        showlegend  = FALSE,
        hovermode   = "x unified",
        annotations = ann,
        xaxis       = list(range = as.character(c(vr$start, vr$end))),
        margin      = list(t = 50, b = 30, l = 70))
  })
  
  #Input Data tab: P, T, Qobs (model-independent)
  output$input_plot <- renderPlotly({
    d  <- ALL_DATA
    vr <- view_range()
    
    p1 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~precip_mm, line = list(color = "#2980b9", width = 1)) |>
      layout(yaxis = list(title = "Precip (mm/d)", nticks = 4, automargin = TRUE))
    
    p2 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~tavg, line = list(color = "#e67e22", width = 1)) |>
      layout(yaxis = list(title = "Air Temp (\u00B0C)", nticks = 4, automargin = TRUE))
    
    p3 <- plot_ly(d, x = ~date) |>
      add_lines(y = ~flow_mm, line = list(color = "#27ae60", width = 1)) |>
      layout(yaxis = list(title = "Obs Q (mm/d)", nticks = 4, automargin = TRUE))
    
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE) |>
      layout(title = list(
        text = "Input Data: Precipitation, Temperature, Observed Discharge",
        font = list(size = 15)),
        showlegend = FALSE,
        hovermode  = "x unified",
        xaxis      = list(range = as.character(c(vr$start, vr$end))),
        margin     = list(t = 50, l = 90))
  })
  
  #Scatter tab: observed vs simulated Q
  output$scatter <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model."))
    
    # Only use days in the viewing window and past the spin-up.
    vr <- view_range()
    vi <- which(res$dates >= vr$start &
                  res$dates <= vr$end  &
                  seq_along(res$Qsim) > res$warmup)
    
    d <- data.frame(obs = res$Qobs[vi], sim = res$Qsim[vi])
    d <- d[!is.na(d$obs) & is.finite(d$sim), ]
    validate(need(nrow(d) > 0, "No data."))
    
    # 1:1 reference line spans the combined data range.
    rng <- range(c(d$obs, d$sim), na.rm = TRUE)
    
    plot_ly(d, x = ~obs, y = ~sim,
            type = "scattergl", mode = "markers",    # scattergl = WebGL for speed
            marker = list(size = 3, color = "#3498db", opacity = .4)) |>
      add_lines(x = rng, y = rng, name = "1:1",
                line = list(color = "black", dash = "dash", width = 1)) |>
      layout(title = list(
        text = "Observed vs Simulated Discharge (viewing window)",
        font = list(size = 15)),
        xaxis      = list(title = "Observed (mm/d)"),
        yaxis      = list(title = "Simulated (mm/d)"),
        showlegend = FALSE,
        margin     = list(t = 50))
  })
  
  #Flow-Duration Curve tab
  output$fdc <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), "Press \u25B6 Run Model."))
    
    vr <- view_range()
    vi <- which(res$dates >= vr$start &
                  res$dates <= vr$end  &
                  seq_along(res$Qsim) > res$warmup)
    
    # Sort high -> low and compute exceedance probability on the x-axis.
    obs_sorted <- sort(res$Qobs[vi][!is.na(res$Qobs[vi])],    decreasing = TRUE)
    sim_sorted <- sort(res$Qsim[vi][is.finite(res$Qsim[vi])], decreasing = TRUE)
    exc_obs <- seq_along(obs_sorted) / length(obs_sorted) * 100
    exc_sim <- seq_along(sim_sorted) / length(sim_sorted) * 100
    
    plot_ly() |>
      add_lines(x = exc_obs, y = obs_sorted, name = "Observed",
                line = list(color = "#2980b9", width = 1.5)) |>
      add_lines(x = exc_sim, y = sim_sorted, name = "Simulated",
                line = list(color = "#e74c3c", width = 1.5)) |>
      layout(title = list(
        text = "Flow-Duration Curve (viewing window)",
        font = list(size = 15)),
        xaxis  = list(title = "Exceedance (%)"),
        yaxis  = list(title = "Discharge (mm/d)", type = "log"),
        legend = list(orientation = "h", y = -0.12),
        margin = list(t = 50))
  })
  
  # TWI Map tab: base-R image + watershed outline + DEM contours
  output$twi_map <- renderPlot({
    par(mar = c(2, 2, 3, 4), bg = "white")
    image(twi_x, twi_y, TWI_IMG,
          col  = hcl.colors(50, "YlGnBu", rev = TRUE),
          asp  = 1, axes = FALSE, xlab = "", ylab = "",
          main = sprintf("Topographic Wetness Index  (FD8, \u03BB = %.1f)", TWI_LAM),
          cex.main = 1.3)
    polygon(WS_POLY[, 1], WS_POLY[, 2],
            border = "#222", lwd = 2.5, col = NA)           # watershed boundary
    contour(dem_x, dem_y, DEM_IMG, add = TRUE, nlevels = 10,
            col = adjustcolor("#333", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)     # elevation contours
    mtext("Contour elevations in meters",
          side = 1, line = 0.5, cex = 0.8, col = "#666")
  }, res = 120)
  
  #TWI Histogram tab
  output$twi_hist <- renderPlot({
    vals <- TWI_VALS[TWI_VALID]
    par(mar = c(4, 4, 3, 1), bg = "white")
    hist(vals, breaks = input$twi_classes,
         col = "#3498db", border = "white",
         main = sprintf("Distribution of TWI values (%d classes)",
                        input$twi_classes),
         xlab = "TWI", ylab = "Cell count", cex.main = 1.3)
    # Vertical reference lines at mean and median.
    abline(v = mean(vals),   col = "#e74c3c", lwd = 2, lty = 2)
    abline(v = median(vals), col = "#27ae60", lwd = 2, lty = 2)
    legend("topright",
           legend = c(sprintf("Mean = %.2f",   mean(vals)),
                      sprintf("Median = %.2f", median(vals))),
           col = c("#e74c3c", "#27ae60"),
           lwd = 2, lty = 2, bty = "n")
  }, res = 110)
  
  # Text summary for the TWI histogram panel.
  output$twi_stats <- renderText({
    vals <- TWI_VALS[TWI_VALID]
    sprintf(
      "TWI Summary\n  Min:    %.2f\n  Max:    %.2f\n  Mean:   %.2f\n  Median: %.2f\n  SD:     %.2f\n  Cells:  %d\n  Lambda: %.2f",
      min(vals), max(vals), mean(vals), median(vals),
      sd(vals), length(vals), TWI_LAM)
  })
  
  # Indices of every day inside the current viewing window.
  view_idx <- reactive({
    res <- rv(); req(!is.null(res))
    vr  <- view_range()
    which(res$dates >= vr$start & res$dates <= vr$end)
  })
  
  # When the viewing window changes, reset the day slider.
  observe({
    vi <- view_idx()
    req(length(vi) > 0)
    updateSliderInput(session, "sat_day",
                      min = 1, max = length(vi), value = 1)
  })
  
  # Click on the discharge plot -> jump sat_day to that date.
  observeEvent(event_data("plotly_click", source = "sat_q_src"), {
    cd <- event_data("plotly_click", source = "sat_q_src")
    if (is.null(cd)) return()
    clicked_date <- as.Date(cd$x)
    vi  <- view_idx(); req(length(vi) > 0)
    res <- rv()
    # Find the index of the day nearest the clicked date.
    pos <- which.min(abs(as.numeric(res$dates[vi] - clicked_date)))
    updateSliderInput(session, "sat_day", value = pos)
  })
  
  # The map itself: shade each cell saturated vs unsaturated.
  output$sat_map <- renderPlot({
    res <- rv()
    validate(need(!is.null(res), "Run model first."))
    
    vi <- view_idx(); req(length(vi) > 0)
    si <- min(input$sat_day %||% 1, length(vi))
    di <- vi[si]
    
    # A cell is saturated when its TWI is >= lambda + Sd/m.
    Sd     <- res$Sd[di]
    thresh <- res$lam + Sd / res$m
    label  <- format(res$dates[di], "%Y-%m-%d")
    
    # Build a 0/1 saturation map with NAs preserved.
    sat_v <- rep(NA_real_, length(TWI_VALS))
    sat_v[TWI_VALID] <- ifelse(TWI_VALS[TWI_VALID] >= thresh, 1, 0)
    sat_mat <- matrix(sat_v, nrow = nrow(TWI_MAT), ncol = ncol(TWI_MAT))
    sat_img <- t(sat_mat[nrow(sat_mat):1, ])
    
    # Percent of the watershed that's saturated today (for the title).
    pct <- round(100 * sum(sat_v == 1, na.rm = TRUE) / sum(TWI_VALID), 1)
    
    par(mar = c(2, 2, 3, 1), bg = "white")
    image(twi_x, twi_y, sat_img,
          col = c("#f5e6c8", "#2166ac"),
          breaks = c(-0.5, 0.5, 1.5),
          asp = 1, axes = FALSE, xlab = "", ylab = "",
          main = sprintf("%s \u2014 %.1f%% saturated", label, pct),
          cex.main = 1.2)
    polygon(WS_POLY[, 1], WS_POLY[, 2],
            border = "#222", lwd = 2.5, col = NA)
    contour(dem_x, dem_y, DEM_IMG, add = TRUE, nlevels = 10,
            col = adjustcolor("#555", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
    legend("bottomright",
           legend = c("Saturated", "Unsaturated"),
           fill   = c("#2166ac",   "#f5e6c8"),
           border = c("#2166ac",   "#f5e6c8"),
           bty    = "n", cex = 0.9)
    mtext("Contour elevations in meters",
          side = 1, line = 0.5, cex = 0.7, col = "#666")
  }, res = 110)
  
  # Discharge plot next to the map — also the source of click events.
  output$sat_q <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), ""))
    
    vi <- view_idx(); req(length(vi) > 0)
    si <- min(input$sat_day %||% 1, length(vi))
    di <- vi[si]
    vr <- view_range()
    
    d <- data.frame(date = res$dates, Qobs = res$Qobs, Qsim = res$Qsim)
    
    p <- plot_ly(d, x = ~date, source = "sat_q_src") |>
      add_lines(y = ~Qobs, name = "Obs",
                line = list(color = "#2980b9", width = 1),
                hovertemplate = "%{x|%Y-%m-%d}<br>Obs: %{y:.2f} mm/d<extra></extra>") |>
      add_lines(y = ~Qsim, name = "Sim",
                line = list(color = "#e74c3c", width = 1),
                hovertemplate = "%{x|%Y-%m-%d}<br>Sim: %{y:.2f} mm/d<extra></extra>") |>
      # Black dot marks the currently selected day.
      add_markers(x = res$dates[di], y = res$Qsim[di], name = "Now",
                  marker = list(color = "black", size = 10),
                  hovertemplate = paste0("<b>",
                                         format(res$dates[di], "%Y-%m-%d"),
                                         "</b><extra></extra>")) |>
      layout(title = list(text = "Discharge (mm/d) — click to jump map",
                          font = list(size = 14)),
             yaxis = list(title = ""),
             xaxis = list(title = "",
                          range = as.character(c(vr$start, vr$end))),
             showlegend = FALSE,
             margin     = list(t = 40, b = 20, l = 40, r = 10),
             hovermode  = "closest")
    
    # event_register has to be the last call — otherwise plotly_click events
    # won't fire for this source.
    event_register(p, "plotly_click")
  })
  
  # Saturated-area percentage timeseries next to the map.
  output$sat_pct <- renderPlotly({
    res <- rv()
    validate(need(!is.null(res), ""))
    
    vi <- view_idx(); req(length(vi) > 0)
    si <- min(input$sat_day %||% 1, length(vi))
    di <- vi[si]
    vr <- view_range()
    
    d <- data.frame(date = res$dates, SF = res$SF * 100)
    
    plot_ly(d, x = ~date) |>
      add_lines(y = ~SF, name = "Sat%",
                line = list(color = "#2166ac", width = 1),
                hovertemplate = "%{x|%Y-%m-%d}<br>Sat: %{y:.1f}%<extra></extra>") |>
      add_markers(x = res$dates[di], y = d$SF[di], name = "Now",
                  marker = list(color = "black", size = 10),
                  hovertemplate = paste0("<b>",
                                         format(res$dates[di], "%Y-%m-%d"),
                                         "</b><extra></extra>")) |>
      layout(title = list(text = "Saturated Catchment (%)",
                          font = list(size = 14)),
             yaxis = list(title = ""),
             xaxis = list(title = "",
                          range = as.character(c(vr$start, vr$end))),
             showlegend = FALSE,
             margin     = list(t = 40, b = 20, l = 40, r = 10),
             hovermode  = "closest")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  LAUNCH
# ═══════════════════════════════════════════════════════════════════════════════
shinyApp(ui, server)