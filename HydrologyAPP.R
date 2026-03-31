# ═══════════════════════════════════════════════════════════════════════════════
# app.R — HBEF Watershed 3 TOPMODEL Explorer
#
# A Shiny app for running TOPMODEL, a rainfall-runoff model that predicts
# streamflow based on watershed topography, soil properties, and climate inputs.
# Built for the Hubbard Brook Experimental Forest (HBEF) Watershed 3.
#
# Samuel Handel, Kaelyn Harvey, Max Hughes
# ═══════════════════════════════════════════════════════════════════════════════

# ── Core libraries ─────────────────────────────────────────────────────────────
library(shiny)      # web app framework
library(bslib)      # bootstrap-based UI themes
library(terra)      # raster/vector spatial data
library(raster)     # package for contour plotting
library(plotly)     # interactive plots
library(dplyr)      # data wrangling
library(lubridate)  # data handling
library(sf)         # vector/shapefile data (watershed boundary)

# WBT is only needed to compute TWI
if (requireNamespace("whitebox", quietly = TRUE)) {
  library(whitebox)
  wbt_init()
  HAS_WHITEBOX <- TRUE
} else {
  HAS_WHITEBOX <- FALSE
}

# ═══════════════════════════════════════════════════════════════════════════════
# 0.  FILE PATHS
# ═══════════════════════════════════════════════════════════════════════════════
DATA_DIR    <- "data"
TWI_F       <- file.path(DATA_DIR, "twi_precomputed.tif")  # TWI raster
STARTUP_F   <- file.path(DATA_DIR, "app_startup.rds")      # all other precomputed objects


# Temporary folder for WBT files
# TEMP_DIR <- file.path(tempdir(), "topmodel_wb")
# dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD PRECOMPUTED DATA
# Replaces all heavy start-up processing. Both files are created by running
# precompute_twi.R and then precompute_data.R locally.
# ═══════════════════════════════════════════════════════════════════════════════
if (!file.exists(STARTUP_F)) {
  stop("Start-up cache not found. Run precompute_twi.R and then precompute_data.R first.")
}
message("Loading precomputed startup cache ...")
cache        <- readRDS(STARTUP_F)
ALL_DATA     <- cache$ALL_DATA       # merged daily precip/flow/temp
TOPIDX       <- cache$TOPIDX         # TWI class table for TOPMODEL
TWI_LAM      <- cache$TWI_LAM        # mean TWI (lambda)
WS_SF_NATIVE <- cache$WS_SF_NATIVE   # watershed boundary polygon
DEM_WS       <- cache$DEM_WS         # clipped DEM for contour lines
rm(cache)  # free memory

# The TWI raster is kept as a separate .tif (not bundled in the .rds)
if (!file.exists(TWI_F)) {
  stop("TWI raster not found. Run precompute_twi.R first.")
}
TWI_RASTER <- rast(TWI_F)
message("  Startup complete.")

# load_daily_data <- function() {
#   # Read & filter each dataset to Watershed 3
#   precip <- read_csv(PRECIP_F, stringsAsFactors = FALSE)
#   precip$DATE <- as.Date(trimws(precip$DATE))
#   precip <- precip[trimws(precip$watershed) == "W3", c("DATE", "Precip")]
#   names(precip) <- c("date", "precip_mm")
# 
#   flow <- read_csv(FLOW_F, stringsAsFactors = FALSE)
#   flow$DATE <- as.Date(trimws(flow$DATE))
#   flow <- flow[flow$WS == 3, c("DATE", "Streamflow")]
#   names(flow) <- c("date", "flow_mm")
#   
#   # Average min/max temperatures to get the daily mean
#   temp <- read_csv(TEMP_F, stringsAsFactors = FALSE)
#   temp$date <- as.Date(trimws(temp$date))
#   temp_avg <- temp |> 
#     group_by(date) |> 
#     summarize(tavg = mean(AVE, na.rm = TRUE), .groups = "drop")
# 
#   # Inner join precip & flow and left join by temp_avg
#   df <- merge(precip, flow, by = "date", all = FALSE)
#   df <- merge(df, temp_avg, by = "date", all.x = TRUE)
#   df <- df[order(df$date), ]

#   # Fill missing temperature values by linear interpolation b/w valid points
#   if (any(is.na(df$tavg))) {
#     idx <- which(!is.na(df$tavg))
#     df$tavg <- approx(idx, df$tavg[idx], seq_len(nrow(df)), rule = 2)$y
#   }
#   df$precip_mm[is.na(df$precip_mm)] <- 0
#   df
# }

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  HAMON PET (mm/day)
# Estimates potential evapotranspiration (PET) from temperature & day length.
# ═══════════════════════════════════════════════════════════════════════════════
hamon_pet <- function(tavg, doy, lat_deg = 43.95) {
  lat_rad  <- lat_deg * pi / 180
  dec      <- 0.4093 * sin(2 * pi * (284 + doy) / 365)           # solar declination
  ws       <- acos(pmax(pmin(-tan(lat_rad) * tan(dec), 1), -1))  # sunset hour angle
  daylight <- 24 * ws / pi                                       # day length (hours)
  esat     <- 0.611 * exp(17.27 * tavg / (tavg + 237.3))         # saturation vapor pressure (kPa)
  pet      <- ifelse(tavg > 0,
                     0.1651 * daylight * esat / (tavg + 273.3) * 29.8, 0)
  pmax(pet, 0)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  TOPMODEL — full simulation with all outputs.
#
# Outputs: simulated flow, baseflow, overland flow, SWE, ET, NSE score, &
# monthly mean deficits for the saturation map animation.
# ═══════════════════════════════════════════════════════════════════════════════
run_topmodel <- function(pars, data, topidx, pet_pre = NULL) {

  n    <- nrow(data)
  P    <- data$precip_mm
  tavg <- data$tavg
  obs  <- data$flow_mm

  # Use the precomputed PET if provided, otherwise calculate from temperature
  if (!is.null(pet_pre)) { pet <- pet_pre
  } else { pet <- hamon_pet(tavg, yday(data$date)) }

  # Parameters 
  qs0     <- pars["qs0"]           # initial subsurface flow (mm/d)
  lnTe    <- pars["lnTe"]          # log of soil transmissivity
  m       <- pars["m"]             # transmissivity decay w/depth
  Sr0     <- pars["Sr0"]           # initial root zone deficit (mm)
  Srmax   <- pars["Srmax"]         # max root zone storage (mm)
  td      <- max(pars["td"], 0.1)  # unsaturated zone time delay (days)
  snow_t  <- pars["snow_t"]        # snow/rain temperature threshold (°C)
  snow_mf <- pars["snow_mf"]       # degree-day melt factor (mm/°C/d)

  # Area-weighted mean TWI (lambda): converts watershed deficit to per-cell thresholds
  twi_vals <- topidx[, 1]
  twi_frac <- topidx[, 2]
  lam      <- sum(twi_vals * twi_frac)  # area-weighted mean TWI

  # Pre-sort TWI thresholds (descending)
  twi_thresh <- m * (twi_vals - lam)
  si  <- order(twi_thresh, decreasing = TRUE)
  tts <- twi_thresh[si]        # sorted thresholds
  tfc <- cumsum(twi_frac[si])  # cumulative area fractions
  nc  <- length(tts)

  # Initialize state variables
  Sd  <- -m * (log(max(qs0, 1e-10)) - lnTe)  # initial saturation deficit (mm)
  if (Sd > 2000) Sd <- 2000; if (Sd < 0) Sd <- 0
  Srz <- Sr0        # root zone storage deficit
  Suz <- 0          # unsaturated zone storage
  SWE <- 0          # snow water equivalent      
  Te  <- exp(lnTe)  # transmissivity at saturation
  td_m <- td / m    # scaled time delay

  # Pre-allocate output arrays
  Qsim <- numeric(n); Qb_out <- numeric(n); Qof_out <- numeric(n)
  SWE_out <- numeric(n); ET_out <- numeric(n); Sd_out <- numeric(n)
  SF_out <- numeric(n)

  # ── Daily time step loop ───────────────────────────────────────────────────
  for (t in seq_len(n)) {
    tv <- tavg[t]
    
    # Snow accumulation & melt (simple degree-day model)
    if (tv <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0        # all precip falls as snow
    } else {
      ma <- snow_mf * (tv - snow_t)       # potential melt
      if (ma > SWE) ma <- SWE             # can't melt more than exists
      SWE <- SWE - ma; rain <- P[t] + ma  # rain + snow melt = effective rainfall
    }
    SWE_out[t] <- SWE

    # Saturated area fraction — binary search on sorted thresholds
    # sf = fraction of watershed that is currently at or above the water table
    if (Sd >= tts[1L]) { sf <- 0
    } else if (Sd < tts[nc]) { sf <- tfc[nc]
    } else {
      lo <- 1L; hi <- nc
      while (hi - lo > 1L) {
        mid <- (lo + hi) %/% 2L
        if (tts[mid] > Sd) lo <- mid else hi <- mid
      }
      sf <- tfc[lo]
    }

    # Saturation-excess overland flow: rain on already-saturated areas
    qof   <- sf * rain
    infil <- rain - qof
    Qof_out[t] <- qof

    # Root zone water balance
    Srz <- Srz - infil                           # infiltration fills root zone
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }  # overflow to unsaturated zone
    if (Srz < Srmax) { ae <- pet[t] * (1 - Srz / Srmax) } else { ae <- 0 }
    Srz <- Srz + ae; if (Srz > Srmax) Srz <- Srmax
    ET_out[t] <- ae

    # Drainage from unsaturated zone to water table (delayed by td)
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd); if (quz > Suz) quz <- Suz; Suz <- Suz - quz
    } else { quz <- 0 }

    # Update watershed-average deficit
    Sd <- Sd - quz; if (Sd < 0) Sd <- 0
    
    # Baseflow: exponential decay from transmissivity
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + 50) qb <- Sd + 50; if (qb < 0) qb <- 0
    Sd <- Sd + qb; if (Sd < 0) Sd <- 0  # deficit increases as water leaves

    Qb_out[t]  <- qb
    Qsim[t]    <- qb + qof  # total flow = baseflow + overland flow
    Sd_out[t]  <- Sd
  }

  Qsim[!is.finite(Qsim)] <- 0
  Qsim <- pmax(Qsim, 0)

  # Nash-Sutcliffe Efficiency (NSE): 1 = perfect, 0 = as good as mean, <0 = poor
  # Skip the first 365 days (warmup), so that the initial conditions don't bias 
  # the score.
  warmup <- min(365, n - 1)
  idx <- (warmup + 1):n
  o <- obs[idx]; s <- Qsim[idx]
  ok <- !is.na(o) & is.finite(s)
  
  if (sum(ok) < 30) { nse <- NA_real_
  } else {
    ss_res <- sum((o[ok] - s[ok])^2)
    ss_tot <- sum((o[ok] - mean(o[ok]))^2)
    nse <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
  }

  list(Qsim = Qsim, Qobs = obs, Qb = Qb_out, Qof = Qof_out,
       SWE = SWE_out, ET = ET_out, Sd = Sd_out, SF = SF_out,
       dates = data$date, nse = nse, warmup = warmup,
       m = m, lam = lam)
  }

# ═══════════════════════════════════════════════════════════════════════════════
# 3b.   LIGHTWEIGHT NSE-ONLY TOPMODEL (optimizer only — no arrays)
#
# Same idea as run_topmodel(), but skips storing the output arrays. Returns
# negative NSE (minimization convention for optim()).
# ═══════════════════════════════════════════════════════════════════════════════
topmodel_nse_only <- function(par_vec, par_names, P, tavg, obs, pet,
                              topidx, n, warmup, lower, upper) {

  if (any(!is.finite(par_vec))) return(1e6)
  if (any(par_vec < lower) || any(par_vec > upper)) return(1e6)

  qs0 <- par_vec[1]; lnTe <- par_vec[2]; m <- par_vec[3]; Sr0 <- par_vec[4]
  Srmax <- par_vec[5]; td <- par_vec[6]; snow_t <- par_vec[7]; snow_mf <- par_vec[8]
  if (td < 0.1) td <- 0.1

  twi_vals <- topidx[, 1]; twi_frac <- topidx[, 2]
  lam <- sum(twi_vals * twi_frac)
  twi_thresh <- m * (twi_vals - lam)
  si  <- order(twi_thresh, decreasing = TRUE)
  tts <- twi_thresh[si]
  tfc <- cumsum(twi_frac[si])
  nc  <- length(tts)

  Sd <- -m * (log(max(qs0, 1e-10)) - lnTe)
  if (Sd > 2000) Sd <- 2000; if (Sd < 0) Sd <- 0
  Srz <- Sr0; Suz <- 0; SWE <- 0
  Te <- exp(lnTe); td_m <- td / m

  obs_after <- obs[(warmup + 1):n]
  obs_ok <- obs_after[!is.na(obs_after)]
  if (length(obs_ok) < 30) return(1e6)
  obs_mean <- mean(obs_ok)

  # Accumulate residuals without storing arrays
  ss_res <- 0; ss_obs <- 0; n_ok <- 0L

  for (t in seq_len(n)) {
    tv <- tavg[t]
    if (tv <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0
    } else {
      ma <- snow_mf * (tv - snow_t)
      if (ma > SWE) ma <- SWE
      SWE <- SWE - ma; rain <- P[t] + ma
    }

    if (Sd >= tts[1L]) { sf <- 0
    } else if (Sd < tts[nc]) { sf <- tfc[nc]
    } else {
      lo <- 1L; hi <- nc
      while (hi - lo > 1L) {
        mid <- (lo + hi) %/% 2L
        if (tts[mid] > Sd) lo <- mid else hi <- mid
      }
      sf <- tfc[lo]
    }

    infil <- (1 - sf) * rain
    Srz <- Srz - infil
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }
    if (Srz < Srmax) { ae <- pet[t] * (1 - Srz / Srmax) } else { ae <- 0 }
    Srz <- Srz + ae; if (Srz > Srmax) Srz <- Srmax

    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd); if (quz > Suz) quz <- Suz; Suz <- Suz - quz
    } else { quz <- 0 }

    Sd <- Sd - quz; if (Sd < 0) Sd <- 0
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + 50) qb <- Sd + 50; if (qb < 0) qb <- 0
    Sd <- Sd + qb; if (Sd < 0) Sd <- 0

    qsim <- qb + sf * rain

    # Only score post-warmup days with valid observations
    if (t > warmup) {
      o_t <- obs[t]
      if (!is.na(o_t) && is.finite(qsim)) {
        ss_res <- ss_res + (o_t - qsim)^2
        ss_obs <- ss_obs + (o_t - obs_mean)^2
        n_ok <- n_ok + 1L
      }
    }
  }

  if (n_ok < 30 || ss_obs <= 0) return(1e6)
  -(1 - ss_res / ss_obs)  # negative NSE (optim() minimizes)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  OPTIMIZER — Nelder-Mead on a calibration window
#
# Tries three parameter starting points and keeps the best result. Uses only 
# recent data (cal_years + 1 year warmup) to keep the run time reasonable.
# ═══════════════════════════════════════════════════════════════════════════════
optimize_params <- function(data, topidx, cal_years = 3) {

  par_names <- c("qs0","lnTe","m","Sr0","Srmax","td","snow_t","snow_mf")
  lower <- c(0.01, -7,   5,  0,   5,  1, -2, 1)
  upper <- c(10,    5, 100, 50, 300, 60,  2, 6)
  names(lower) <- par_names
  names(upper) <- par_names

  # Subset to calibration window (most recent cal_years + 1yr warmup)
  n_full   <- nrow(data)
  cal_days <- cal_years * 365 + 365
  if (n_full > cal_days) {
    cal_data <- data[(n_full - cal_days + 1):n_full, ]
  } else { cal_data <- data }

  # Precompute vectors once & pass into each optim() call
  cal_n      <- nrow(cal_data)
  cal_P      <- cal_data$precip_mm
  cal_tavg   <- cal_data$tavg
  cal_obs    <- cal_data$flow_mm
  cal_pet    <- hamon_pet(cal_tavg, yday(cal_data$date))
  cal_warmup <- min(365, cal_n - 1)

  obj_fast <- function(par_vec) {
    tryCatch(
      topmodel_nse_only(par_vec, par_names, cal_P, cal_tavg, cal_obs,
                        cal_pet, topidx, cal_n, cal_warmup, lower, upper),
      error = function(e) 1e6
    )
  }

  # Three starting points to avoid getting stuck in a local minimum
  starts <- list(
    c(qs0 = 1.0, lnTe =  1.0, m = 30, Sr0 =  2, Srmax = 100, td = 10, snow_t =  0.0, snow_mf = 3.0),
    c(qs0 = 2.0, lnTe =  2.0, m = 20, Sr0 =  1, Srmax =  60, td =  5, snow_t =  1.0, snow_mf = 2.0),
    c(qs0 = 3.0, lnTe = -2.0, m = 40, Sr0 = 10, Srmax = 150, td = 20, snow_t = -0.5, snow_mf = 3.5)
  )

  best_val  <- Inf
  best_pars <- starts[[1]]

  for (i in seq_along(starts)) {
    if (!is.null(progress_fn))
      progress_fn(i, length(starts), sprintf("Start %d of %d \u2026", i, length(starts)))
    opt <- tryCatch(
      optim(starts[[i]], obj_fast, method = "Nelder-Mead",
            control = list(maxit = 400, reltol = 1e-6)),
      error = function(e) NULL)
    if (!is.null(opt) && is.finite(opt$value) && opt$value < best_val) {
      best_val <- opt$value; best_pars <- opt$par
    }
  }

  names(best_pars) <- par_names
  pmax(pmin(best_pars, upper), lower) # clamp to valid parameter bounds
}

# # ═══════════════════════════════════════════════════════════════════════════════
# # 5.  PRECOMPUTE AT STARTUP
# #
# # These objects are created once the app launches and are shared across all
# # user sessions (global scope in Shiny = computed before server() runs).
# # ═══════════════════════════════════════════════════════════════════════════════
# message("Loading daily data ...")
# ALL_DATA <- load_daily_data()
# message(sprintf("  %d records: %s to %s",
#                 nrow(ALL_DATA), min(ALL_DATA$date), max(ALL_DATA$date)))
# 
# message("Computing TWI ...")
# TWI_PRECOMPUTED <- file.path(DATA_DIR, "twi_precomputed.tif")
# if (file.exists(TWI_PRECOMPUTED)) {
#   # Fast path: precomputed file exists (always the case on shinyapps.io)
#   message("  Loading precomputed TWI ...")
#   TWI_RASTER <- rast(TWI_PRECOMPUTED)
# } else {
#   # Slow path: compute from scratch locally (requires WBT)
#   if (!HAS_WHITEBOX) stop("WhiteBox is not installed and there is no precomputed TWI found.
#                           Run precompute_twi.R first.")
#   message("  Computing via WhiteBox Tools (FD8) ...")
#   TWI_RASTER <- compute_twi_raster()
#   writeRaster(TWI_RASTER, TWI_PRECOMPUTED, overwrite = TRUE)
#   message("  Saved precomputed TWI to ", TWI_PRECOMPUTED)
# }
# 
# # Build the TOPIDX (16 TWI classes → area fractions)
# TOPIDX     <- make_topidx_classes(TWI_RASTER, n_classes = 16)
# TWI_LAM    <- sum(TOPIDX[, 1] * TOPIDX[, 2])  # mean TWI (lambda)
# message("  Lambda (mean TWI): ", round(TWI_LAM, 2))
# 
# # Load watershed boundary & re-project to match the TWI raster CRS
# WS_SF_NATIVE <- st_read(SHP_F, quiet = TRUE)
# WS_SF_NATIVE <- st_transform(WS_SF_NATIVE, crs(TWI_RASTER))
# 
# # Clip DEM to watershed boundary for drawing elevation contours on maps
# DEM_WS <- mask(crop(rast(DEM_F),
#                      vect(st_transform(WS_SF_NATIVE, crs(rast(DEM_F))))),
#                vect(st_transform(WS_SF_NATIVE, crs(rast(DEM_F)))))
# 
# message("  Startup complete.")

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: snap a value to the nearest slider step
# ═══════════════════════════════════════════════════════════════════════════════
snap_to_step <- function(val, min_val, max_val, step) {
  snapped <- round((val - min_val) / step) * step + min_val
  max(min(snapped, max_val), min_val)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  SHINY UI
# ═══════════════════════════════════════════════════════════════════════════════
# Custom CSS for the NSE score box & sidebar parameter groups
app_css <- "
  body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
  .nse-box {
    text-align: center; padding: 14px; border-radius: 10px;
    font-size: 28px; font-weight: 700; margin-bottom: 12px;
  }
  .nse-good  { background: #d4edda; color: #155724; border: 2px solid #28a745; }
  .nse-ok    { background: #fff3cd; color: #856404; border: 2px solid #ffc107; }
  .nse-bad   { background: #f8d7da; color: #721c24; border: 2px solid #dc3545; }
  .pg { background: #f8f9fa; border-radius: 8px;
        padding: 12px 14px 4px 14px; margin-bottom: 10px; }
  .pg h6 { margin-top: 0; color: #495057; font-weight: 600; }
  .action-btn { width: 100%; margin-top: 6px; margin-bottom: 6px; }
  .opt-params { font-size: 12px; color: #495057; margin-top: 4px;
                padding: 6px 8px; background: #e9ecef; border-radius: 6px;
                font-family: monospace; line-height: 1.6; }
  .pg .help-block { font-size: 11px; color: #888;
                    margin-top: 2px; margin-bottom: 2px; line-height: 1.3; }
  .param-section h6 { cursor: pointer; }
"

ui <- page_sidebar(
  title = "HBEF Watershed 3 \u2014 TOPMODEL Explorer",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  tags$head(tags$style(HTML(app_css))),

  # ── Left sidebar: simulation controls & parameter sliders ────────────────
  sidebar = sidebar(
    width = 330,

    # Year range selector: filters ALL_DATA to the chosen window
    sliderInput("year_range", "Viewing Period",
                min   = as.integer(format(min(ALL_DATA$date), "%Y")),
                max   = as.integer(format(max(ALL_DATA$date), "%Y")),
                value = c(2000, 2001),
                step  = 1, sep = "",
                ticks = TRUE),
    helpText("Select the time window to view. The model runs on the selected period."),

    actionButton("run_model", "\u25B6  Run Model",
                 class = "btn-success action-btn"),

    # NSE score box (green/yellow/red based on the value) & an optimized parameter readout
    uiOutput("nse_display"),
    uiOutput("opt_params_display"),

    # Subsurface parameters: control baseflow & drainage
    div(class = "pg",
      h6("\u2193 Subsurface Flow Parameters"),
      
      helpText("Initial subsurface flow at time zero. Sets the starting saturation deficit."),
      sliderInput("qs0",  "qs0 \u2013 init. flow (mm/d)",
                  min = 0.01, max = 10, value = 5.5, step = 0.01),
      
      helpText("Log of soil transmissivity when fully saturated. Controls baseflow magnitude \u2014 higher = more baseflow."),
      sliderInput("lnTe", "ln(Te) \u2013 transmissivity",
                  min = -7, max = 5, value = 1.9, step = 0.1),
      
      helpText("Rate transmissivity declines with depth. Small m = flashy response; large m = slow drainage."),
      sliderInput("m",    "m \u2013 decay param (mm)",
                  min = 5, max = 100, value = 19, step = 1),
      
      helpText("Time delay for water to drain through the unsaturated zone to the water table."),
      sliderInput("td",   "td \u2013 unsat. delay (days)",
                  min = 1, max = 60, value = 3, step = 1),
    ),

    # Root zone parameters: control ET & soil moisture
    div(class = "pg",
      h6("\U0001F331 Root Zone  & Evapotranspiration"),
      
      helpText("Initial root zone storage deficit. How dry the soil is at the start of simulation."),
      sliderInput("Sr0",   "Sr0 \u2013 init. deficit (mm)",
                  min = 0, max = 50, value = 2.5, step = 0.5),
      helpText("Maximum water the root zone can hold before draining. Controls ET \u2014 larger = more ET in summer."),
      sliderInput("Srmax", "Srmax \u2013 max capacity (mm)",
                  min = 5, max = 300, value = 65, step = 5),
    ),

    # Snow parameters: degree-day accumulation & melt
    div(class = "pg",
      h6("\u2744 Snow Accumulation & Melt"),
      
      helpText("Temperature threshold: below this, precipitation falls as snow."),
      sliderInput("snow_t",  "Tmelt (\u00B0C)",
                  min = -2, max = 2, value = -0.7, step = 0.1),
      helpText("Degree-day melt factor. mm of snowmelt per degree above Tmelt per day."),
      sliderInput("snow_mf", "Melt coeff (mm/\u00B0C/d)",
                  min = 1, max = 6, value = 3.4, step = 0.1),
    ),

    actionButton("optimize", "\u26A1  Optimize Parameters",
                 class = "btn-primary action-btn"),
    helpText("Automatically finds the best parameter values by minimizing the
             difference between observed and simulated streamflow. Uses the 
             Nelder-Mead algorithm, a derivative-free method that works well 
             for non-smooth hydrological models.")
  ),

  # ── Main panel: tabbed plots ───────────────────────────────────────────────
  navset_card_tab(
    nav_panel("Model Output",    plotlyOutput("model_output", height = "700px")),
    nav_panel("Observed vs Simulated Discharge",
              plotlyOutput("scatter", height = "520px")),
    nav_panel("Flow Duration",   plotlyOutput("fdc", height = "520px")),
    nav_panel("TWI Map",         plotOutput("twi_map", height = "560px")),
    nav_panel("Saturation Map",
              fluidRow(
                column(8,
                       sliderInput("sat_day", "Day", min = 1, max = 730, value = 1,
                                   step = 7, width = "100%",
                                   animate = animationOptions(interval = 1000, loop = FALSE)),
                       plotOutput("sat_map", height = "460px")
                ),
                column(4,
                       plotlyOutput("sat_discharge", height = "220px"),
                       plotlyOutput("sat_fraction", height = "220px")
                )
              ))
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  SHINY SERVER
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # Reactive values shared across observers & outputs
  model_result   <- reactiveVal(NULL)  # stores the last run_topmodel() result
  opt_params_txt <- reactiveVal(NULL)  # stores formatted optimized param string

  # Filter ALL_DATA to the selected year range
  sim_data <- reactive({
    req(input$year_range)
    start_date <- as.Date(paste0(input$year_range[1], "-01-01"))
    end_date   <- as.Date(paste0(input$year_range[2], "-12-31"))
    d <- ALL_DATA[ALL_DATA$date >= start_date & ALL_DATA$date <= end_date, ]
    validate(need(nrow(d) > 400, "Select a wider date range (> 1 year)."))
    d
  })

  # Bundle current slider values into a named parameter vector
  current_pars <- reactive({
    c(qs0 = input$qs0, lnTe = input$lnTe, m = input$m,
      Sr0 = input$Sr0, Srmax = input$Srmax, td = input$td,
      snow_t = input$snow_t, snow_mf = input$snow_mf)
  })

  # ── Run Model button ────────────────────────────────────────────────────────
  observeEvent(input$run_model, {
    d <- tryCatch(sim_data(), error = function(e) NULL)
    
    if (is.null(d)) {
      showNotification("Check date range.", type = "error"); return()
      }
    res <- tryCatch(run_topmodel(current_pars(), d, TOPIDX),
                    error = function(e) {
                      showNotification(paste("Error:", e$message),
                                       type = "error"); NULL
                      })
    
    if (!is.null(res)) {
      model_result(res)
      nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
      showNotification(sprintf("NSE = %s", nse_str), type = "message", duration = 3)
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)

  # ── Auto-run on first load so that the plots aren't blank at start-up ────────
  ran_once <- reactiveVal(FALSE)
  observe({
    req(!ran_once())
    d <- tryCatch(sim_data(), error = function(e) NULL)
    if (!is.null(d)) {
      res <- tryCatch(run_topmodel(isolate(current_pars()), d, TOPIDX),
                      error = function(e) NULL)
      if (!is.null(res)) model_result(res)
      ran_once(TRUE)
    }
  })

  # ── Optimize button ─────────────────────────────────────────────────────────
  observeEvent(input$optimize, {
    data <- tryCatch(sim_data(), error = function(e) NULL)
    if (is.null(data)) {
      showNotification("Check date range.", type = "error"); return()
      }
    
    showNotification("Optimizing \u2026 Start 1 of 3",
                     type = "message", duration = NULL, id = "opt_prog")
    
    best <- optimize_params(data, TOPIDX, progress_fn = function(i, n, msg) {
      showNotification(sprintf("Optimizing \u2026 %s", msg),
                       type = "message", duration = NULL, id = "opt_prog")
    })
    
    removeNotification("opt_prog")

    # Push the optimized values back to the sliders
    updateSliderInput(session, "qs0",     value = snap_to_step(best["qs0"],     0.01, 10,  0.01))
    updateSliderInput(session, "lnTe",    value = snap_to_step(best["lnTe"],    -7,   5,   0.1))
    updateSliderInput(session, "m",       value = snap_to_step(best["m"],        5,   100,  1))
    updateSliderInput(session, "td",      value = snap_to_step(best["td"],       1,   60,   1))
    updateSliderInput(session, "Sr0",     value = snap_to_step(best["Sr0"],      0,   50,   0.5))
    updateSliderInput(session, "Srmax",   value = snap_to_step(best["Srmax"],    5,   300,  5))
    updateSliderInput(session, "snow_t",  value = snap_to_step(best["snow_t"],  -2,   2,    0.1))
    updateSliderInput(session, "snow_mf", value = snap_to_step(best["snow_mf"],  1,   6,    0.1))

    # Display exact (un-snapped) optimized values for reference
    opt_params_txt(sprintf(
      "qs0=%.3f  lnTe=%.2f  m=%.1f\nSr0=%.1f  Srmax=%.1f  td=%.1f\nTmelt=%.2f  Melt=%.2f",
      best["qs0"], best["lnTe"], best["m"],
      best["Sr0"], best["Srmax"], best["td"],
      best["snow_t"], best["snow_mf"]))

    # Re-run the full simulation with the exact optimized parameters
    res <- tryCatch(run_topmodel(best, data, TOPIDX), error = function(e) NULL)
    if (!is.null(res)) {
      model_result(res)
      nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
      showNotification(sprintf("Optimized! NSE = %s", nse_str),
                       type = "message", duration = 5)
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # ── NSE score box: green ≥ 0.5, yellow ≥ 0, red < 0 ───────────────────────
  output$nse_display <- renderUI({
    res <- model_result()
    if (is.null(res) || is.na(res$nse) || !is.finite(res$nse)) {
      div(class = "nse-box nse-bad", "NSE: \u2014")
      } else {
        cls <- if (res$nse >= 0.5) "nse-good"
        else if (res$nse >= 0) "nse-ok" else "nse-bad"
        div(class = paste("nse-box", cls), sprintf("NSE = %.3f", res$nse))
        }
    })

  # ── Show exact optimized parameter values below the NSE box ────────────────
  output$opt_params_display <- renderUI({
    txt <- opt_params_txt()
    if (is.null(txt)) return(NULL)
    div(class = "opt-params",
        tags$strong("Optimized values:"), tags$br(),
        tags$pre(style = "margin:2px 0 0 0; white-space:pre-wrap;", txt))
  })

  # ════════════════════════════════════════════════════════════════════════════
  # MODEL OUTPUT: stacked subplots with shared x-axis
  # ════════════════════════════════════════════════════════════════════════════
  
  # Model output plots
  output$model_output <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    
    d <- data.frame(date = res$dates, Qobs = res$Qobs, Qsim = res$Qsim,
                    Qof = res$Qof,    Qb = res$Qb,     ET = res$ET,
                    SWE = res$SWE,    Sd = res$Sd,     SF = res$SF * 100)
    nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
    
    # 1. Discharge
    p1 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qobs, name = "Observed", line = list(color = "#2980b9", width = 1.2)) |> 
      add_lines(y = ~Qsim, name = "Simulated", line = list(color = "#e74c3c", width = 1.2)) |> 
      layout(yaxis = list(title = "Q (mm/d)"),
             annotations = list(text = paste0("NSE = ", nse_str),
                                x = 0.02, y = 0.95, xref = "paper", yref = "paper",
                                showarrow = FALSE, font = list(size = 12)))
    
    # 2. Surface runoff
    p2 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qof, name = "Surface Runoff", line = list(color = "#e67e22", width = 1)) |> 
      layout(yaxis = list(title = "Qof (mm/d)"))
    
    # 3. Subsurface flow
    p3 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qb, name = "Subsurface Flow", line = list(color = "#8e44ad", width = 1)) |> 
      layout(yaxis = list(title = "Qb (mm/d)"))
    
    # 4. ET
    p4 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~ET, name = "Actual ET", line = list(color = "#27ae60", width = 1)) |> 
      layout(yaxis = list(title = "ET (mm/d)"))
    
    # 5. SWE
    p5 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~SWE, name = "SWE", line = list(color = "#3498db", width = 1)) |> 
      layout(yaxis = list(title = "SWE (mm)"))
    
    # 6. Storage deficit
    p6 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~Sd, name = "Storage Deficit", line = list(color = "#c0392b", width = 1)) |> 
      layout(yaxis = list(title = "Sd (mm)"))
    
    # 7. Saturated fraction
    p7 <- plot_ly(d, x = ~date) |> 
      add_lines(y = ~SF, name = "Saturated %", line = list(color = "#2166ac", width = 1)) |> 
      layout(yaxis = list(title = "Sat (%)"),
             xaxis = list(title = ""))
    
    subplot(p1, p2, p3, p4, p5, p6, p7,
            nrows = 7, shareX = TRUE, titleY = TRUE,
            heights = c(0.22, 0.11, 0.11, 0.11, 0.15, 0.15, 0.15)) |> 
      layout(title = list(text = "Model Output", font = list(size = 16)),
             showlegend = FALSE,
             hovermode = "x unified",
             margin = list(t = 50, b = 30, l = 60))
    })

  # Observed vs simulated scatter: perfect model would follow the 1:1 line
  output$scatter <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    idx <- (res$warmup + 1):length(res$Qsim)
    d <- data.frame(obs = res$Qobs[idx], sim = res$Qsim[idx])
    d <- d[!is.na(d$obs) & is.finite(d$sim), ]
    validate(need(nrow(d) > 0, "No valid data."))
    rng <- range(c(d$obs, d$sim), na.rm = TRUE)
    plot_ly(d, x = ~obs, y = ~sim, type = "scattergl", mode = "markers",
            marker = list(size = 3, color = "#3498db", opacity = 0.4)) |> 
      add_lines(x = rng, y = rng, name = "1:1",
                line = list(color = "black", dash = "dash", width = 1)) |> 
      layout(title = list(text = "Observed vs Simulated Discharge",
                          font = list(size = 15)),
             xaxis = list(title = "Observed (mm/d)"),
             yaxis = list(title = "Simulated (mm/d)"),
             showlegend = FALSE, margin = list(t = 50))
  })

  # Flow Duration Curve: log-scale exceedance plot; shows high/low flow fit
  output$fdc <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    idx   <- (res$warmup + 1):length(res$Qsim)
    obs_v <- sort(res$Qobs[idx][!is.na(res$Qobs[idx])], decreasing = TRUE)
    sim_v <- sort(res$Qsim[idx][is.finite(res$Qsim[idx])], decreasing = TRUE)
    exc_o <- seq_along(obs_v) / length(obs_v) * 100
    exc_s <- seq_along(sim_v) / length(sim_v) * 100
    plot_ly() |> 
      add_lines(x = exc_o, y = obs_v, name = "Observed",
                line = list(color = "#2980b9", width = 1.5)) |> 
      add_lines(x = exc_s, y = sim_v, name = "Simulated",
                line = list(color = "#e74c3c", width = 1.5)) |> 
      layout(title = list(text = "Flow-Duration Curve", font = list(size = 15)),
             xaxis = list(title = "Exceedance (%)"),
             yaxis = list(title = "Discharge (mm/d)", type = "log"),
             legend = list(orientation = "h", y = -0.12), margin = list(t = 50))
  })

  # Snow Water Equivalent and actual ET on dual y-axes
  output$snow_et <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    d <- data.frame(date = res$dates, SWE = res$SWE, ET = res$ET)
    plot_ly(d, x = ~date) |> 
      add_lines(y = ~SWE, name = "SWE (mm)",
                line = list(color = "#3498db", width = 1.3),
                yaxis = "y") |> 
      add_lines(y = ~ET,  name = "Actual ET (mm/d)",
                line = list(color = "#27ae60", width = 1.2),
                yaxis = "y2") |> 
      layout(title = list(text = "Snow Water Equivalent & Actual ET",
                          font = list(size = 15)),
             xaxis = list(title = ""),
             yaxis  = list(title = "SWE (mm)", side = "left",
                           titlefont = list(color = "#3498db"),
                           tickfont  = list(color = "#3498db")),
             yaxis2 = list(title = "Actual ET (mm/d)", side = "right",
                           overlaying = "y", showgrid = FALSE,
                           titlefont = list(color = "#27ae60"),
                           tickfont  = list(color = "#27ae60")),
             legend = list(orientation = "h", y = -0.12),
             hovermode = "x unified", margin = list(t = 50))
  })

  # Static TWI map: shows landscape wetness potential (fixed, not model-dependent)
  output$twi_map <- renderPlot({
    par(mar = c(2, 2, 3, 4), bg = "white")
    plot(TWI_RASTER,
         col = hcl.colors(50, "YlGnBu", rev = TRUE),
         axes = FALSE,
         main = sprintf("Topographic Wetness Index  (FD8, \u03BB = %.1f)", TWI_LAM),
         cex.main = 1.3)
    plot(st_geometry(WS_SF_NATIVE), add = TRUE,
         border = "#222222", lwd = 2.5, col = NA)

    # Elevation contour lines
    contour(raster::raster(DEM_WS), add = TRUE,
            nlevels = 10, col = adjustcolor("#333333", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
    mtext("Contour elevation is in meters", side = 1, line = 0.5,
          cex = 0.8, col = "#666666")
  }, res = 120)

  # Update the saturation map day slider range when a new model result is available
  observe({
    res <- model_result()
    req(!is.null(res))
    n_days <- length(res$Sd)
    updateSliderInput(session, "sat_day", min = 1, max = n_days, value = 1)
  })

  # Animated saturation map: binary saturated/unsaturated based on a daily
  # mean deficit vs. the TWI threshold for each cell. Cells with TWI >= 
  # threshold are predicted to be at or above the water table (saturated).
  output$sat_map <- renderPlot({
    res <- model_result()
    validate(need(!is.null(res), "Run or optimize the model first."))
    
    day_idx <- input$sat_day
    req(day_idx)
    n_days <- length(res$Sd)
    if (day_idx > n_days) day_idx <- n_days
    
    # Convert daily deficit to a spatial saturation threshold
    Sd  <- res$Sd[day_idx]
    m   <- res$m
    lam <- res$lam
    thresh <- lam + Sd / m
    label <- format(res$dates[day_idx], "%Y-%m-%d")
    
    # Classify each raster cell as saturated (1) or unsaturated (0)
    twi_v <- values(TWI_RASTER)
    if (is.matrix(twi_v)) twi_v <- twi_v[, 1]
    sat_v <- rep(NA_real_, length(twi_v))
    valid <- !is.na(twi_v)
    sat_v[valid] <- ifelse(twi_v[valid] >= thresh, 1, 0)
    
    sat_rast <- TWI_RASTER
    values(sat_rast) <- sat_v
    
    pct_sat <- round(100 * sum(sat_v == 1, na.rm = TRUE) / sum(valid), 1)
    
    par(mar = c(2, 2, 3, 1), bg = "white")
    plot(sat_rast,
         col = c("#f5e6c8", "#2166ac"),
         breaks = c(-0.5, 0.5, 1.5),
         legend = FALSE, axes = FALSE,
         main = sprintf("%s  \u2014  %.1f%% saturated", label, pct_sat),
         cex.main = 1.2)
    plot(st_geometry(WS_SF_NATIVE), add = TRUE,
         border = "#222222", lwd = 2.5, col = NA)

    # Elevation contour lines for terrain context
    contour(raster::raster(DEM_WS), add = TRUE,
            nlevels = 10, col = adjustcolor("#555555", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
    legend("bottomright",
           legend = c("Saturated", "Unsaturated"),
           fill = c("#2166ac", "#f5e6c8"),
           border = c("#2166ac", "#f5e6c8"),
           bty = "n", cex = 0.9)
    mtext("Contour elevations in meters", side = 1, line = 0.5,
          cex = 0.7, col = "#666666")
  }, res = 110)

  # Discharge time series with marker showing current day
  output$sat_discharge <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), ""))
    
    day_idx <- input$sat_day
    req(day_idx)
    n_days <- length(res$Sd)
    if (day_idx > n_days) day_idx <- n_days
    
    d <- data.frame(date = res$dates, Qobs = res$Qobs, Qsim = res$Qsim)
    
    plot_ly(d, x = ~date) |> 
      add_lines(y = ~Qobs, name = "Observed", line = list(color = "#2980b9", width = 1)) |> 
      add_lines(y = ~Qsim, name = "Simulated", line = list(color = "#e74c3c", width = 1)) |> 
      add_markers(x = res$dates[day_idx], y = res$Qsim[day_idx],
                  name = "Current", marker = list(color = "black", size = 8)) |> 
      layout(title = list(text = "Discharge", font = list(size = 12)),
             yaxis = list(title = "mm/d", titlefont = list(size = 10)),
             xaxis = list(title = ""),
             showlegend = FALSE,
             margin = list(t = 30, b = 20, l = 50, r = 10),
             hovermode = "x unified")
  })
  
  # Saturated fraction time series with marker
  output$sat_fraction <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), ""))
    
    day_idx <- input$sat_day
    req(day_idx)
    n_days <- length(res$Sd)
    if (day_idx > n_days) day_idx <- n_days
    
    d <- data.frame(date = res$dates, SF = res$SF * 100)
    
    plot_ly(d, x = ~date) |> 
      add_lines(y = ~SF, name = "Saturated %", line = list(color = "#2166ac", width = 1)) |> 
      add_markers(x = res$dates[day_idx], y = d$SF[day_idx],
                  name = "Current", marker = list(color = "black", size = 8)) |> 
      layout(title = list(text = "Saturated Catchment Area", font = list(size = 12)),
             yaxis = list(title = "%", titlefont = list(size = 10)),
             xaxis = list(title = ""),
             showlegend = FALSE,
             margin = list(t = 30, b = 20, l = 50, r = 10),
             hovermode = "x unified")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  LAUNCH
# ═══════════════════════════════════════════════════════════════════════════════
shinyApp(ui, server)

