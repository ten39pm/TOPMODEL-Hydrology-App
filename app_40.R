library(shiny)
library(bslib)
library(terra)
library(raster)
library(plotly)
library(dplyr)
library(lubridate)
library(sf)
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
DATA_DIR <- "data"
PRECIP_F <- file.path(DATA_DIR, "DailyWatershed.csv")
FLOW_F   <- file.path(DATA_DIR, "HBEF_DailyStreamflow_1956-2024.csv")
TEMP_F   <- file.path(DATA_DIR, "HBEF_air_temp_daily.csv")
DEM_F    <- file.path(DATA_DIR, "dem_10m.tif")
SHP_F    <- file.path(DATA_DIR, "Watershed3HB.shp")

TEMP_DIR <- file.path(tempdir(), "topmodel_wb")
dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD & MERGE DAILY DATA
# ═══════════════════════════════════════════════════════════════════════════════
load_daily_data <- function() {
  precip <- read.csv(PRECIP_F, stringsAsFactors = FALSE)
  precip$DATE <- as.Date(trimws(precip$DATE))
  precip <- precip[trimws(precip$watershed) == "W3", c("DATE", "Precip")]
  names(precip) <- c("date", "precip_mm")

  flow <- read.csv(FLOW_F, stringsAsFactors = FALSE)
  flow$DATE <- as.Date(trimws(flow$DATE))
  flow <- flow[flow$WS == 3, c("DATE", "Streamflow")]
  names(flow) <- c("date", "flow_mm")

  temp <- read.csv(TEMP_F, stringsAsFactors = FALSE)
  temp$date <- as.Date(trimws(temp$date))
  temp_avg <- temp %>%
    group_by(date) %>%
    summarise(tavg = mean(AVE, na.rm = TRUE), .groups = "drop")

  df <- merge(precip, flow, by = "date", all = FALSE)
  df <- merge(df, temp_avg, by = "date", all.x = TRUE)
  df <- df[order(df$date), ]

  if (any(is.na(df$tavg))) {
    idx <- which(!is.na(df$tavg))
    df$tavg <- approx(idx, df$tavg[idx], seq_len(nrow(df)), rule = 2)$y
  }
  df$precip_mm[is.na(df$precip_mm)] <- 0
  df
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  COMPUTE TWI via WhiteBox Tools
#
#     Pipeline:
#       1. Clip DEM to watershed boundary
#       2. Breach depressions (least-cost)
#       3. FD8 specific catchment area (WhiteBox — multi-directional flow)
#       4. Slope in degrees (WhiteBox)
#       5. TWI = ln(SCA / tan(slope))
# ═══════════════════════════════════════════════════════════════════════════════
compute_twi_raster <- function() {
  dem <- rast(DEM_F)
  ws  <- st_read(SHP_F, quiet = TRUE)
  ws_v <- vect(st_transform(ws, crs(dem)))
  dem_ws <- mask(crop(dem, ws_v), ws_v)

  dem_clip_f <- file.path(TEMP_DIR, "dem_clip.tif")
  breached_f <- file.path(TEMP_DIR, "dem_breached.tif")
  sca_f      <- file.path(TEMP_DIR, "sca.tif")
  slope_f    <- file.path(TEMP_DIR, "slope.tif")
  twi_f      <- file.path(TEMP_DIR, "twi.tif")

  writeRaster(dem_ws, dem_clip_f, overwrite = TRUE)

  # Breach depressions to ensure drainage connectivity
  wbt_breach_depressions_least_cost(
    dem    = dem_clip_f,
    output = breached_f,
    dist   = 10,
    fill   = TRUE
  )

  # FD8 specific catchment area (distributes flow to all downslope neighbors)
  wbt_fd8_flow_accumulation(
    dem      = breached_f,
    output   = sca_f,
    out_type = "specific contributing area"
  )

  # Slope in degrees via WhiteBox
  wbt_slope(
    dem    = breached_f,
    output = slope_f,
    units  = "degrees"
  )

  # Compute TWI = ln(SCA / tan(slope))
  sca <- rast(sca_f)
  slp <- rast(slope_f)

  # Convert slope to radians, clamp to avoid division by zero
  slp_rad <- slp * pi / 180
  tan_b   <- tan(slp_rad)
  tan_b[tan_b < 0.001] <- 0.001

  twi <- log(sca / tan_b)
  mask(twi, dem_ws)
}

make_topidx_classes <- function(twi_raster, n_classes = 16) {
  vals <- values(twi_raster, na.rm = TRUE)[, 1]
  brks <- seq(min(vals), max(vals), length.out = n_classes + 1)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  cnts <- hist(vals, breaks = brks, plot = FALSE)$counts
  frac <- cnts / sum(cnts)
  cbind(twi = mids, frac = frac)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  HAMON PET (mm/day)
# ═══════════════════════════════════════════════════════════════════════════════
hamon_pet <- function(tavg, doy, lat_deg = 43.95) {
  lat_rad  <- lat_deg * pi / 180
  dec      <- 0.4093 * sin(2 * pi * (284 + doy) / 365)
  ws       <- acos(pmax(pmin(-tan(lat_rad) * tan(dec), 1), -1))
  daylight <- 24 * ws / pi
  esat     <- 0.611 * exp(17.27 * tavg / (tavg + 237.3))
  pet      <- ifelse(tavg > 0,
                     0.1651 * daylight * esat / (tavg + 273.3) * 29.8, 0)
  pmax(pet, 0)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4.  TOPMODEL — full output version
# ═══════════════════════════════════════════════════════════════════════════════
run_topmodel <- function(pars, data, topidx, pet_pre = NULL) {

  n    <- nrow(data)
  P    <- data$precip_mm
  tavg <- data$tavg
  obs  <- data$flow_mm

  if (!is.null(pet_pre)) { pet <- pet_pre
  } else { pet <- hamon_pet(tavg, yday(data$date)) }

  qs0     <- pars["qs0"]
  lnTe    <- pars["lnTe"]
  m       <- pars["m"]
  Sr0     <- pars["Sr0"]
  Srmax   <- pars["Srmax"]
  td      <- max(pars["td"], 0.1)
  snow_t  <- pars["snow_t"]
  snow_mf <- pars["snow_mf"]

  twi_vals <- topidx[, 1]
  twi_frac <- topidx[, 2]
  lam      <- sum(twi_vals * twi_frac)

  twi_thresh <- m * (twi_vals - lam)
  si  <- order(twi_thresh, decreasing = TRUE)
  tts <- twi_thresh[si]
  tfc <- cumsum(twi_frac[si])
  nc  <- length(tts)

  Sd  <- -m * (log(max(qs0, 1e-10)) - lnTe)
  if (Sd > 2000) Sd <- 2000; if (Sd < 0) Sd <- 0
  Srz <- Sr0; Suz <- 0; SWE <- 0
  Te  <- exp(lnTe); td_m <- td / m

  Qsim <- numeric(n); Qb_out <- numeric(n); Qof_out <- numeric(n)
  SWE_out <- numeric(n); ET_out <- numeric(n); Sd_out <- numeric(n)

  for (t in seq_len(n)) {
    tv <- tavg[t]
    if (tv <= snow_t) {
      SWE <- SWE + P[t]; rain <- 0
    } else {
      ma <- snow_mf * (tv - snow_t)
      if (ma > SWE) ma <- SWE
      SWE <- SWE - ma; rain <- P[t] + ma
    }
    SWE_out[t] <- SWE

    # Saturated area fraction — binary search on sorted thresholds
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

    # Saturation-excess overland flow
    qof   <- sf * rain
    infil <- rain - qof
    Qof_out[t] <- qof

    # Root zone
    Srz <- Srz - infil
    if (Srz < 0) { Suz <- Suz - Srz; Srz <- 0 }
    if (Srz < Srmax) { ae <- pet[t] * (1 - Srz / Srmax) } else { ae <- 0 }
    Srz <- Srz + ae; if (Srz > Srmax) Srz <- Srmax
    ET_out[t] <- ae

    # Unsaturated zone drainage
    if (Suz > 0 && Sd > 0) {
      quz <- Suz / (td_m * Sd); if (quz > Suz) quz <- Suz; Suz <- Suz - quz
    } else { quz <- 0 }

    # Update deficit
    Sd <- Sd - quz; if (Sd < 0) Sd <- 0

    # Baseflow
    qb <- Te * exp(-Sd / m)
    if (qb > Sd + 50) qb <- Sd + 50; if (qb < 0) qb <- 0
    Sd <- Sd + qb; if (Sd < 0) Sd <- 0

    Qb_out[t] <- qb
    Qsim[t]   <- qb + qof
    Sd_out[t]  <- Sd
  }

  Qsim[!is.finite(Qsim)] <- 0
  Qsim <- pmax(Qsim, 0)

  # NSE with 365-day warmup
  warmup <- min(365, n - 1)
  idx <- (warmup + 1):n
  o <- obs[idx]; s <- Qsim[idx]
  ok <- !is.na(o) & is.finite(s)
  if (sum(ok) < 30) {
    nse <- NA_real_
  } else {
    ss_res <- sum((o[ok] - s[ok])^2)
    ss_tot <- sum((o[ok] - mean(o[ok]))^2)
    nse <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
  }

  # Compute mean Sd per month for saturation map animation
  months_vec <- format(data$date, "%Y-%m")
  monthly_Sd <- tapply(Sd_out, months_vec, mean)
  monthly_labels <- names(monthly_Sd)

  list(Qsim = Qsim, Qobs = obs, Qb = Qb_out, Qof = Qof_out,
       SWE = SWE_out, ET = ET_out, dates = data$date,
       nse = nse, warmup = warmup,
       final_Sd = Sd, m = m, lam = lam,
       monthly_Sd = as.numeric(monthly_Sd),
       monthly_labels = monthly_labels)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4b. LIGHTWEIGHT NSE-ONLY TOPMODEL (optimizer only — no arrays)
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
  -(1 - ss_res / ss_obs)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5.  OPTIMISER — Nelder-Mead on 3-yr calibration window
# ═══════════════════════════════════════════════════════════════════════════════
optimise_params <- function(data, topidx, cal_years = 3) {

  par_names <- c("qs0","lnTe","m","Sr0","Srmax","td","snow_t","snow_mf")
  lower <- c(0.01, -7,   5,  0,   5,  1, -2, 1)
  upper <- c(10,    5, 100, 50, 300, 60,  2, 6)
  names(lower) <- par_names
  names(upper) <- par_names

  # Use last cal_years + 1yr warmup for calibration
  n_full   <- nrow(data)
  cal_days <- cal_years * 365 + 365
  if (n_full > cal_days) {
    cal_data <- data[(n_full - cal_days + 1):n_full, ]
  } else {
    cal_data <- data
  }

  # Precompute vectors once
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

  # 3 smart starts only — fast and effective
  starts <- list(
    c(qs0 = 1.0, lnTe =  1.0, m = 30, Sr0 =  2, Srmax = 100, td = 10, snow_t =  0.0, snow_mf = 3.0),
    c(qs0 = 2.0, lnTe =  2.0, m = 20, Sr0 =  1, Srmax =  60, td =  5, snow_t =  1.0, snow_mf = 2.0),
    c(qs0 = 3.0, lnTe = -2.0, m = 40, Sr0 = 10, Srmax = 150, td = 20, snow_t = -0.5, snow_mf = 3.5)
  )

  best_val  <- Inf
  best_pars <- starts[[1]]

  for (s in starts) {
    opt <- tryCatch(
      optim(s, obj_fast, method = "Nelder-Mead",
            control = list(maxit = 400, reltol = 1e-6)),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value) && opt$value < best_val) {
      best_val  <- opt$value
      best_pars <- opt$par
    }
  }

  names(best_pars) <- par_names
  pmax(pmin(best_pars, upper), lower)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6.  PRECOMPUTE AT STARTUP
# ═══════════════════════════════════════════════════════════════════════════════
message("Loading daily data ...")
ALL_DATA <- load_daily_data()
message(sprintf("  %d records: %s to %s",
                nrow(ALL_DATA), min(ALL_DATA$date), max(ALL_DATA$date)))

message("Computing TWI ...")
TWI_PRECOMPUTED <- file.path(DATA_DIR, "twi_precomputed.tif")
if (file.exists(TWI_PRECOMPUTED)) {
  # Fast path: load precomputed TWI (works on shinyapps.io without WhiteBox)
  message("  Loading precomputed TWI ...")
  TWI_RASTER <- rast(TWI_PRECOMPUTED)
} else {
  # Slow path: compute via WhiteBox (local only)
  if (!HAS_WHITEBOX) stop("WhiteBox not installed and no precomputed TWI found. Run locally first.")
  message("  Computing via WhiteBox Tools (FD8) ...")
  TWI_RASTER <- compute_twi_raster()
  # Save for next time / web deployment
  writeRaster(TWI_RASTER, TWI_PRECOMPUTED, overwrite = TRUE)
  message("  Saved precomputed TWI to ", TWI_PRECOMPUTED)
}
TOPIDX     <- make_topidx_classes(TWI_RASTER, n_classes = 16)
TWI_LAM    <- sum(TOPIDX[, 1] * TOPIDX[, 2])
message("  Lambda (mean TWI): ", round(TWI_LAM, 2))

WS_SF_NATIVE <- st_read(SHP_F, quiet = TRUE)
WS_SF_NATIVE <- st_transform(WS_SF_NATIVE, crs(TWI_RASTER))

message("  Startup complete.")

# Clipped DEM for contour lines on maps
DEM_WS <- mask(crop(rast(DEM_F),
                     vect(st_transform(WS_SF_NATIVE, crs(rast(DEM_F))))),
               vect(st_transform(WS_SF_NATIVE, crs(rast(DEM_F)))))

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: snap value to nearest slider step
# ═══════════════════════════════════════════════════════════════════════════════
snap_to_step <- function(val, min_val, max_val, step) {
  snapped <- round((val - min_val) / step) * step + min_val
  max(min(snapped, max_val), min_val)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7.  SHINY UI
# ═══════════════════════════════════════════════════════════════════════════════
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
  .pg .help-block { font-size: 11px; color: #888; margin-top: -4px;
                    margin-bottom: 10px; line-height: 1.3; }
"

ui <- page_sidebar(
  title = "HBEF Watershed 3 \u2014 TOPMODEL Explorer",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  tags$head(tags$style(HTML(app_css))),

  sidebar = sidebar(
    width = 330,

    sliderInput("year_range", "Simulation Period",
                min   = as.integer(format(min(ALL_DATA$date), "%Y")),
                max   = as.integer(format(max(ALL_DATA$date), "%Y")),
                value = c(2000, 2001),
                step  = 1, sep = "",
                ticks = TRUE),

    actionButton("run_model", "\u25B6  Run Model",
                 class = "btn-success action-btn"),

    uiOutput("nse_display"),
    uiOutput("opt_params_display"),

    div(class = "pg",
      h6("\u2193 Subsurface"),
      sliderInput("qs0",  "qs0 \u2013 init. flow (mm/d)",
                  min = 0.01, max = 10, value = 5.5, step = 0.01),
      helpText("Initial subsurface flow at time zero. Sets the starting saturation deficit."),
      sliderInput("lnTe", "ln(Te) \u2013 transmissivity",
                  min = -7, max = 5, value = 1.9, step = 0.1),
      helpText("Log of soil transmissivity when fully saturated. Controls baseflow magnitude \u2014 higher = more baseflow."),
      sliderInput("m",    "m \u2013 decay param (mm)",
                  min = 5, max = 100, value = 19, step = 1),
      helpText("Rate transmissivity declines with depth. Small m = flashy response; large m = slow drainage."),
      sliderInput("td",   "td \u2013 unsat. delay (days)",
                  min = 1, max = 60, value = 3, step = 1),
      helpText("Time delay for water to move through the unsaturated zone to the water table.")
    ),

    div(class = "pg",
      h6("\U0001F331 Root Zone"),
      sliderInput("Sr0",   "Sr0 \u2013 init. deficit (mm)",
                  min = 0, max = 50, value = 2.5, step = 0.5),
      helpText("Initial root zone storage deficit. How dry the soil is at the start of simulation."),
      sliderInput("Srmax", "Srmax \u2013 max capacity (mm)",
                  min = 5, max = 300, value = 65, step = 5),
      helpText("Maximum water the root zone can hold before draining. Controls ET \u2014 larger = more ET in summer.")
    ),

    div(class = "pg",
      h6("\u2744 Snow"),
      sliderInput("snow_t",  "Tmelt (\u00B0C)",
                  min = -2, max = 2, value = -0.7, step = 0.1),
      helpText("Temperature threshold: below this, precipitation falls as snow."),
      sliderInput("snow_mf", "Melt coeff (mm/\u00B0C/d)",
                  min = 1, max = 6, value = 3.4, step = 0.1),
      helpText("Degree-day melt factor. mm of snowmelt per degree above Tmelt per day.")
    ),

    actionButton("optimise", "\u26A1  Optimise Parameters",
                 class = "btn-primary action-btn"),
    helpText("Nelder-Mead on 3-yr window (~5\u201320 sec)")
  ),

  navset_card_tab(
    nav_panel("Hydrograph",      plotlyOutput("hydrograph", height = "520px")),
    nav_panel("Scatter Plot",    plotlyOutput("scatter",    height = "520px")),
    nav_panel("Flow Duration",   plotlyOutput("fdc",        height = "520px")),
    nav_panel("Snow & ET",       plotlyOutput("snow_et",    height = "520px")),
    nav_panel("TWI Map",         plotOutput("twi_map",      height = "560px")),
    nav_panel("Saturation Map",
              fluidRow(
                column(12,
                  sliderInput("sat_month", "Month",
                              min = 1, max = 12, value = 1,
                              step = 1, width = "100%",
                              animate = animationOptions(interval = 800, loop = TRUE)),
                  plotOutput("sat_map", height = "480px")
                )
              ))
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# 8.  SHINY SERVER
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  model_result   <- reactiveVal(NULL)
  opt_params_txt <- reactiveVal(NULL)

  sim_data <- reactive({
    req(input$year_range)
    start_date <- as.Date(paste0(input$year_range[1], "-01-01"))
    end_date   <- as.Date(paste0(input$year_range[2], "-12-31"))
    d <- ALL_DATA[ALL_DATA$date >= start_date & ALL_DATA$date <= end_date, ]
    validate(need(nrow(d) > 400, "Select a wider date range (> 1 year)."))
    d
  })

  current_pars <- reactive({
    c(qs0 = input$qs0, lnTe = input$lnTe, m = input$m,
      Sr0 = input$Sr0, Srmax = input$Srmax, td = input$td,
      snow_t = input$snow_t, snow_mf = input$snow_mf)
  })

  # ── RUN MODEL ──
  observeEvent(input$run_model, {
    d <- tryCatch(sim_data(), error = function(e) NULL)
    if (is.null(d)) {
      showNotification("Check date range.", type = "error"); return()
    }
    res <- tryCatch(
      run_topmodel(current_pars(), d, TOPIDX),
      error = function(e) {
        showNotification(paste("Error:", e$message), type = "error"); NULL
      }
    )
    if (!is.null(res)) {
      model_result(res)
      nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
      showNotification(sprintf("NSE = %s", nse_str), type = "message", duration = 3)
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)

  # ── AUTO-RUN at startup ──
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

  # ── OPTIMISE ──
  observeEvent(input$optimise, {
    showNotification("Optimising (5-yr calibration window) \u2026",
                     type = "message", duration = NULL, id = "opt_note")

    data <- tryCatch(sim_data(), error = function(e) NULL)
    if (is.null(data)) {
      removeNotification("opt_note")
      showNotification("Check date range.", type = "error"); return()
    }

    best <- optimise_params(data, TOPIDX, cal_years = 5)
    removeNotification("opt_note")

    # Snap values to slider steps
    updateSliderInput(session, "qs0",     value = snap_to_step(best["qs0"],     0.01, 10,  0.01))
    updateSliderInput(session, "lnTe",    value = snap_to_step(best["lnTe"],    -7,   5,   0.1))
    updateSliderInput(session, "m",       value = snap_to_step(best["m"],        5,   100,  1))
    updateSliderInput(session, "td",      value = snap_to_step(best["td"],       1,   60,   1))
    updateSliderInput(session, "Sr0",     value = snap_to_step(best["Sr0"],      0,   50,   0.5))
    updateSliderInput(session, "Srmax",   value = snap_to_step(best["Srmax"],    5,   300,  5))
    updateSliderInput(session, "snow_t",  value = snap_to_step(best["snow_t"],  -2,   2,    0.1))
    updateSliderInput(session, "snow_mf", value = snap_to_step(best["snow_mf"],  1,   6,    0.1))

    # Show exact optimised values
    opt_params_txt(sprintf(
      "qs0=%.3f  lnTe=%.2f  m=%.1f\nSr0=%.1f  Srmax=%.1f  td=%.1f\nTmelt=%.2f  Melt=%.2f",
      best["qs0"], best["lnTe"], best["m"],
      best["Sr0"], best["Srmax"], best["td"],
      best["snow_t"], best["snow_mf"]
    ))

    # Run full record with EXACT optimised params
    res <- tryCatch(run_topmodel(best, data, TOPIDX), error = function(e) NULL)
    if (!is.null(res)) {
      model_result(res)
      nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
      showNotification(sprintf("Optimised! NSE = %s (full period)", nse_str),
                       type = "message", duration = 5)
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # ── NSE DISPLAY ──
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

  # ── OPTIMISED PARAMS READOUT ──
  output$opt_params_display <- renderUI({
    txt <- opt_params_txt()
    if (is.null(txt)) return(NULL)
    div(class = "opt-params",
        tags$strong("Optimised values:"), tags$br(),
        tags$pre(style = "margin:2px 0 0 0; white-space:pre-wrap;", txt))
  })

  # ════════════════════════════════════════════════════════════════════════════
  # PLOTS
  # ════════════════════════════════════════════════════════════════════════════

  output$hydrograph <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    d <- data.frame(date = res$dates, Observed = res$Qobs, Simulated = res$Qsim)
    nse_str <- if (is.finite(res$nse)) sprintf("%.3f", res$nse) else "N/A"
    plot_ly(d, x = ~date) %>%
      add_lines(y = ~Observed,  name = "Observed",
                line = list(color = "#2980b9", width = 1.2)) %>%
      add_lines(y = ~Simulated, name = "Simulated",
                line = list(color = "#e74c3c", width = 1.2)) %>%
      layout(title = list(text = paste0("Daily Streamflow  |  NSE = ", nse_str),
                          font = list(size = 15)),
             xaxis = list(title = ""),
             yaxis = list(title = "Streamflow (mm/day)"),
             legend = list(orientation = "h", y = -0.12),
             hovermode = "x unified", margin = list(t = 50))
  })

  output$scatter <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    idx <- (res$warmup + 1):length(res$Qsim)
    d <- data.frame(obs = res$Qobs[idx], sim = res$Qsim[idx])
    d <- d[!is.na(d$obs) & is.finite(d$sim), ]
    validate(need(nrow(d) > 0, "No valid data."))
    rng <- range(c(d$obs, d$sim), na.rm = TRUE)
    plot_ly(d, x = ~obs, y = ~sim, type = "scattergl", mode = "markers",
            marker = list(size = 3, color = "#3498db", opacity = 0.4)) %>%
      add_lines(x = rng, y = rng, name = "1:1",
                line = list(color = "black", dash = "dash", width = 1)) %>%
      layout(title = list(text = "Observed vs Simulated", font = list(size = 15)),
             xaxis = list(title = "Observed (mm/d)"),
             yaxis = list(title = "Simulated (mm/d)"),
             showlegend = FALSE, margin = list(t = 50))
  })

  output$fdc <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    idx   <- (res$warmup + 1):length(res$Qsim)
    obs_v <- sort(res$Qobs[idx][!is.na(res$Qobs[idx])], decreasing = TRUE)
    sim_v <- sort(res$Qsim[idx][is.finite(res$Qsim[idx])], decreasing = TRUE)
    exc_o <- seq_along(obs_v) / length(obs_v) * 100
    exc_s <- seq_along(sim_v) / length(sim_v) * 100
    plot_ly() %>%
      add_lines(x = exc_o, y = obs_v, name = "Observed",
                line = list(color = "#2980b9", width = 1.5)) %>%
      add_lines(x = exc_s, y = sim_v, name = "Simulated",
                line = list(color = "#e74c3c", width = 1.5)) %>%
      layout(title = list(text = "Flow-Duration Curve", font = list(size = 15)),
             xaxis = list(title = "Exceedance (%)"),
             yaxis = list(title = "Streamflow (mm/d)", type = "log"),
             legend = list(orientation = "h", y = -0.12), margin = list(t = 50))
  })

  output$snow_et <- renderPlotly({
    res <- model_result()
    validate(need(!is.null(res), "Press \u25B6 Run Model to see results."))
    d <- data.frame(date = res$dates, SWE = res$SWE, ET = res$ET)
    plot_ly(d, x = ~date) %>%
      add_lines(y = ~SWE, name = "SWE (mm)",
                line = list(color = "#3498db", width = 1.3),
                yaxis = "y") %>%
      add_lines(y = ~ET,  name = "Actual ET (mm/d)",
                line = list(color = "#27ae60", width = 1.2),
                yaxis = "y2") %>%
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

  # ── TWI MAP — shows the underlying topographic wetness index ──
  output$twi_map <- renderPlot({
    par(mar = c(2, 2, 3, 4), bg = "white")
    plot(TWI_RASTER,
         col  = hcl.colors(50, "YlGnBu", rev = TRUE),
         axes = FALSE,
         main = sprintf("Topographic Wetness Index  (FD8, \u03BB = %.1f)", TWI_LAM),
         cex.main = 1.3)
    plot(st_geometry(WS_SF_NATIVE), add = TRUE,
         border = "#222222", lwd = 2.5, col = NA)

    # Elevation contour lines
    contour(raster::raster(DEM_WS), add = TRUE,
            nlevels = 10, col = adjustcolor("#333333", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)
  }, res = 120)

  # ── SATURATION MAP — binary saturated/unsaturated from model output ──
  # ── Update month slider range when model results change ──
  observe({
    res <- model_result()
    req(!is.null(res))
    n_months <- length(res$monthly_Sd)
    updateSliderInput(session, "sat_month",
                      min = 1, max = n_months, value = n_months)
  })

  output$sat_map <- renderPlot({
    res <- model_result()
    validate(need(!is.null(res), "Run or optimise the model first."))

    month_idx <- input$sat_month
    req(month_idx)
    n_months <- length(res$monthly_Sd)
    if (month_idx > n_months) month_idx <- n_months

    Sd  <- res$monthly_Sd[month_idx]
    m   <- res$m
    lam <- res$lam
    thresh <- lam + Sd / m
    label <- res$monthly_labels[month_idx]

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
         col    = c("#f5e6c8", "#2166ac"),
         breaks = c(-0.5, 0.5, 1.5),
         legend = FALSE, axes = FALSE,
         main   = sprintf("%s  \u2014  %.1f%% saturated", label, pct_sat),
         cex.main = 1.3)
    plot(st_geometry(WS_SF_NATIVE), add = TRUE,
         border = "#222222", lwd = 2.5, col = NA)

    # Elevation contour lines for terrain context
    contour(raster::raster(DEM_WS), add = TRUE,
            nlevels = 10, col = adjustcolor("#555555", alpha.f = 0.5),
            lwd = 0.6, drawlabels = TRUE, labcex = 0.6)

    legend("bottomright",
           legend = c("Saturated", "Unsaturated"),
           fill   = c("#2166ac", "#f5e6c8"),
           border = c("#2166ac", "#f5e6c8"),
           bty = "n", cex = 1.1)
  }, res = 120)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9.  LAUNCH
# ═══════════════════════════════════════════════════════════════════════════════
shinyApp(ui, server)

