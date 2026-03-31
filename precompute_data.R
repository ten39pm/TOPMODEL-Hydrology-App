# ═══════════════════════════════════════════════════════════════════════════════
# precompute_data.R
#
# PURPOSE: This file should be run ONCE locally before deploying to shinyapps.io.
# It saves all the heavy pre-processing results to Data/ so that the app loads 
# them instantly at start-up without re-computing anything.
# ═══════════════════════════════════════════════════════════════════════════════
  
# ── Libraries ──────────────────────────────────────────────────────────────────
library(terra)      # package for raster & vector spatial data
library(sf)         # package for vector/shapefile data (watershed boundary)
library(dplyr)      # data wrangling
library(lubridate)  # data handling
library(whitebox)   # package for geospatial data analysis

# Initialize WBT
wbt_init()

# File paths
DATA_DIR <- "data"
PRECIP_F <- file.path(DATA_DIR, "DailyWatershed.csv")
FLOW_F   <- file.path(DATA_DIR, "HBEF_DailyStreamflow_1956-2024.csv")
TEMP_F   <- file.path(DATA_DIR, "HBEF_air_temp_daily.csv")
DEM_F    <- file.path(DATA_DIR, "dem_10m.tif")
SHP_F    <- file.path(DATA_DIR, "Watershed3HB.shp")
TWI_F    <- file.path(DATA_DIR, "twi_precomputed.tif")  # from precompute_twi.R
OUT_F    <- file.path(DATA_DIR, "app_startup.rds")      # output of this script

# ── Step 1: Load & merge daily climate/flow data ─────────────────────────────
message("Loading daily data ...")

# Read & filter each dataset to Watershed 3
precip <- read.csv(PRECIP_F, stringsAsFactors = FALSE)
precip$DATE <- as.Date(trimws(precip$DATE))
precip <- precip[trimws(precip$watershed) == "W3", c("DATE", "Precip")]
names(precip) <- c("date", "precip_mm")

flow <- read.csv(FLOW_F, stringsAsFactors = FALSE)
flow$DATE <- as.Date(trimws(flow$DATE))
flow <- flow[flow$WS == 3, c("DATE", "Streamflow")]
names(flow) <- c("date", "flow_mm")

# Average min/max temperatures to get the daily mean
temp <- read.csv(TEMP_F, stringsAsFactors = FALSE)
temp$date <- as.Date(trimws(temp$date))
temp_avg <- temp |>
  group_by(date) |>
  summarise(tavg = mean(AVE, na.rm = TRUE), .groups = "drop")

# Inner join precip & flow and left join by temp_avg
df <- merge(precip, flow, by = "date", all = FALSE)
df <- merge(df, temp_avg, by = "date", all.x = TRUE)
df <- df[order(df$date), ]

# Fill missing temperature values by linear interpolation b/w valid points
if (any(is.na(df$tavg))) {
  idx <- which(!is.na(df$tavg))
  df$tavg <- approx(idx, df$tavg[idx], seq_len(nrow(df)), rule = 2)$y
}
df$precip_mm[is.na(df$precip_mm)] <- 0
ALL_DATA <- df

message(sprintf("  %d records: %s to %s",
                nrow(ALL_DATA), min(ALL_DATA$date), max(ALL_DATA$date)))

# ── Step 2: Load precomputed TWI (run precompute_twi.R first if missing) ───────
message("Loading TWI raster ...")
if (!file.exists(TWI_F)) {
  stop("TWI raster not found. Please run precompute_twi.R first, then re-run this script.")
}
TWI_RASTER <- rast(TWI_F)

# ── Step 3: Build TOPMODEL topographic index table ────────────────────────────
# Bins TWI values into 16 classes; each class gets a midpoint TWI value
# and the fraction of watershed area it represents.
message("Building topographic index table ...")
vals <- values(TWI_RASTER, na.rm = TRUE)[, 1]
brks <- seq(min(vals), max(vals), length.out = 17)  # 16 bins = 17 break points
mids <- (brks[-1] + brks[-length(brks)]) / 2
cnts <- hist(vals, breaks = brks, plot = FALSE)$counts
frac <- cnts / sum(cnts)
TOPIDX  <- cbind(twi = mids, frac = frac)
TWI_LAM <- sum(TOPIDX[, 1] * TOPIDX[, 2])           # area-weighted mean TWI (lambda)

message("  Lambda (mean TWI): ", round(TWI_LAM, 2))

# ── Step 4: Load watershed boundary & re-project to match the TWI raster CRS ─────────────────────────────
message("Loading watershed shapefile ...")
WS_SF_NATIVE <- st_read(SHP_F, quiet = TRUE)
WS_SF_NATIVE <- st_transform(WS_SF_NATIVE, crs(TWI_RASTER))

# ── Step 5: Clip DEM to watershed boundary for drawing elevation contours on maps ───────────────────
message("Clipping DEM to watershed ...")
ws_vect <- vect(st_transform(WS_SF_NATIVE, crs(rast(DEM_F))))
DEM_WS  <- mask(crop(rast(DEM_F), ws_vect), ws_vect)

# ── Step 6: Bundle everything and save to a single .rds file ──────────────────
# The app loads this one file at start-up instead of re-running all of the above.
message("Saving startup cache to ", OUT_F, " ...")
startup_cache <- list(
  ALL_DATA     = ALL_DATA,      # merged daily precip/flow/temp data frame
  TOPIDX       = TOPIDX,        # TWI class table (midpoints + area fractions)
  TWI_LAM      = TWI_LAM,       # mean TWI (lambda)
  WS_SF_NATIVE = WS_SF_NATIVE,  # watershed boundary in TWI CRS
  DEM_WS       = DEM_WS         # watershed-clipped DEM for contour plotting
)
saveRDS(startup_cache, OUT_F)
message("Done! You can now deploy to shinyapps.io.")
message("Files needed in Data/: twi_precomputed.tif, app_startup.rds")