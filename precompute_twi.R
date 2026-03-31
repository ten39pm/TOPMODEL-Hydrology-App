# ═══════════════════════════════════════════════════════════════════════════════
# precompute_twi.R
#
# PURPOSE: This file pre-calculates the Topographic Wetness Index (TWI) from a
# raw Digital Elevation Model (DEM) and saves the result as a raster file. This
# step only needs to be run ONCE locally. It creates data/twi_precomputed.tif 
# so that the app doesn't run WhiteboxTools (WBT) every time.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Libraries ──────────────────────────────────────────────────────────────────
library(terra)     # package for raster & vector spatial data
library(sf)        # package for vector/shapefile data (watershed boundary)
library(whitebox)  # package for geospatial data analysis

# Initialize WBT
wbt_init()

# ── File paths ─────────────────────────────────────────────────────────────────
# Input: 10-meter resolution DEM for the Hubbard Brook Watershed 3
DEM_F <- "Data/dem_10m.tif"

# Input: watershed boundary polygon
SHP_F <- "Data/Watershed3HB.shp"

# Output: final pre-computed TWI raster
OUT_F <- "Data/twi_precomputed.tif"

# Temporary working directory for files produced by WBT
TEMP_DIR <- file.path(tempdir(), "twi_precompute")
dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Step 1: Load and clip DEM to watershed boundary ───────────────────────────
message("Loading DEM and watershed ...")

# Read the DEM raster into R as a SpatRaster object
dem <- rast(DEM_F)

# Read the watershed shapefile as an sf (simple features) object
ws  <- st_read(SHP_F, quiet = TRUE)

# Re-project the watershed polygon to match the DEM's coordinate reference
# system (CRS)
ws_v <- vect(st_transform(ws, crs(dem)))

# Clip (crop and mask) the DEM to the watershed boundary
dem_ws <- mask(crop(dem, ws_v), ws_v)

# Define file paths for the rasters that WBT will read/write
dem_clip_f <- file.path(TEMP_DIR, "dem_clip.tif")       # clipped DEM
breached_f <- file.path(TEMP_DIR, "dem_breached.tif")   # breached DEM
sca_f      <- file.path(TEMP_DIR, "sca.tif")      # specific contributing area
slope_f    <- file.path(TEMP_DIR, "slope.tif")    # slope in degrees

# Write the clipped DEM to file for WBT
writeRaster(dem_ws, dem_clip_f, overwrite = TRUE)

# ── Step 2: Running WBT ─────────────────────────────────────────────────
message("Breaching depressions ...")
wbt_breach_depressions_least_cost(
  dem = dem_clip_f, 
  output = breached_f, 
  dist = 10, 
  fill = TRUE
)

message("FD8 flow accumulation ...")
wbt_fd8_flow_accumulation(
  dem = breached_f, 
  output = sca_f, 
  out_type = "specific contributing area"
)

message("Computing slope ...")
wbt_slope(dem = breached_f, 
          output = slope_f, 
          units = "degrees")

message("Computing TWI ...")

# load the WBT output rasters back into R
sca <- rast(sca_f)
slp <- rast(slope_f)

# Convert slope from degrees to radians
slp_rad <- slp * pi / 180

# Compute the slope gradient
tan_b <- tan(slp_rad)

# Prevent infinite TWI on flat areas
tan_b[tan_b < 0.001] <- 0.001

# TWI formula
twi <- log(sca / tan_b)

# Mask the TWI to the watershed boundary
twi <- mask(twi, dem_ws)

# ── Step 3: Save output ────────────────────────────────────────────────────────
# Write the final TWI raster to the Data/ folder.
writeRaster(twi, OUT_F, overwrite = TRUE)

message("Done! Saved to: ", OUT_F)
message("You can now deploy to shinyapps.io without WhiteBox.")

# ── Step 4: Create TOPIDX Function ─────────────────────────────────────────────
# Function that turns TWI values into histogram-style classes for TOPMODEL's
# topographic index look-up table (TOPIDX)
make_topidx_classes <- function(twi_raster, n_classes = 16) {
  vals <- values(twi_raster, na.rm = TRUE)[, 1]
  brks <- seq(min(vals), max(vals), length.out = n_classes + 1)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  cnts <- hist(vals, breaks = brks, plot = FALSE)$counts
  frac <- cnts / sum(cnts)
  cbind(twi = mids, frac = frac)
}
