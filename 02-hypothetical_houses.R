library(rgdal)
library(raster)

# Set working directory
wd <- getwd()
setwd(wd)

# Read land use grid and extract urban regions
LULC <- raster("./Inputs/lulc.asc")
roads <- readOGR("./Inputs/Roads/roads.shp") # Selinsgrove roads
rivers <- readOGR("./Inputs/Rivers/rivers.shp")
LULC_no_roads <- mask(LULC, roads, inverse = T)
LULC_no_roads_no_rivers <- mask(LULC_no_roads, rivers, inverse = T)
urban <- function(x) {
  ifelse (x == 22 | x == 23 | x == 24, 1, NA)
}
urban_region <- calc(LULC_no_roads_no_rivers, fun = urban)
selinsgrove <-
  readOGR("./Inputs/selinsgrove_shapefile/Selinsgrove.shp")
urban_extract = mask(urban_region, selinsgrove)

# Extract urban regions of DEM
dem <- raster("./LISFLOOD/dem10.asc")
shp <- rasterToPolygons(urban_region)
dem_urban <- mask(dem, shp)
dem_urban_selinsgrove <- mask(dem_urban, selinsgrove)
# Select 2000 hypothetical houses (grid cells) in the urban region
set.seed(1)
houses <- sampleRandom(dem_urban_selinsgrove, size = 2000, asRaster = TRUE)

writeRaster(houses,
            "./Inputs/houses.asc",
            format = "ascii",
            overwrite = TRUE)
writeRaster(
  dem_urban_selinsgrove,
  "./Inputs/dem_urban.asc",
  format = "ascii",
  overwrite = TRUE
)
