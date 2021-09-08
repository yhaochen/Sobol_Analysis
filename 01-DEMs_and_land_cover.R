library(FedData)
library(rgdal)
library(raster)

# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)

# Read Selinsgrove's boundaries shapefile
clip_area <- readOGR("./Inputs/clip_area/clip_area.shp")

# Download and save USGS 10m DEM for Selinsgrove
dem <- get_ned(template = clip_area,
               res = '13',
               label = 'DEM')
crs <- crs(clip_area)
DEM10 <- projectRaster(dem, res = 10, crs = crs)
DEM10 <- round(DEM10, digits = 2)
DEM30 <- aggregate(DEM10, fact = 30 / 10, fun = mean)
DEM30 <- round(DEM30, digits = 2)
DEM50 <- aggregate(DEM10, fact = 50 / 10, fun = mean)
DEM50 <- round(DEM50, digits = 2)

# Clipped land use land cover map for Selinsgrove
LULC <- raster("./Inputs/NLCD/NLCD2016_clip")
LULC_proj <- projectRaster(LULC,
                           res = 30,
                           crs = crs,
                           method = "ngb")

writeRaster(DEM10,
            "./LISFLOOD/dem10.asc",
            format = "ascii",
            overwrite = T)
writeRaster(DEM30,
            './LISFLOOD/dem30.asc',
            format = "ascii",
            overwrite = T)
writeRaster(DEM50,
            './LISFLOOD/dem50.asc',
            format = "ascii",
            overwrite = T)
writeRaster(LULC_proj,
            "./Inputs/lulc.asc",
            format = "ascii",
            overwrite = T)
