library(rgdal)
library(raster)
library(geosphere)

# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)


# Read shapefiles of river bank and houses for sale
river_bank <-
  readOGR("./Inputs/Susquehanna_west_bank/west_bank.shp") #river bank
houses_for_sale <-
  readOGR("./Inputs/houses_for_sale/houses_for_sale.shp") #

# Calculate distance from river and elevation above river for houses
points <-
  as.matrix(data.frame(houses_for_sale$X, houses_for_sale$Y))
line <- spTransform(river_bank, crs(houses_for_sale))
distance_from_river <- dist2Line(points, line)
elevation_above_river <-
  houses_for_sale$Z - 126 #average river bank elevation = 126 m

# Multivariate linear regression to estimate unit house price based on elevation and distance from river
unit_house_price <- as.numeric(houses_for_sale$Price) /
  (as.numeric(houses_for_sale$Living_Are) * 0.3048 ^ 2) #unit house price in USD/sq.m
data <-
  data.frame(distance_from_river, elevation_above_river, unit_house_price)
fit <- lm(unit_house_price ~ distance + elevation_above_river, data)
intercept <-
  fit$coefficients[1] # unit_house_price=c_dist*distance+c_elev*elevation+intercept
c_dist <- fit$coefficients[2]
c_elev <- fit$coefficients[3]

# Estimate unit price of hypothetical houses
houses <- raster("./Inputs/houses.asc")
houses_points <- as.data.frame(houses, xy = T)
utmcoor <-
  SpatialPoints(cbind(houses_points$x, houses_points$y),
                proj4string = CRS("+proj=utm +zone=18"))
longlatcoor <- spTransform(utmcoor, CRS("+proj=longlat"))
houses_points$x <- coordinates(longlatcoor)[, 1]
houses_points$y <- coordinates(longlatcoor)[, 2]
houses_points <- houses_points[complete.cases(houses_points), ]
xy <- houses_points[, 1:2]
house_distance <- as.data.frame(dist2Line(xy, line))
xy <- as.data.frame(houses, xy = T)
xy <- xy[complete.cases(xy), ]
xyz = cbind(xy[, 1:2], z = house_distance$distance)
distance_raster <- rasterFromXYZ(xyz)

elevation_raster <- houses - 126
values(elevation_raster)[values(elevation_raster) < 0] = 0

house_price = intercept + (c_dist * distance_raster) + (c_elev * elevation_raster)
writeRaster(house_price,
            "./Inputs/house_price.asc",
            format = "ascii",
            overwrite = T)
