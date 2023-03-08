library(terra)
library(leastcostpath)
library(tmap)

#### LOAD AND PROCESS DATA ####

ita <- terra::rast("./Data/DEM/ITA_DEM_250m.tif")

ext <- sf::st_read("./Data/Outline/gadm41_ITA_1.shp")
ext <- ext[!ext$NAME_1 %in% c("Sicily", "Sardegna"),]
ext <- sf::st_transform(ext, sf::st_crs(ita))

ita2 <- terra::mask(ita, ext)
ita2[ita2 < 0] <- NA

terra::writeRaster(ita2, "./Data/DEM/ITA_DEM_250m_mask.tif", overwrite = TRUE)

roads <- sf::st_read("./Data/RR_Italy/Roman_Roads_Italy.shp")
roads <- sf::st_intersection(roads, ext)

sf::st_write(roads, "./Data/RR_Italy/Roman_Roads_Italy2.shp", overwite = TRUE)

RR_sites <- read.csv("./Data/RS_Italy/Database_ROMURBITAL.csv")
RR_sites <- RR_sites[RR_sites$Geog_certainty %in% c("certain", "likely") & RR_sites$Occ_Date_Cert %in% c("High", "Medium") & !RR_sites$Occ_Date_U %in% c("-9999", "-") & !RR_sites$Occ_Date_L %in% c("999", "-"),]
RR_sites$Occ_Date_U <- as.numeric(RR_sites$Occ_Date_U)
RR_sites$Occ_Date_L <- as.numeric(RR_sites$Occ_Date_L)
RR_sites$Occ_Date_diff <- abs(RR_sites$Occ_Date_U - RR_sites$Occ_Date_L)
RR_sites <- sf::st_as_sf(RR_sites, coords = c("X", "Y"), crs = "wgs84")
RR_sites <- sf::st_transform(RR_sites, sf::st_crs(roads))

sf::st_write(RR_sites, "./Data/RS_italy/Database_ROMURBITAL2.shp", overwrite = TRUE)

#### CREATE PRE-ROAD COST SURFACE ####

slope_RI_cs <- leastcostpath::create_slope_cs(x = ita, cost_function = "tobler", neighbours = 16)

slope_RI_rast <- leastcostpath::rasterise(slope_RI_cs)
slope_RI_rast <- terra::mask(slope_RI_rast, ext)
slope_RI_rast[is.infinite(slope_RI_rast)] <- NA
terra::writeRaster(slope_RI_rast, "./Output/cs1.tif", overwrite = TRUE)

#### UPDATE PRE-ROAD COST SURFACE TO INCORPORATE ESTIMATED SPEEDS WHEN USING ROADS ####

# assuming that the roads were used for vehicular transportation, the estimated speed is 36km a day / 3.6km an hour (assuming 10 hour workday). We need to convert to miles to metres per second to make comparable to Tobler's hiking function values. Need to divide by distance between cells to account for distance between cells

slope_RI_cs1 <- slope_RI_cs
slope_RI_cs1 <- leastcostpath::update_values(x = slope_RI_cs1, sf = roads, 
                                          FUN = function(x) { replace(x = x, values = (3.6/2.237)/max(terra::res(ita)))})

slope_RI_rast1 <- leastcostpath::rasterise(slope_RI_cs1)
slope_RI_rast1 <- terra::mask(slope_RI_rast1, ext)
slope_RI_rast1[is.infinite(slope_RI_rast1)] <- NA
terra::writeRaster(slope_RI_rast1, "./Output/cs2.tif", overwrite = TRUE)

#### CREATE ACCUMULATED COST SURFACES FROM ALL SETTLEMENTS ####

accum1 <- ita
accum2 <- ita

accum1[] <- 0
accum2[] <- 0

for(i in 1:nrow(RR_sites)) {
  
  print(i)
  accum_1 <- leastcostpath::create_accum_cost(x = slope_RI_cs, origin = RR_sites[i,])
  accum1 <- accum1  + (accum_1 / 3600)
  
  accum_2 <- leastcostpath::create_accum_cost(x = slope_RI_cs1, origin = RR_sites[i,])
  accum2 <- accum2  + (accum_2 / 3600)

}

accum1 <- accum1 / nrow(RR_sites)
accum2 <- accum2 / nrow(RR_sites)

accum1 <- terra::mask(accum1, ext)
accum2 <- terra::mask(accum2, ext)

accum1[is.infinite(accum1)] <- NA
accum2[is.infinite(accum2)] <- NA

terra::writeRaster(accum1, "./Output/accum1.tif")
terra::writeRaster(accum2, "./Output/accum2.tif")

accum3 <- abs((accum2-accum1) / (accum2) * 100)

terra::writeRaster(accum3, "./Output/accum3.tif")

#### CALCULATE METHOD SCHEMATIC PLOTS ####

normalise_raster <- function(x){(x-terra::minmax(x)[1])/(terra::minmax(x)[2]-terra::minmax(x)[1])}

dem <- terra::rast("./Data/DEM/FS_240.tif")
dem <- normalise_raster(dem) * 10

ext<- c(-500, -495, 495, 500)
dem2 <- terra::crop(dem, ext)

slope_cs1 <- leastcostpath::create_slope_cs(x = dem2, cost_function = "tobler", neighbours = 4)

pt <- sf::st_as_sf(data.frame(x = -499.5, y = 495.5), coords = c("x", "y"), crs = terra::crs(dem))

plot1 <- tm_shape(dem2) + 
  tm_raster(palette = grey.colors(10), style = "cont", title = "Elevation (m)") + 
  tm_shape(pt) + 
  tm_dots(col = "black", size = 2) + 
  tm_layout(legend.outside = TRUE, main.title = "A") + 
  tm_add_legend(type = "symbol", col = "black", labels = "Origin")

tmap::tmap_save(plot1, "./Output/Figures/plot1.png", dpi = 300, width = 3, height = 3)

accum_cs1 <- leastcostpath::create_accum_cost(x = slope_cs1, origin = pt)

plot2 <- tm_shape(accum_cs1) + 
  tm_raster(palette = grey.colors(10), style = "cont", title = "Accumulated\ncost (seconds)") + 
  tm_shape(pt) + 
  tm_dots(col = "black", size = 2) + 
  tm_layout(legend.outside = TRUE, main.title = "B")  + 
  tm_add_legend(type = "symbol", col = "black", labels = "Origin")

tmap::tmap_save(plot2, "./Output/Figures/plot2.png", dpi = 300, width = 3, height = 3)

ext2 <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(terra::ext(dem2))))
sf::st_crs(ext2) <- terra::crs(dem)

pt1 <- sf::st_as_sf(data.frame(x = -499.5, y = 499.5), coords = c("x", "y"), crs = terra::crs(dem))
pt2 <- sf::st_as_sf(data.frame(x = -495.5, y = 498.5), coords = c("x", "y"), crs = terra::crs(dem))

pt3 <- rbind(pt1, pt2)

lcps1 <- leastcostpath::create_lcp(x = slope_cs1, origin = pt, destination = pt3, cost_distance = TRUE)

plot3 <- tm_shape(dem2) + 
  tm_raster(palette = grey.colors(10), style = "cont", title = "Elevation (m)") + 
  tm_shape(lcps1) + 
  tm_lines(col = "black") + 
  tm_shape(pt) + 
  tm_dots(col = "black", size = 2) + 
  tm_layout(legend.outside = TRUE, main.title = "C") + 
  tm_add_legend(type = "symbol", col = "black", labels = "Origin") + 
  tm_add_legend(type = "line", col = "black", labels = "Road")

tmap::tmap_save(plot3, "./Output/Figures/plot3.png", dpi = 300, width = 3, height = 3)

slope_cs2 <- slope_cs1

# convert 49km per hour to metres per second by dividing by 3.6. Divided by maximum resolution of DEM (1000m) to make units comparable with speed calculated using leastcostpath (leastcostpath accounts for distance between cells when calculating speed)
slope_cs2 <- leastcostpath::update_values(x = slope_cs2, sf = lcps1, 
                                          FUN = function(x) { replace(x = x, values = (49/3.6)/max(terra::res(dem2)))})

accum_cs2 <- leastcostpath::create_accum_cost(x = slope_cs2, origin = pt)

plot4 <- tm_shape(accum_cs2) + 
  tm_raster(palette = grey.colors(10), style = "cont", title = "Accumulated\ncost (seconds)") + 
  tm_shape(pt) + 
  tm_dots(col = "black", size = 2) + 
  tm_layout(legend.outside = TRUE, main.title = "D") + 
  tm_add_legend(type = "symbol", col = "black", labels = "Origin")

tmap::tmap_save(plot4, "./Output/Figures/plot4.png", dpi = 300, width = 3, height = 3)

lcps2 <- leastcostpath::create_lcp(x = slope_cs2, origin = pt, destination = pt3, cost_distance = TRUE)
lcps2$cost_distance <- round(lcps2$cost_distance, 2)
lcps2$cost_distance <- factor(lcps2$cost_distance)

plot5 <- tm_shape(dem2) + 
  tm_raster(col = "white") +
  tm_shape(lcps2) + 
  tm_lines(col = "cost_distance", title.col = "Accumulated\ncost (seconds)", palette = grey.colors(2), lwd = 3) + 
  tm_shape(pt) + 
  tm_dots(col = "black", size = 2) + 
  tm_layout(legend.outside = TRUE, main.title = "E") + 
  tm_add_legend(type = "symbol", col = "black", labels = "Origin") + 
  tm_add_legend(type = "line", col = "black", labels = "Road")

tmap::tmap_save(plot5, "./Output/Figures/plot5.png", dpi = 300, width = 3, height = 3)












