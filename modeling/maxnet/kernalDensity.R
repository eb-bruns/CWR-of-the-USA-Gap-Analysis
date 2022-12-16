###
# Function to generate a kernal density map of the occurence data for CWR NA.
# This is a visualization tool that will help in understanding the potential
# spatial biasing in the occurrence dataset. I want to be able to highlight
# areas of high sampling density and areas of low density. A possible extension
# of this it attempting to smooth the varability acroos the area,
# by adjusting the number of points in a region, though that will be a
# complicated process.
# 20200414
# dan.carver@carverd.com
###

kernalDensity <- function(species){
  # EBB: new spatialEco (v2.0) has sf.kde for kernal density, not sp.kde,
  # so I've updated to sf and we'll see if everything else runs as expected
  #cleanPoints_sf <- st_as_sf(cleanPoints)
  #thrshold_rast <- rast(thrshold)
  k2 <- spatialEco::sp.kde(x = cleanPoints, newdata = thrshold, standardize = TRUE)
  raster::writeRaster(x = k2, filename = paste0(sp_dir, "/modeling/kernalDensity.tif"),overwrite=TRUE)

  # potential work flow for test idea
  # extract kernal values back to cleanPoints
  # filter out a portion of the points with high kernal density
  # re run KDE
  # repeat until a range the distribution of the kernal density values reachs
  # a specific level of normal? or some other metric(range of KDE values)
  # rerun the modeling methodology and compare outputs 
  }
