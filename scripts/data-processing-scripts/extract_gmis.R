# Script to extract impervious surface data (GMIS) from rasters 
# NOTE: Rasters are not provided in repository but can be downloaded from:
# https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1

# Load required packages
library(tidyverse)
library(sp)
library(raster)

# Load Lat/Longs
latlongs <- read_csv('data-clean/experimentalData_popMeans.csv') %>% 
  dplyr::select(Population, Latitude, Longitude)

# Load GMIS raster datasets
gmis1000_raster <- raster('gis/gmis/17T_gmis_impervious_surface_percentage_geographic_1000m.tif')
gmis30_raster <- raster('gis/gmis/17T_gmis_impervious_surface_percentage_geographic_30m.tif')

# Create spatial points dataframes
spdf <- SpatialPointsDataFrame(coords = latlongs %>% dplyr::select(Longitude, Latitude), proj4string = gmis1000_raster@crs, data = latlongs)

# Extract GMIS data for populations
gmis1000_data <- raster::extract(x = gmis1000_raster, y = spdf, method = 'simple')
gmis30_data <- raster::extract(x = gmis30_raster, y = spdf, method = 'simple')

# Create df and write to disk
gmis_out <- latlongs %>% 
  mutate(gmis1000 = gmis1000_data,
         gmis30 = gmis30_data,
         gmis30 = ifelse(gmis30 == 200, 0, gmis30))

write_csv(gmis_out, 'data-clean/gmis.csv')
