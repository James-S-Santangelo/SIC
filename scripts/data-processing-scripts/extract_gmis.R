library(raster)

# Load population-mean dataframe and retrieve lat/longs
latLongs <- read_csv("data-clean/experimentalData_popMeans.csv") %>% 
  dplyr::select(Population, Latitude, Longitude)

# Load GMIS dataset
imperv <- raster("gis/gmis/gmis_impervious_surface_percentage_geographic_30m.tif")

# Create spatial-points dataframe
spImperv <- SpatialPointsDataFrame(latLongs[,3:2], proj4string=imperv@crs, latLongs)

# Extract GMIS value for each population
# If gmis == 200, set to 0, as per GMIS README
percent_imperv <- raster::extract(imperv, spImperv, df=TRUE) %>% 
  mutate(gmis = ifelse(gmis_impervious_surface_percentage_geographic_30m == 200, 0, gmis_impervious_surface_percentage_geographic_30m)) %>% 
  dplyr::select(-gmis_impervious_surface_percentage_geographic_30m) %>% 
  rename("Population" = "ID") %>% 
  mutate(Population = as.character(Population))

write_csv(percent_imperv, "data-clean/gmis.csv")

