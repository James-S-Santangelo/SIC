library(raster)

# Load population-mean dataframe and retrieve lat/longs
latLongs <- read_csv("data-clean/experimentalData_popMeans.csv") %>% 
  dplyr::select(Population, Latitude, Longitude, Distance_to_core) %>% 
  mutate(Population = as.character(Population))

# Load GMIS dataset
gmis30 <- raster("gis/gmis/17T_gmis_impervious_surface_percentage_geographic_30m.tif")
gmis250 <- raster("gis/gmis/17T_gmis_impervious_surface_percentage_geographic_250m.tif")
gmis1000 <- raster("gis/gmis/17T_gmis_impervious_surface_percentage_geographic_1000m.tif")

# Create spatial-points dataframe
spImperv <- SpatialPointsDataFrame(latLongs[,3:2], proj4string=gmis30@crs, latLongs)

# Extract GMIS value for each population
# If gmis == 200, set to 0, as per GMIS README
gmis_df <- raster::extract(gmis30, spImperv, df=TRUE) %>% 
  left_join(., raster::extract(gmis250, spImperv, df=TRUE), by = "ID") %>% 
  left_join(., raster::extract(gmis1000, spImperv, df=TRUE), by = "ID") %>% 
  rename("gmis30" = "X17T_gmis_impervious_surface_percentage_geographic_30m",
         "gmis250" = "X17T_gmis_impervious_surface_percentage_geographic_250m",
         "gmis1000" = "X17T_gmis_impervious_surface_percentage_geographic_1000m") %>% 
  mutate(gmis30 = ifelse(gmis30 == 200, 0, gmis30)) %>% 
  rename("Population" = "ID") %>% 
  mutate(Population = as.character(Population)) %>% 
  left_join(., latLongs %>% dplyr::select(Population, Distance_to_core), by = "Population") %>% 
  gather(., key = "gmis", value = "value", gmis30:gmis1000)

write_csv(gmis_df, "data-clean/gmis.csv")

