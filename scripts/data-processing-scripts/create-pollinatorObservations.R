library(tidyverse)

# Load in pollinator observation data
pollinatorObservations <- read_csv("data-raw/pollinatorObservations.csv")

# Remove rows where wind ruined pollinator observations
pollinatorObservations <- pollinatorObservations %>%
  filter(!grepl("Breeze|Wind",Comments)) %>% 
  mutate(Population = as.character(Population))

# Load in data with lat/longs
latLongs <- read_csv("data-raw/populationLatLongs.csv") %>% 
  mutate(Population = as.character(Population))

# Add population lat/long data to seed/flower ratio dataframe
pollinatorObservations <- pollinatorObservations %>%
  left_join(., latLongs, by = "Population")

# Add impervious survafe values to data 
pollinatorObservations <- pollinatorObservations %>% 
  left_join(., gmis, by = c("Population", "Latitude", "Longitude"))

# Use haversine formulation to add distance to urban core
source("scripts/haversine.R")

# Create variables to hold lat/long for Toronto urban core. 
# Core assumed to be Young/Dundas square similar to Thompson et al., (2016)
lat_city <- 43.6561
long_city <- -79.3803

# Add distances to dataset
pollinatorObservations <- pollinatorObservations %>%
  mutate(Distance_to_core = haversine(Longitude, Latitude, long_city, lat_city))

# Write data
write_csv(pollinatorObservations, "data-clean/pollinatorObservations.csv")


