# Script to generat both family-mean and population-mean pollination efficiency datasets for field-collected plants

library(tidyverse)

# Load in flwr/seed data from field-collected inflorescences
flwrSeedFieldData <- read_csv("data-raw/flwrSeedRatio_fieldPlants.csv", na = c("", "NA"))

# Create clean data with all individuals
seedFlwrRatio_cleaned <- flwrSeedFieldData %>%
  
  # Add in data first first individual in population 1 that was lost due to interruption
  # at field site.
  mutate(Num.Flwrs = replace(Num.Flwrs, Population == 1 & Individual == 1, 20),
         Num.Seeds = replace(Num.Seeds, Population == 1 & Individual == 1, 21),
         Comments = replace(Comments, Population == 1 & Individual == 1, NA)) %>%
  
  # Keep only rows with no comments. Comments represent missing/poor observations (N = 4)
  filter(is.na(Comments)) %>%
  
  # Define number of seeds per flower
  mutate(Seeds_per_flower = Num.Seeds / Num.Flwrs,
         Population = as.character(Population))

# Load in data with lat/longs
latLongs <- read_csv("data-raw/populationLatLongs.csv") %>% 
  mutate(Population = as.character(Population))

# Add population lat/long data to seed/flower ratio dataframe
seedFlwrRatio_cleaned <- seedFlwrRatio_cleaned %>%
  left_join(., latLongs, by = "Population")

# Add impervious survafe values to data 
seedFlwrRatio_cleaned <- seedFlwrRatio_cleaned %>% 
  left_join(., gmis, by = c("Population", "Latitude", "Longitude"))

# Use haversine formulation to add distance to urban core
source("scripts/haversine.R")

# Create variables to hold lat/long for Toronto urban core. 
# Core assumed to be Young/Dundas square similar to Thompson et al., (2016)
lat_city <- 43.6561
long_city <- -79.3803

# Add distances to dataset
seedFlwrRatio_cleaned <- seedFlwrRatio_cleaned %>%
  mutate(Distance_to_core = haversine(Longitude, Latitude, long_city, lat_city))

# Write clean data to disk
write.csv(seedFlwrRatio_cleaned, "data-clean/flwrSeedRatio_fieldPlants.csv")
  
# Creat population mean dataset
seedFlwrRatio_popMeans <- seedFlwrRatio_cleaned %>%
  
  # Group by population
  group_by(Population, Latitude, Longitude, Distance_to_core) %>%
  
  # Calculte mean number of seeds per flower
  summarize(Num_Seeds = mean(Num.Seeds, na.rm = TRUE),
            Num_Flowers = mean(Num.Flwrs, na.rm = TRUE),
            Seeds_per_flower = mean(Seeds_per_flower, na.rm = TRUE))
  
# Write data
write_csv(seedFlwrRatio_popMeans, "data-clean/seedFlowerRatio_popMeans_fieldData.csv")
