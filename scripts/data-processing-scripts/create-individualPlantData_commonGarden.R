library(tidyverse)

# Load in experimental data from Commond Garden
experimentalData <- read.csv("data-raw/experimentalData_commonGarden.csv", na = c("", "NA", "#NUM!"))

# Load in data with number of flowers and seeds for infividual plants
flwrSeedData <- read_csv("data-raw/flwrSeedRatio_commonGardenPlants.csv", na = c("", "NA", "#NUM!")) %>%
  select(label, Row, Column, Num_flwr_R1, Num_Seeds_R1, Num_flwr_R2, Num_Seeds_R2) %>%
  filter_all(any_vars(!is.na(.)))

# Load in data with masses for bags and envelopes that held reproductive and vegetative biomass
bagEnv_Masses <- read_csv("data-raw/bagMasses.csv")

# Summarrize mass of empty bags and envelopes
bagEnv_Masses <- bagEnv_Masses %>%
  group_by(Type) %>%
  summarize(meanMass = mean(Mass))

# Extract mean mass of bags and evelopes
bagMass <- bagEnv_Masses %>% filter(Type == "Bag") %>% select(meanMass) %>% pull()
envMass <- bagEnv_Masses %>% filter(Type == "Env") %>% select(meanMass) %>% pull()

# Rename, calculate trait means, select columns.
experimentalData_modified <- experimentalData %>%
  
  # Rename columns
  rename("Seed" = "seeds",
         "Days_to_flower" = "Days_to_first_flower",
         "Num_flwrs1" = "Number.of.flowers_inflorescences.1",
         "Num_flwrs2" = "Number.of.flowers_inflorescences.2",
         "Num_Inf" = "Number.of.inflorescences", 
         "Petiole_length1" = "heigth_leaf1",
         "Petiole_length2" = "heigth_leaf2",
         "Petiole_length3" = "heigth_leaf3",         
         "Peduncle_length1" = "Inflor_height1",
         "Peduncle_length2" = "Inflor_height2",
         "Peduncle_length3" = "Inflor_height3",
         "Reprod_biomass_withEnv" = "weigth_flowers.bag",
         "Veget_biomass_withBag" = "weigth_plant.bag") %>%
  
  # Merge with flwr/seed count dataset
  left_join(., flwrSeedData, by = c("label", "Row", "Column")) %>%
  
  # Define number of seeds per flower
  mutate(Seeds_per_flower1 = Num_Seeds_R1 / Num_flwr_R1,
         Seeds_per_flower2 = Num_Seeds_R2 / Num_flwr_R2) %>%
  
  # Calculate means of traits with multiple measurements for individual plants
  mutate(Reprod_biomass = Reprod_biomass_withEnv - envMass,
         Veget_biomass = Veget_biomass_withBag - bagMass,
         Avg_bnr_wdth = select(., starts_with("Bnr_wdth")) %>% rowMeans(na.rm = TRUE),
         Avg_bnr_lgth = select(., starts_with("Bnr_length")) %>% rowMeans(na.rm = TRUE),
         Avg_petiole_lgth = select(., starts_with("Petiole")) %>% rowMeans(na.rm = TRUE),
         Avg_peducle_lgth = select(., starts_with("Peduncle")) %>% rowMeans(na.rm = TRUE),
         Avg_num_flwrs = select(., starts_with("Num_flwrs")) %>% rowMeans(na.rm = TRUE),
         Avg_leaf_wdth= select(., starts_with("width_leaf")) %>% rowMeans(na.rm = TRUE),
         Avg_leaf_lgth = select(., starts_with("length_leaf")) %>% rowMeans(na.rm = TRUE),
         Avg_stolon_thick = select(., starts_with("Width_stolon")) %>% rowMeans(na.rm = TRUE),
         Avg_seeds_per_flower = select(., starts_with("Seeds_per")) %>% rowMeans(na.rm = TRUE)) %>%
  
  # Select columns that will be used for analyses
  select(Population, Family_ID, Seed, label, Time_to_germination, Row, Column, 
         Days_to_flower, Num_Inf, HCN_Results, Reprod_biomass, Veget_biomass,
         Avg_bnr_wdth, Avg_bnr_lgth, Avg_petiole_lgth, Avg_peducle_lgth, 
         Avg_num_flwrs, Avg_leaf_wdth, Avg_leaf_lgth, Avg_stolon_thick,
         Avg_seeds_per_flower) %>%
  
  # Reprod_biomass and num_flwrs should be 0 if plant never flowered, not NA
  mutate(Reprod_biomass = ifelse(Num_Inf == 0, 0, Reprod_biomass),
         Avg_num_flwrs = ifelse(Num_Inf == 0, 0, Avg_num_flwrs)) %>%
  
  # Replace NaN with NA
  na_if("NaN")
  
# Load in data with lat/longs
latLongs <- read_csv("data-raw/populationLatLongs.csv")

# Add population lat/long data to experimental plants dataframe
experimentalData_modified <- experimentalData_modified %>%
  left_join(., latLongs, by = "Population")

# Use haversine formulation to add distance to urban core and distance to common garden 
source("scripts/haversine.R")

# Create variables to hold lat/long for Toronto urban core. 
# Core assumed to be Young/Dundas square similar to Thompson et al., (2016)
lat_city <- 43.6561
long_city <- -79.3803

# Create variables for lat/long of common garden at UofT Mississauga
lat_cg <- 43.549394
long_cg <- -79.662514

# Add distances to dataset
experimentalData_modified <- experimentalData_modified %>%
  mutate(Distance_to_core = haversine(Longitude, Latitude, long_city, lat_city),
         Distance_to_cg = haversine(Longitude, Latitude, long_cg, lat_cg))

# Write clean data
write_csv(experimentalData_modified, "data-clean/experimentalData_individualPlants.csv")
