library(tidyverse)

# Load in experimental data on individual plants
indPlantData <- read_csv("data-clean/experimentalData_individualPlants.csv")

# Create population-mean dataset
popMeans <- indPlantData %>%
  group_by(Population, Latitude, Longitude, Distance_to_core, Distance_to_cg) %>%
  summarize(Num_Plants = n(),
            Germination = mean(Time_to_germination, na.rm = TRUE),
            Days_to_flower = mean(Days_to_flower, na.rm = TRUE),
            Num_Inf = mean(Num_Inf, na.rm = TRUE),
            Reprod_biomass = mean(Reprod_biomass, na.rm = TRUE),
            Veget_biomass = mean(Veget_biomass, na.rm = TRUE),
            Bnr_wdth = mean(Avg_bnr_wdth, na.rm = TRUE),
            Bnr_lgth = mean(Avg_bnr_lgth, na.rm = TRUE),
            Petiole_lgth = mean(Avg_petiole_lgth, na.rm = TRUE),
            Peducle_lgth = mean(Avg_peducle_lgth, na.rm = TRUE),
            Num_flwrs = mean(Avg_num_flwrs, na.rm = TRUE),
            Leaf_wdth = mean(Avg_leaf_wdth, na.rm = TRUE),
            Leaf_lgth = mean(Avg_leaf_lgth, na.rm = TRUE),
            Stolon_thick = mean(Avg_stolon_thick, na.rm = TRUE),
            Seeds_per_flower = mean(Avg_seeds_per_flower, na.rm = TRUE),
            Num_flwrs = mean(Avg_num_flwrs, na.rm = TRUE),
            Num_Cyano = sum(HCN_Results, na.rm = TRUE),
            FreqHCN = Num_Cyano / Num_Plants)

# Write data
write_csv(popMeans, "data-clean/experimentalData_popMeans.csv")
