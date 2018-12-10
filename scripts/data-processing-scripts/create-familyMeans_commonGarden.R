library(tidyverse)

# Load in clean dataset
experimental_data <- read_csv("data-clean/experimentalData_individualPlants.csv")

# Create family means dataset
familyMeans <- experimental_data %>%
  group_by(Population, Family_ID, Distance_to_core, Distance_to_cg) %>%
  summarize(Time_to_germination = mean(Time_to_germination, na.rm = TRUE),
            Days_to_flower = mean(Days_to_flower, na.rm = TRUE),
            Num_Inf = mean(Num_Inf, na.rm = TRUE),
            Reprod_biomass = mean(Reprod_biomass, na.rm = TRUE),
            Veget_biomass = mean(Veget_biomass, na.rm = TRUE),
            Avg_bnr_wdth = mean(Avg_bnr_wdth, na.rm = T),
            Avg_bnr_lgth = mean(Avg_bnr_lgth, na.rm = T),
            Avg_peducle_lgth = mean(Avg_peducle_lgth, na.rm = T),
            Avg_num_flwrs = mean(Avg_num_flwrs, na.rm = T),
            Avg_leaf_wdth = mean(Avg_leaf_wdth, na.rm = T),
            Avg_leaf_lgth = mean(Avg_leaf_lgth, na.rm = T),
            Avg_petiole_lgth = mean(Avg_petiole_lgth, na.rm = T),
            Avg_stolon_thick = mean(Avg_stolon_thick, na.rm = T),
            Avg_seeds_per_flower = mean(Avg_seeds_per_flower, na.rm = TRUE),
            n_HCN = sum(HCN_Results),
            total_plants = n(),
            freqHCN = n_HCN / total_plants)

# Write family means to disk
write_csv(familyMeans, "data-clean/experimentalData_familyMeans.csv")
