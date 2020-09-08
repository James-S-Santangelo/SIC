# Script to generate family-mean trait data from data collected on individual plants

library(tidyverse)

# Load in clean dataset
experimental_data <- read_csv("data-clean/experimentalData_individualPlants.csv")

# Create family means dataset
familyMeans <- experimental_data %>%
  group_by(Family, Latitude, Longitude, gmis1000, gmis30, Distance_to_core, Distance_to_cg) %>%
  summarize(num_plants = n(),
            Time_to_germination = mean(Time_to_germination, na.rm = TRUE),
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
            Avg_seeds_per_inf = mean(Avg_seeds_per_inf, na.rm = TRUE),
            sex_asex = Reprod_biomass / Veget_biomass,
            n_HCN = sum(HCN_Results),
            total_plants = n(),
            freqHCN = n_HCN / total_plants,
            Time_to_germination_C = mean(Time_to_germination_C, na.rm = TRUE),
            Days_to_flower_C = mean(Days_to_flower_C, na.rm = TRUE),
            Num_Inf_C = mean(Num_Inf_C, na.rm = TRUE),
            Reprod_biomass_C = mean(Reprod_biomass_C, na.rm = TRUE),
            Veget_biomass_C = mean(Veget_biomass_C, na.rm = TRUE),
            Avg_bnr_wdth_C = mean(Avg_bnr_wdth_C, na.rm = TRUE),
            Avg_bnr_lgth_C = mean(Avg_bnr_lgth_C, na.rm = TRUE),
            Avg_petiole_lgth_C = mean(Avg_petiole_lgth_C, na.rm = TRUE),
            Avg_peducle_lgth_C = mean(Avg_peducle_lgth_C, na.rm = TRUE),
            Avg_num_flwrs_C = mean(Avg_num_flwrs_C, na.rm = TRUE),
            Avg_leaf_wdth_C = mean(Avg_leaf_wdth_C, na.rm = TRUE),
            Avg_leaf_lgth_C = mean(Avg_leaf_lgth_C, na.rm = TRUE),
            Avg_stolon_thick_C = mean(Avg_stolon_thick_C, na.rm = TRUE)) %>% 
  separate(Family, sep = "_", into = c("Population", "Family_ID"), remove = FALSE) %>% 
  dplyr::select(-Family_ID) %>% 
  ungroup() %>% 
  mutate(freqHCN_C = freqHCN / mean(freqHCN, na.rm = TRUE))

# Write family means to disk
write_csv(familyMeans, "data-clean/experimentalData_familyMeans.csv")
