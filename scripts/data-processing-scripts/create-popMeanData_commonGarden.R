library(tidyverse)

# Load in experimental data on individual plants
familyMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")

# Create population-mean dataset
popMeans <- familyMeans %>%
  group_by(Population, Latitude, Longitude, gmis, Distance_to_core, Distance_to_cg) %>%
  summarize(Num_Plants = n(),
            Time_to_germination = mean(Time_to_germination, na.rm = TRUE),
            Days_to_flower = mean(Days_to_flower, na.rm = TRUE),
            Num_Inf = mean(Num_Inf, na.rm = TRUE),
            Reprod_biomass = mean(Reprod_biomass, na.rm = TRUE),
            Veget_biomass = mean(Veget_biomass, na.rm = TRUE),
            Avg_bnr_wdth = mean(Avg_bnr_wdth, na.rm = TRUE),
            Avg_bnr_lgth = mean(Avg_bnr_lgth, na.rm = TRUE),
            Avg_petiole_lgth = mean(Avg_petiole_lgth, na.rm = TRUE),
            Avg_peducle_lgth = mean(Avg_peducle_lgth, na.rm = TRUE),
            Avg_num_flwrs = mean(Avg_num_flwrs, na.rm = TRUE),
            Avg_leaf_wdth = mean(Avg_leaf_wdth, na.rm = TRUE),
            Avg_leaf_lgth = mean(Avg_leaf_lgth, na.rm = TRUE),
            Avg_stolon_thick = mean(Avg_stolon_thick, na.rm = TRUE),
            Avg_seeds_per_flower = mean(Avg_seeds_per_flower, na.rm = TRUE),
            Avg_seeds_per_inf = mean(Avg_seeds_per_inf, na.rm = TRUE),
            Avg_num_flwrs = mean(Avg_num_flwrs, na.rm = TRUE),
            freqHCN = mean(freqHCN, na.rm = TRUE),
            sex_asex = Reprod_biomass / Veget_biomass,
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
            Avg_stolon_thick_C = mean(Avg_stolon_thick_C, na.rm = TRUE),
            freqHCN_C = mean(freqHCN_C, na.rm = TRUE)) 

write_csv(popMeans, "data-clean/experimentalData_popMeans.csv")
