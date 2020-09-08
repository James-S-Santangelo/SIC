# Metadata for datasets
##### Manuscript: Multivariate phenotypic trait divergence along an urbanization gradient (tentative title)
##### Authors: James S. Santangelo, Ken A. Thompson, Carole Advenard, and L. Ruth Rivkin
##### Journal:

This file contains the metadata for individual columns/variables for each of the [clean datasets](/data-clean) used throughout the analyses_

### Data measured on individual plants in the common garden. [experimentalData_individualPlants.csv](/data-clean/experimentalData_individualPlants.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population fro, which plant was collected_ Ranges from 1 to 27 | Integer |
| Family_ID | Plant family within population from which plant originates. Ranges from 1 to 8 | Integer |
| Seed | Seed of plant from plant family_ Ranges from 1 to 5 | Integer |
| label | Label of plant at field site_ Concatenation of "Population", "Family_ID", and "Seed" | Character/String |
| Time_to_germination | Number of days to germination of seed | Integer |
| Row | Row position of plant in common garden. Range: 1-21 | Integer |
| Column | Column position of plant in common garden. Range: 1-32 | Integer |
| Days_to_flower | Number of days to opening of first flower | Integer |
| Num_Inf | Number of inflorescences produced by plant | Integer |
| HCN_Results | Presence (1) or absence (0) of HCN. | Integer |
| Reprod_biomass | Total biomass (g) of reproductive tissue produced by plant | Float |
| Veget_biomass| Total biomass (g) of vegetative biomass produced by plants | Float |
| Avg_bnr_wdth| Average width (mm) of banner petals | Float |
| Avg_bnr_lgth | Average length (mm) of banner petals | Float |
| Avg_petiole_lgth| Average length (mm) of petiole (tissue from stolon to leaf) | Float |
| Avg_peduncle_lgth | Average length (mm) of peduncle (tissue from stolon to inflorescence | Float |
| Avg_num_flwrs | Average number of flowers produced by plant | Float |
| Avg_leaf_wdth | Average width of leaflets (mm) produced by plant | Float |
| Avg_leaf_lgth | Average length of leaflets (mm) produced by plant | Float |
| Avg_stolon_thick | Average thickness (i.e. diameter, mm) of stolons produced by plant | Float |
| Avg_seeds_per_flower | Average number of seeds produced per flower | Float |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| gmis1000 | Percent impervious surface using a 1000 meter buffer | Float |
| gmis30 | Percent impervious surface using a 30 meter buffer | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |
| Time_to_germination_C : Avg_stolon_thick_C | Same as traits above but standardized by dividing by experiment-wise means | Float |
| sex_asex | Ratio of `Veget_biomass` and `Reprod_biomass` | Float |
| Family | Unique family identifier. Concatenation of `Family_ID` and `Populations` | Character/String |

### Family mean dataset (i.e. mean traits across all seeds within a plant family within populations) for common garden. [experimentalData_familyMeans.csv](/data-clean/experimentalData_familyMeans.csv)

| Column | Description | Type |
|--------|-------------|------|
| Family | Unique family identifier | Character/String |
| Population  | Population fro, which plant was collected_ Ranges from 1 to 27 | Integer |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| gmis1000 | Percent impervious surface using a 1000 meter buffer | Float |
| gmis30 | Percent impervious surface using a 30 meter buffer | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |
| Num_plants | Number of plants in `Family` | Integer |
| Time_to_germination | Family-mean time to germination (days) | Float |
| Days_to_flower | Family-mean days to opening of first flower | Float |
| Num_Inf | Family-mean number of inflorescences produced by plant | Float |
| Reprod_biomass | Family-mean total biomass (g) of reproductive tissue | Float |
| Veget_biomass| Family-mean total biomass (g) of vegetative biomass | Float |
| Avg_bnr_wdth| Family-mean width (mm) of banner petals | Float |
| Avg_bnr_lgth | Family-mean length (mm) of banner petals | Float |
| Avg_peduncle_lgth | Family-mean length (mm) of peduncle (tissue from stolon to inflorescence) | Float |
| Avg_num_flwrs | Family-mean number of flowers | Float |
| Avg_leaf_wdth | Family-mean width of leaflets (mm) | Float |
| Avg_leaf_lgth | Family-mean length of leaflets (mm) | Float |
| Avg_petiole_lgth | Family-mean length (mm) of petiole (tissue from stolon to leaf) | Float |
| Avg_stolon_thick | Family-mean thickness (i.e. diameter, mm) of stolons | Float |
| Avg_seeds_per_flower | Family-mean number of seeds produced per flower | Float |
| sex_asex | Family-mean ratio of vegetative to sexual reproductive biomass | Float |
| n_HCN | Number of cyanogenic plants in `Family` | Float |
| total_plants | Number of plants in `Family`. Same as `Num_plants` | Float |
| freqHCN | Family-mean frequency of HCN | Float |
| Time_to_germination_C : freq_HCN_C | Same as traits above but standardized by dividing by experiment-wise means | Float |

### Population mean dataset (i.e. mean of family-means for all traits within populations) for common garden. [experimentalData_popMeans.csv](/data-clean/experimentalData_popMeans.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| gmis1000 | Percent impervious surface using a 1000 meter buffer | Float |
| gmis30 | Percent impervious surface using a 30 meter buffer | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |
| Num_families | Number of plant families in the population | Integer |
| Time_to_germination | Population-mean number of days to germination of seed | Float |
| Days_to_flower | Population-mean number of days to opening of first flower | Float |
| Num_Inf | Population-mean number of inflorescences produced by plant | Integer |
| Reprod_biomass | Population-mean total biomass (g) of reproductive tissue produced by plant | Float |
| Veget_biomass| Population-mean total biomass (g) of vegetative biomass produced by plants | Float |
| Avg_bnr_wdth| Population-mean width (mm) of banner petals | Float |
| Avg_bnr_lgth | Population-mean length (mm) of banner petals | Float |
| Avg_petiole_lgth| Population-mean length (mm) of petiole (tissue from stolon to leaf) | Float |
| Avg_peduncle_lgth | Population-mean length (mm) of peduncle (tissue from stolon to inflorescence | Float |
| Avg_num_flwrs | Population-mean number of flowers produced by plant | Float |
| Avg_leaf_wdth | Population-mean width of leaflets (mm) produced by plant | Float |
| Avg_leaf_lgth | Population-mean length of leaflets (mm) produced by plant | Float |
| Avg_stolon_thick | Population-mean thickness (i.e. diameter, mm) of stolons produced by plant | Float |
| Avg_seeds_per_flower | Population-mean number of seeds produced per flower | Float |
| freqHCN | Population-mean frequency of HCN | Float |
| sex_asex | Population-mean ratio of vegetative to sexual reproductive biomass | Float |
| Time_to_germination_C : freq_HCN_C | Same as traits above but standardized by dividing by experiment-wise means | Float |

### Pollinator observation dataset. [pollinatorObservations.csv](/data-clean/pollinatorObservations.csv)

| Column | Description | Type |
|--------|-------------|------|
| Date | Date on which observations were performed. Day/Month/Year. | Character/String |
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Num_Inf | Number of inflorescences in square meter plot in which observations were being performed | Integer |
| Pollinator | Number for the pollinator entering the plot (e.g. 1 for first pollinator, 2 for second, etc.) | Integer |
| Num_Visit | Number of inflorescences in plot visited by pollinator | Integer |
| Morph | Pollinator morphospecies/taxon | Character/String |
| Observer | Name of person permorming observations on plot | Character/String |
| Comments | Comments/notes associated with pollinator observations | Character/String |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| gmis1000 | Percent impervious surface using a 1000 meter buffer | Float |
| gmis30 | Percent impervious surface using a 30 meter buffer | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |

### Population mean pollination efficiency in the field dataset. [seedFlowerRatio_popMeans_fieldData.csv](/data-clean/seedFlowerRatio_popMeans_fieldData.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Num_Seeds  | Average number of seeds produced by plants in the population | Float |
| Num_Flowers  | Average number of flowers produced by plants in the population | Float |
| Seeds_per_flower  | Average number of seeds per flower produced by plants in the populations | Float |

### Pollination efficiency in the field dataset for individual inflorescences. [flwrSeedRatio_fieldPlants_cleaned.csv](/data-clean/flwrSeedRatio_fieldPlants.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Individual  | Individual collected from population. Ranges from 1 to 20 | Integer |
| Num.Flwrs  | Number of flowers produced by individual | Float |
| Num.Seeds  | Number of seeds produced by individual | Float |
| Comments | Comments/notes individual inflorescence | Character/String |
| Seeds.per.flower  | Number of seeds per flower produced by individual | Float |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |

### Plant germination dataset [germination.csv](/data-clean/germination.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Maternal_plant  | ID of maternal plant from which seeds were collected | Integer |
| Family_ID  | Same as `Maternal_plant` but generated new IDs | Float |
| Seed | Seed ID from Maternal plant | Integer |
| Date_planted | Date seed was planted in soil | Character/String |
| Germination_data  | Date seed germinated | Character/String |
| Time to germination | Number of days from `Date_planted` to `Germination_date` | Integer |
| Is.germinated | Whether plant did (1) of did not (0) germinate | Float |

### Impervious surface cover dataset [gmis.csv](/data-clean/gmis.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| gmis1000 | Percent impervious surface using a 1000 meter buffer | Float |
| gmis30 | Percent impervious surface using a 30 meter buffer | Float |

### Dataset with a couple "environmental" covariates [enviroData.csv](/data-clean/enviroData.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Imperv | Percent impervious surface using a 1000 meter buffer. Same as `gmis1000` in other datasets | Float |
| popDens | Human population density extracted using 1000 meter buffer | Float |