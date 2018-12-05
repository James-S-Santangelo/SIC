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
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |

### Family mean dataset (i.e. mean traits across all seeds within a plant family within populations) for common garden. [experimentalData_familyMeans.csv](/data-clean/experimentalData_familyMeans.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population fro, which plant was collected_ Ranges from 1 to 27 | Integer |
| Family_ID | Plant family within population from which plant originates. Ranges from 1 to 8 | Integer |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |
| Time_to_germination | Number of days to germination of seed | Integer |
| Days_to_flower | Average number of days to opening of first flower | Integer |
| Num_Inf | Average number of inflorescences produced by plant | Integer |
| Reprod_biomass | Average total biomass (g) of reproductive tissue produced by plant | Float |
| Veget_biomass| Average total biomass (g) of vegetative biomass produced by plants | Float |
| Avg_bnr_wdth| Average width (mm) of banner petals | Float |
| Avg_bnr_lgth | Average length (mm) of banner petals | Float |
| Avg_peduncle_lgth | Average length (mm) of peduncle (tissue from stolon to inflorescence | Float |
| Avg_num_flwrs | Average number of flowers produced by plant | Float |
| Avg_leaf_wdth | Average width of leaflets (mm) produced by plant | Float |
| Avg_leaf_lgth | Average length of leaflets (mm) produced by plant | Float |
| Avg_petiole_lgth| Average length (mm) of petiole (tissue from stolon to leaf) | Float |
| Avg_stolon_thick | Average thickness (i.e. diameter, mm) of stolons produced by plant | Float |
| Avg_seeds_per_flower | Average number of seeds produced per flower | Float |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |

### Population mean dataset (i.e. mean traits across all seeds within populations, pooled across families) for common garden. [experimentalData_popMeans.csv](/data-clean/experimentalData_popMeans.csv)

| Column | Description | Type |
|--------|-------------|------|
| Population  | Population from which plant was collected_ Ranges from 1 to 27 | Integer |
| Latitude | Decimal degree latitude of population from which plant was collected | Float |
| Longitude| Decimal degree longitude of population from which plant was collected | Float |
| Distance_to_core  | Distance of plant population to city center (km) | Float |
| Distance_to_cg | Distance of plant population to common garden at the University of Toronto Mississauga | Float |
| Num_Plants | Number of plants in the population | Integer |
| Germination | Average number of days to germination of seed | Integer |
| Days_to_flower | Average number of days to opening of first flower | Integer |
| Num_Inf | Average number of inflorescences produced by plant | Integer |
| Reprod_biomass | Average total biomass (g) of reproductive tissue produced by plant | Float |
| Veget_biomass| Average total biomass (g) of vegetative biomass produced by plants | Float |
| Bnr_wdth| Average width (mm) of banner petals | Float |
| Bnr_lgth | Average length (mm) of banner petals | Float |
| Petiole_lgth| Average length (mm) of petiole (tissue from stolon to leaf) | Float |
| Peduncle_lgth | Average length (mm) of peduncle (tissue from stolon to inflorescence | Float |
| Num_flwrs | Average number of flowers produced by plant | Float |
| Leaf_wdth | Average width of leaflets (mm) produced by plant | Float |
| Leaf_lgth | Average length of leaflets (mm) produced by plant | Float |
| Stolon_thick | Average thickness (i.e. diameter, mm) of stolons produced by plant | Float |
| Seeds_per_flower | Average number of seeds produced per flower | Float |

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

### Pollination efficiency in the field dataset for individual inflorescences. [seedFlwrRatio_fieldPlants_cleaned.csv](/data-clean/seedFlwrRatio_fieldPlants_cleaned.csv)

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

