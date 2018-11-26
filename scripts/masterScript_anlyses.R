
# Load require packages
library(lme4)
library(lmerTest)

#### TRAIT CHANGE WITH URBANIZATION: ALL DATA ####

# Load in family means dataset
commonGardenData <- read_csv("data-clean/experimentalData_individualPlants.csv")


## MODELS ##

# GERMINATION TIME
germTime <- lmer(Time_to_germination ~ Distance_to_core*HCN_Results + 
                   (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Residuals right skewed and resids vs. fitted shows fanning
germTime_resids <- residuals(germTime)
hist(germTime_resids)
plot(germTime)

# Summarize. Distance to core significant
summary(germTime)

# DAYS TO FIRST FLOWER
firstFlower <- lmer(Days_to_flower ~ Distance_to_core*HCN_Results + 
                      (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good
firstFlower_resids <- residuals(firstFlower)
hist(firstFlower_resids)
plot(firstFlower)

# Summarize. No effects.
summary(firstFlower)

# NUMBER OF INFLORESCENCES
numInf <- lmer(Num_Inf ~ Distance_to_core*HCN_Results + 
                      (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good. Some fanning.
numInf_resids <- residuals(numInf)
hist(numInf_resids)
plot(numInf)

# Summarize. No effects.
summary(numInf)

# REPRODUCTIVE BIOMASS
repBio <- lmer(Reprod_biomass ~ Distance_to_core*HCN_Results + 
                 (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good. Some fanning.
repBio_resids <- residuals(repBio)
hist(repBio_resids)
plot(repBio)

# Summarize
summary(repBio)

# VEGETATIVE BIOMASS. Singular fit issue. No variance assigned to effect of Family.
vegBio <- lmer(Veget_biomass ~ Distance_to_core*HCN_Results + 
                 (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
vegBio_resids <- residuals(vegBio)
hist(vegBio_resids)
plot(vegBio)

# Summarize. Disntance significant.
summary(vegBio)

# BANNER WIDTH
bnrWdth <- lmer(Avg_bnr_wdth ~ Distance_to_core*HCN_Results + 
                 (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
bnrWdth_resids <- residuals(bnrWdth)
hist(bnrWdth_resids)
plot(bnrWdth)

# Summarize. Disntance, HCN, and interaction significant
summary(bnrWdth)

# BANNER LENGTH
bnrLgth <- lmer(Avg_bnr_lgth ~ Distance_to_core*HCN_Results + 
                  (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
bnrLgth_resids <- residuals(bnrLgth)
hist(bnrLgth_resids)
plot(bnrLgth)

# Summarize. Disntance, HCN, and interaction significant
summary(bnrLgth)

# PETIOLE LENGTH
petLgth <- lmer(Avg_petiole_lgth ~ Distance_to_core*HCN_Results + 
                  (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
petLgth_resids <- residuals(petLgth)
hist(petLgth_resids)
plot(petLgth)

# Summarize. No effects
summary(petLgth)

# PEDUNCLE LENGTH. Singular fit issue. No variance assigned to effect of Family.
pedLgth <- lmer(Avg_peducle_lgth ~ Distance_to_core*HCN_Results + 
                  (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
pedLgth_resids <- residuals(pedLgth)
hist(pedLgth_resids)
plot(pedLgth)

# Summarize. No effects
summary(pedLgth)

# NUMBER OF FLOWERS
numFlwr <- lmer(Avg_num_flwrs ~ Distance_to_core*HCN_Results + 
                  (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good. Slight left skew.
numFlwr_resids <- residuals(numFlwr)
hist(numFlwr_resids)
plot(numFlwr)

# Summarize. No effects
summary(numFlwr)

# LEAF WIDTH
leafWdth <- lmer(Avg_leaf_wdth ~ Distance_to_core*HCN_Results + 
                  (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
leafWdth_resids <- residuals(leafWdth)
hist(leafWdth_resids)
plot(leafWdth)

# Summarize. No effects
summary(leafWdth)

# LEAF LENGTH
leafLgth <- lmer(Avg_leaf_lgth ~ Distance_to_core*HCN_Results + 
                   (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
leafLgth_resids <- residuals(leafLgth)
hist(leafLgth_resids)
plot(leafLgth)

# Summarize. Distance x HCN interaction marginal.
summary(leafLgth)

# STOLON THICKNESS
stolThick <- lmer(Avg_stolon_thick ~ Distance_to_core*HCN_Results + 
                   (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
stolThick_resids <- residuals(stolThick)
hist(stolThick_resids)
plot(stolThick)

# Summarize. No effects.
summary(stolThick)

# STOLON THICKNESS
seedPerFlwr <- lmer(Avg_seeds_per_flower ~ Distance_to_core*HCN_Results + 
                    (1|Population/Family_ID), data = commonGardenData)

# Plot histogram of residuals and resudials vs. fitted.
# Looks good.
seedPerFlwr_resids <- residuals(seedPerFlwr)
hist(seedPerFlwr_resids)
plot(seedPerFlwr)

# Summarize. Distance marginal.
summary(seedPerFlwr)

#### FIGURES ####

## TRAIT CHANGES WITH URBANIZATION: ALL DATA ##

