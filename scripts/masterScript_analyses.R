# install packages
# install.packages("lme4")
# install.packages("lmerTest")

# Load require packages
library(lme4)
library(lmerTest)

#### MODELS FOR TRAIT CHANGE WITH URBANIZATION: ALL DATA ####

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

# Summarize. Distance significant.
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

#Theme used for plots throughout script
ng1=theme(aspect.ratio=0.7,panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1),
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black",size=15),
          axis.title=element_text(color="black",size=1),
          axis.title.y=element_text(vjust=2,size=17),
          axis.title.x=element_text(vjust=0.1,size=17),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
          strip.background = element_rect(colour="black"),
          legend.position = "top", legend.direction="vertical",
          legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))

## TRAIT CHANGES WITH URBANIZATION: ALL DATA ##

# Plotting effects that were significant from linear models above

# GERMINATION TIME
# Relationship driven by high germination time outliers in urban pops.
germTime_plot <- ggplot(commonGardenData, aes(x = Distance_to_core, y = Time_to_germination)) +
  geom_point(size = 1.5, colour = "black", position = position_dodge(width = 0.1)) +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
germTime_plot

ggsave("analysis/figures/traitChanges_allData/germinationTime_allData.pdf", 
       plot = germTime_plot, width = 5, height = 5, unit = "in", dpi = 600)

# Germination time pattern easier to see if visualizing pop-means
germTime_popMeanPlot <- commonGardenData %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Time_to_germination, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
germTime_popMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/germinationTime_popMeans.pdf", 
       plot = germTime_popMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

# BIOMASS
vegBio_plot <- ggplot(commonGardenData, aes(x = Distance_to_core, y = Veget_biomass)) +
  geom_point(size = 1.5, colour = "black", position = position_dodge(width = 0.1)) +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
vegBio_plot

ggsave("analysis/figures/traitChanges_allData/vegetativeBiomass_allData.pdf", 
       plot = vegBio_plot, width = 5, height = 5, unit = "in", dpi = 600)

# Biomass pattern easier to see if visualizing pop-means
vegBio_popMeanPlot <- commonGardenData %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Veget_biomass, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
vegBio_popMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/vegetativeBiomass_popMeans.pdf", 
       plot = vegBio_popMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

# BANNER WIDTH
commonGardenData$HCN_Results <- as.factor(commonGardenData$HCN_Results)
bnrWdth_plot <- ggplot(commonGardenData, aes(x = Distance_to_core, y = Avg_bnr_wdth, group = HCN_Results)) +
  geom_point(size = 1.5, aes(colour = HCN_Results, shape = HCN_Results), position = position_dodge(width = 0.1)) +
  geom_smooth(method = "lm", size = 2.0, aes(linetype = HCN_Results, colour = HCN_Results), se = FALSE) +
  ng1
bnrWdth_plot

ggsave("analysis/figures/traitChanges_allData/bannerWidth-HCN_allData.pdf", 
       plot = bnrWdth_plot, width = 6, height = 7, unit = "in", dpi = 600)

# Population mean banner width for HCN+ plants
bnrWdth_CyanPopMeanPlot <- commonGardenData %>%
  filter(HCN_Results == "1") %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Avg_bnr_wdth, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
bnrWdth_CyanPopMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/bannerWidth-Cyan_popMeans.pdf", 
       plot = bnrWdth_CyanPopMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

# Population mean banner width for HCN- plants
bnrWdth_AcyanPopMeanPlot <- commonGardenData %>%
  filter(HCN_Results == "0") %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Avg_bnr_wdth, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
bnrWdth_AcyanPopMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/bannerWidth-Acyan_popMeans.pdf", 
       plot = bnrWdth_AcyanPopMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

# BANNER LENGTH
bnrLgth_plot <- ggplot(commonGardenData, aes(x = Distance_to_core, y = Avg_bnr_lgth, group = HCN_Results)) +
  geom_point(size = 1.5, aes(colour = HCN_Results, shape = HCN_Results), position = position_dodge(width = 0.1)) +
  geom_smooth(method = "lm", size = 2.0, aes(linetype = HCN_Results, colour = HCN_Results), se = FALSE) +
  ng1
bnrLgth_plot

ggsave("analysis/figures/traitChanges_allData/bannerLength-HCN_allData.pdf", 
       plot = bnrLgth_plot, width = 6, height = 7, unit = "in", dpi = 600)

# Population mean banner width for HCN+ plants
bnrLgth_CyanPopMeanPlot <- commonGardenData %>%
  filter(HCN_Results == "1") %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Avg_bnr_lgth, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
bnrLgth_CyanPopMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/bannerLength-Cyan_popMeans.pdf", 
       plot = bnrLgth_CyanPopMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

# Population mean banner width for HCN- plants
bnrLgth_AcyanPopMeanPlot <- commonGardenData %>%
  filter(HCN_Results == "0") %>%
  group_by(Population, Distance_to_core) %>%
  summarize(mean = mean(Avg_bnr_lgth, na.rm = TRUE)) %>%
  ggplot(., aes(x = Distance_to_core, y = mean)) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
bnrLgth_AcyanPopMeanPlot

ggsave("analysis/figures/traitChanges_popMeans/bannerLength-Acyan_popMeans.pdf", 
       plot = bnrLgth_AcyanPopMeanPlot, width = 5, height = 5, unit = "in", dpi = 600)

#### TRAIT CORRELATIONS ####

# Will assess trait correlations using family means

familyMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")
traits <- familyMeans %>%
  select(Time_to_germination, Days_to_flower, Num_Inf, Reprod_biomass,
         Veget_biomass, Avg_bnr_wdth, Avg_bnr_lgth, Avg_leaf_lgth,
         Avg_peducle_lgth, Avg_num_flwrs, Avg_leaf_wdth, Avg_petiole_lgth,
         Avg_stolon_thick, Avg_seeds_per_flower)

## FUNCTIONS ##

##Function for changing upper panel in 'pairs' correlation matrix to show correlation
#coefficients and p-values
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use = "complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

#Function to add least squares regression line to lower panel of 'pairs' scatterplot matrix
lsline = function(x,y) {
  points(x,y,pch=".")
  abline(lsfit(x,y),col="blue")
}

# Trait correlation table
pairs(traits, upper.panel = panel.cor,lower.panel=lsline)

# With a few exceptions, traits are not strongly correlated
# Bnr length and bnr width are correlated (r = 0.76) and both traits show the same
# patten with urbanization. Should these be condensed into a single traits (i.e. flower size).
# Same is true for leaft length and leaf width.

#### SEEDS PER FLOWER: FIELD DATA ####

seedFlwrFieldData <- read_csv("data-clean/flwrSeedRatio_fieldPlants_cleaned.csv") %>%
  select(-X1, -Comments)

seedPerFlwrField <- lmer(Seeds_per_flower ~ Distance_to_core + (1|Population), 
                         data = seedFlwrFieldData)

# Plot histogram of residuals and resudials vs. fitted.
# Hist looks good. Some fanning in fitted plot.
seedPerFlwrField_resids <- residuals(seedPerFlwrField)
hist(seedPerFlwrField_resids)
plot(seedPerFlwrField)

# Summarize. Distance marginal
summary(seedPerFlwrField)


#### POLLINATOR OBSERVATIONS ####
