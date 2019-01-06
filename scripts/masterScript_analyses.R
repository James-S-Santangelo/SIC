# install packages (these packages are loaded in this script)
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("car")

## Required by other scripts. Will eventually be migrated.
# install.packages("psych")
# install.packages("cowplot")
# install.packages("RcppEigen")
# install.packages("mnormt")
# install.packages("BH")
# install.packages("plogr")

# Commands used to initialize packrat. 
# packrat::init(restart = TRUE, enter = FALSE, infer.dependencies = FALSE)
# packrat::on()
# .libPaths()
# packrat::status()
# packrat::snapshot()

# Load require packages
library(lme4)
library(lmerTest)
library(tidyverse)
library(vegan)
library(car)

#### MULTIVARIATE TRAIT CHANGE WITH URBANIZATION: POPULATION MEANS ####

# Load in family mean dataset
popMeans <- read_csv("data-clean/experimentalData_popMeans.csv")

# Subset popMeans datafor use in RDA
popMeans_forRDA <- popMeans %>%
  select(-Seeds_per_flower, -Num_Cyano) %>%
  select(Population, Distance_to_core, Germination:FreqHCN)

# Perform RDA with multiple traits as response, distance as sole predictors
rdaPop <- rda(popMeans_forRDA %>% 
                     select(Germination:FreqHCN) ~ 
             popMeans_forRDA$Distance_to_core,
         scale = TRUE, na.action = "na.omit")
summary(rdaPop)

# Permutation based test of significance of distance term in RDA
set.seed(42)
anova.cca(rdaPop, by = "term", permutations = 1000)

# Extract canonical coefficients of traits onto RDA1 (i.e. 'species scores')
species_scores <- scores(rdaPop)$species[,"RDA1"]

# Calculate cline_max according to Stock et al. 
# Traits are first standardized by dividing by experiment-wide mean.
# Standardized trait values are multiplied by their cannonical regression
# coefficients on RDA. This is done for each trait and then summed. 
popMeans <- popMeans %>%
  mutate(clinemax = 
           (Germination / mean(Germination)) * species_scores["Germination"] +
           (Days_to_flower / mean(Days_to_flower)) * species_scores["Days_to_flower"] +
           (Num_Inf / mean(Num_Inf)) * species_scores["Num_Inf"] +
           (Reprod_biomass / mean(Reprod_biomass)) * species_scores["Reprod_biomass"] +
           (Veget_biomass / mean(Veget_biomass)) * species_scores["Veget_biomass"] +
           (Bnr_wdth / mean(Bnr_wdth)) * species_scores["Bnr_wdth"] +
           (Bnr_lgth / mean(Bnr_lgth)) * species_scores["Bnr_lgth"] +
           (Petiole_lgth / mean(Petiole_lgth)) * species_scores["Petiole_lgth"] +
           (Peducle_lgth / mean(Peducle_lgth)) * species_scores["Peducle_lgth"] +
           (Num_flwrs / mean(Num_flwrs)) * species_scores["Num_flwrs"] +
           (Leaf_wdth / mean(Leaf_wdth)) * species_scores["Leaf_wdth"] +
           (Leaf_lgth / mean(Leaf_lgth)) * species_scores["Leaf_lgth"] +
           (Stolon_thick / mean(Stolon_thick)) * species_scores["Stolon_thick"] +
           (FreqHCN / mean(FreqHCN)) * species_scores["FreqHCN"])

# Model testing for cline in cline_max
clineMax_mod <- lm(clinemax ~ Distance_to_core, data = popMeans)
summary(clineMax_mod)

#### UNIVARIATE TRAIT CHANGE WITH URBANIZATION: POPULATION MEANS ####

# Sig
germMod <- lm(Germination ~ Distance_to_core, data = popMeans)
summary(germMod)
plot(germMod)
hist(residuals(germMod)) # No transformation needed

# Marg
ffMod <- lm(Days_to_flower^4 ~ Distance_to_core, data = popMeans)
summary(ffMod)
plot(ffMod)
hist(residuals(ffMod)) # Residuals improved following ^4 transformation

# Sig
vegMod <- lm(Veget_biomass^2 ~ Distance_to_core, data = popMeans)
summary(vegMod)
plot(vegMod)
hist(residuals(vegMod)) # Normality improved following ^2 transformation

# Marg
bwMod <- lm(Bnr_wdth ~ Distance_to_core, data = popMeans)
summary(bwMod)
plot(bwMod)
hist(residuals(bwMod)) #  No transformation needed

# Sig
blMod <- lm(Bnr_lgth ~ Distance_to_core, data = popMeans)
summary(blMod)
plot(blMod)
hist(residuals(blMod)) # No transformation needed

# Marg
HCNMod <- lm(sqrt(FreqHCN) ~ Distance_to_core, data = popMeans)
summary(HCNMod)
plot(HCNMod)
hist(residuals(HCNMod)) # Square root transformation improves normality

# NS
stMod <- lm(Stolon_thick ~ Distance_to_core, data = popMeans)
summary(stMod)
plot(stMod)
hist(residuals(stMod)) # No transformation needed

# NS
sexInvestMod <- lm(Reprod_biomass / Veget_biomass ~ Distance_to_core, data = popMeans)
summary(sexInvestMod)
plot(sexInvestMod)
hist(residuals(sexInvestMod)) # No transformation needed

#### POLLINATOR OBSERVATIONS ####

# load data
polldata <- read.csv('data-clean/pollinatorObservations.csv')

# analysis of pollinator observations
# impute zeros for unobserved polls
complete_polldata <- polldata %>% 
  group_by(Population,Morph, Distance_to_core) %>% 
  complete(Population, Morph, fill = list(Num_Visit = 0)) %>%
  filter(Morph != "Syrphid")

# Summarize data by pollinator functional group
# Used for visit/infl analysis
visitshare_polldata <- complete_polldata %>% 
  group_by(Population, Distance_to_core, Morph) %>% 
  # calculate number of individuals observed per site per morph; number of inflorescence, number of visit per infl.
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>%  
  mutate(Visits_per_Inf = Num_Visit / Num_Inf) %>% 
  # replace NA in visits/inf with zero
  mutate(Visits_per_Inf = replace_na(Visits_per_Inf, 0))

# Summarize daya by population only (merge morphs)
# Used for total abundance analysis
popmeans_polldata <- polldata %>% 
  group_by(Population, Distance_to_core) %>% 
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>% 
  mutate(Visits_per_Inf = Num_Visit / Num_Inf)

# Square root transform visits per inflorescence to improve normality of model residuals
visitshare_polldata$sqRootVisInf <- sqrt(visitshare_polldata$Visits_per_Inf)

# Multiple regression analysis on log(visit/inf) with distance, morph, and their interaction. 
pollVisit <- lm(sqRootVisInf ~ Distance_to_core + Morph + Distance_to_core:Morph,
         data = visitshare_polldata)
summary(pollVisit)
Anova(pollVisit, type = "III") #Type 3 for interpreting interaction

# Diagnostic plots
hist(residuals(pollVisit))
par(mfrow = c(2,2))
plot(pollVisit)

#### SEEDS PER FLOWER ####

# Load in seeds per flower from field-collected inflorescences
seedFlwrFieldData <- read_csv("data-clean/flwrSeedRatio_fieldPlants.csv") %>%
  select(-X1, -Comments) %>%
  group_by(Population, Distance_to_core) %>%
  summarize(Seeds_per_flower = mean(Seeds_per_flower)) 

# Retrieve seeds per flower from common garden plants
seedFlwrComGarData <- popMeans %>%
  select(Population, Distance_to_core, Seeds_per_flower)

# Combine both dataframes to run single linear model
seedsPerFlwr <- bind_rows(seedFlwrFieldData, seedFlwrComGarData,
                          .id = "source") %>%
  mutate(source = fct_recode(source, "field" = "1", "common garden" = "2"))

# Model testing for differences in seeds per flower among field and CG plants
seedPerFlwrField <- lm(Seeds_per_flower ~ Distance_to_core*source, 
                         data = seedsPerFlwr)
summary(seedPerFlwrField)
Anova(seedPerFlwrField, type = "III")
plot(seedPerFlwrField) # Loogs good
hist(residuals(seedPerFlwrField)) # Loogs good

# Means
seedsPerFlwr %>%
  group_by(source) %>%
  summarize(mean = mean(Seeds_per_flower))

# Get betas for individual sources
seedMods <- seedsPerFlwr %>%
  group_by(source) %>%
  do(mod = lm(Seeds_per_flower ~ Distance_to_core, data = .)) %>%
  broom::tidy(mod, seedMods)

#### CORRELATION OF DISTANCE, IMPERV, AND POP DENS ####

enviroData <- read_csv("enviroData.csv") %>%
  left_join(., popMeans %>% select(Population, Distance_to_core))

cor(enviroData$Distance_to_core, enviroData$Imperv)
cor(enviroData$Distance_to_core, enviroData$popDens, use = "complete.obs")

#### FIGURES: MAIN TEXT ####

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


## MULTIVARIATE TRAIT CHANGE ##

## FIGURE 2 ##

# figure 2A #

# Extract RDA1 and PC1 site scores
df_sites  <- data.frame(scores(rdaPop, display = "sites", scaling = "sites")[,1:2]) %>%
  rownames_to_column(var = "Population") %>%
  mutate(Population = seq(1:27)) %>%
  left_join(., popMeans_forRDA %>% select(Population, Distance_to_core))

# Extract RDA1 and PC1 species scores
df2_species  <- data.frame(scores(rdaPop, display = "species", scaling = "species")[,1:2])     # loadings for PC1 and PC2
row.names(df_sites) <- seq(1:27)

# Plot site scores along first 2 axes. Colour points by distance
rda_plot <- ggplot(df_sites, aes(x = RDA1, y = PC1)) + 
  geom_point(size = 4, shape = 21, colour = "black", aes(fill = Distance_to_core)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient(low = "white", high = "black",
                      limits = c(0, 50), breaks = seq(0, 50, 10)) +
  guides(fill = guide_colorbar(ticks = FALSE,
                               barwidth = 1.5, 
                               barheight = 15)) +
  ng1 + theme(legend.position = "right", 
              legend.direction="vertical",
              legend.text = element_text(size=15), 
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm")) 

# Add arrows to RDA plot. 
rda_triplot <- rda_plot +
  geom_segment(data = df2_species, aes(x = 0, xend = RDA1, y=0, yend = PC1), 
               color = "black", size = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(data = df2_species, 
            aes(x = RDA1, y = PC1, label = rownames(df2_species),
                hjust = 0.5 * (1 - sign(RDA1)), vjust = 0.5 * (1 - sign(PC1))), 
            color = "black", size = 2.5) +
  xlab("RDA1 (8.6%)") + ylab("PC1 (35%)")
rda_triplot

ggsave("analysis/figures/main-text/figure2A_RDA-triplot.pdf", 
       plot = rda_triplot, width = 8, height = 6, unit = "in", dpi = 600)

# figure 2B #

clineMax_plot <- ggplot(popMeans, aes(x = Distance_to_core, y = clinemax)) +
  geom_point(size = 4, colour = "black") +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  ng1
clineMax_plot

ggsave("analysis/figures/main-text/figure2B_clineMax_by_distance.pdf", 
       plot = clineMax_plot, width = 8, height = 6, unit = "in", dpi = 600)

## FIGURE 3 ##

# figure 3A #

cols <- c("BB"="#FF0000","HB"="#F2AD00","SB"="#5BBCD6")
linetype <- c("BB"="dashed","HB"="dotted","SB"="dotdash")
plotPoll <-
  ggplot(visitshare.polldata, aes(x = Distance_to_core, y = Visits_per_Inf)) +
  labs(x = "Source population distance\nto urban centre (km)",
       y = "Number of pollinator visits\n per inflorescence") +
  geom_line(data = visitshare.polldata %>% filter(Morph == "Honey Bee"),
            stat = "smooth",
            method = "loess",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "HB",
                colour = "HB")) +
  geom_line(data = visitshare.polldata %>% filter(Morph == "Bumble Bee"),
            stat = "smooth",
            method = "loess",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "BB",
                colour = "BB")) +
  geom_line(data = visitshare.polldata %>% filter(Morph == "Sweat Bee"),
            stat = "smooth",
            method = "loess",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "SB",
                colour = "SB")) +
  geom_smooth(method = "lm", size = 2, colour = "black", se = FALSE) +
  coord_cartesian(xlim = c(2, 47), ylim = c(0, 2.5)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  scale_x_continuous(breaks = seq(0, 40, 10)) +
  scale_color_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                     values = cols) +
  scale_linetype_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                        values = linetype) +
  ng1 + theme(legend.key.size = unit(1.5, "cm"),
              legend.position = "top",
              legend.direction = "horizontal")
plotPoll

ggsave("analysis/figures/main-text/figure3A_visitsPerInf_by_Distance.pdf", 
       plot = plotPoll, width = 8, height = 6, unit = "in", dpi = 600)

# figure 3B #

linetypeSeedFlwrs <- c("Field" = "dashed", "CG" = "dotted")
seedPerFlwr_plot <- ggplot(seedsPerFlwr, aes(x = Distance_to_core, y = Seeds_per_flower)) +
  geom_smooth(method = "lm", size = 2.0, colour = "black", se = FALSE) +
  geom_point(size = 4, alpha = 0.5, colour = "black", aes(shape = source, fill = source)) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("black", "white")) +
  geom_line(data = seedsPerFlwr %>% filter(source == "field"),
            stat = "smooth", 
            method = "lm",
            alpha = 0.5,
            size = 1.5,
            aes(linetype = "Field")) +
  geom_line(data = seedsPerFlwr %>% filter(source == "common garden"),
            stat = "smooth", 
            method = "lm",
            alpha = 0.5,
            size = 1.5,
            aes(linetype = "CG")) +
  scale_linetype_manual(name = "", labels = c("Field", "Common Garden"),
                        values = linetypeSeedFlwrs) +
  ng1 + theme(legend.key.size = unit(1.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"))
seedPerFlwr_plot

ggsave("analysis/figures/main-text/figure3B_seeds_per_flower.pdf", 
       plot = seedPerFlwr_plot, width = 8, height = 6, unit = "in", dpi = 600)

#### SUPPLEMENTARY MATERIALS ####

#### GENETIC VARIATION ACROSS FAMILIES ####

##Function for calculating coefficient of genotypic variance
#x is a variance dataframe and y is trait mean
CVg <- function(x,y) {
  100 * (sqrt(x$vcov[1])/y)
}

# Load in family means dataset
# Create unique family ID column
commonGardenData <- read_csv("data-clean/experimentalData_individualPlants.csv") %>%
  mutate(Family_unique = interaction(Population, Family_ID))

# Run models with family as a random effect

# GERMINATION TIME
germTime <- lmer(Time_to_germination ~ (1|Family_unique), 
                 data = commonGardenData, REML = TRUE)
summary(germTime)
pvalGermTime <- round(ranova(germTime)["Pr(>Chisq)"][2, 1], 3)
varGermTime <- as.data.frame(VarCorr(germTime), comp = "Variance")
meanGermTime <-  mean(commonGardenData$Time_to_germination, na.rm = TRUE)
CvgGermTime <- CVg(varGermTime, meanGermTime)

# DAYS TO FIRST FLOWER
firstFlower <- lmer(Days_to_flower ~ (1|Family_unique),
                   data = commonGardenData, REML = TRUE)

summary(firstFlower)
pvalFirstFlower <- round(ranova(firstFlower)["Pr(>Chisq)"][2, 1], 3)
varFirstFlower <- as.data.frame(VarCorr(firstFlower), comp = "Variance")
meanFirstFlower <-  mean(commonGardenData$Days_to_flower, na.rm = TRUE)
CvgFirstFlower <- CVg(varFirstFlower, meanFirstFlower)

# NUMBER OF INFLORESCENCES
numInf <- lmer(Num_Inf ~ (1|Family_unique), 
               data = commonGardenData, REML = TRUE)

summary(numInf)
pvalNumInf <- round(ranova(numInf)["Pr(>Chisq)"][2, 1], 3)
varNumInf <- as.data.frame(VarCorr(numInf), comp = "Variance")
meanNumInf <-  mean(commonGardenData$Num_Inf, na.rm = TRUE)
CvgNumInf <- CVg(varNumInf, meanNumInf)

# REPRODUCTIVE BIOMASS
repBio <- lmer(Reprod_biomass ~ (1|Family_unique), 
               data = commonGardenData, REML = TRUE)

summary(repBio)
pvalRepBio <- round(ranova(repBio)["Pr(>Chisq)"][2, 1], 3)
varRepBio <- as.data.frame(VarCorr(repBio), comp = "Variance")
meanRepBio <-  mean(commonGardenData$Reprod_biomass, na.rm = TRUE)
CvgRepBio <- CVg(varRepBio, meanRepBio)

# VEGETATIVE BIOMASS
vegBio <- lmer(Veget_biomass ~ (1|Family_unique), 
               data = commonGardenData, REML = TRUE)

summary(vegBio)
pvalVegBio <- round(ranova(vegBio)["Pr(>Chisq)"][2, 1], 3)
varVegBio <- as.data.frame(VarCorr(vegBio), comp = "Variance")
meanVegBio <-  mean(commonGardenData$Veget_biomass, na.rm = TRUE)
CvgVegBio <- CVg(varVegBio, meanVegBio)

# BANNER WIDTH
bnrWdth <- lmer(Avg_bnr_wdth ~ (1|Family_unique), 
               data = commonGardenData, REML = TRUE)

summary(bnrWdth)
pvalBnrWdth <- round(ranova(bnrWdth)["Pr(>Chisq)"][2, 1], 3)
varBnrWdth <- as.data.frame(VarCorr(bnrWdth), comp = "Variance")
meanBnrWdth <-  mean(commonGardenData$Avg_bnr_wdth, na.rm = TRUE)
CvgBnrWdth <- CVg(varBnrWdth, meanBnrWdth)

# BANNER WIDTH
bnrLgth <- lmer(Avg_bnr_lgth ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(bnrLgth)
pvalBnrLgth <- round(ranova(bnrLgth)["Pr(>Chisq)"][2, 1], 3)
varBnrLgth <- as.data.frame(VarCorr(bnrLgth), comp = "Variance")
meanBnrLgth <-  mean(commonGardenData$Avg_bnr_lgth, na.rm = TRUE)
CvgBnrLgth <- CVg(varBnrLgth, meanBnrLgth)

# PETIOLE LENGTH
petLgth <- lmer(Avg_petiole_lgth ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(petLgth)
pvalPetLgth <- round(ranova(petLgth)["Pr(>Chisq)"][2, 1], 3)
varPetLgth <- as.data.frame(VarCorr(petLgth), comp = "Variance")
meanPetLgth <-  mean(commonGardenData$Avg_petiole_lgth, na.rm = TRUE)
CvgPetLgth <- CVg(varPetLgth, meanPetLgth)

# PEDUNCLE LENGTH
pedLgth <- lmer(Avg_peducle_lgth ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(pedLgth)
pvalPedLgth <- round(ranova(pedLgth)["Pr(>Chisq)"][2, 1], 3)
varPedLgth <- as.data.frame(VarCorr(pedLgth), comp = "Variance")
meanPedLgth <-  mean(commonGardenData$Avg_peducle_lgth, na.rm = TRUE)
CvgPedLgth <- CVg(varPedLgth, meanPedLgth)

# LEAF WIDTH
leafWdth <- lmer(Avg_leaf_wdth ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(leafWdth)
pvalLeafWdth <- round(ranova(leafWdth)["Pr(>Chisq)"][2, 1], 3)
varLeafWdth <- as.data.frame(VarCorr(leafWdth), comp = "Variance")
meanLeafWdth <-  mean(commonGardenData$Avg_leaf_wdth, na.rm = TRUE)
CvgLeafWdth <- CVg(varLeafWdth, meanLeafWdth)

# LEAF LENGTH
leafLgth <- lmer(Avg_leaf_lgth ~ (1|Family_unique), 
                 data = commonGardenData, REML = TRUE)

summary(leafLgth)
pvalLeafLgth <- round(ranova(leafLgth)["Pr(>Chisq)"][2, 1], 3)
varLeafLgth <- as.data.frame(VarCorr(leafLgth), comp = "Variance")
meanLeafLgth <-  mean(commonGardenData$Avg_leaf_lgth, na.rm = TRUE)
CvgLeafLgth <- CVg(varLeafLgth, meanLeafLgth)

# STOLON THICKNESS
stolThick <- lmer(Avg_stolon_thick ~ (1|Family_unique), 
                 data = commonGardenData, REML = TRUE)

summary(stolThick)
pvalStolThick <- round(ranova(stolThick)["Pr(>Chisq)"][2, 1], 3)
varStolThick <- as.data.frame(VarCorr(stolThick), comp = "Variance")
meanStolThick <-  mean(commonGardenData$Avg_stolon_thick, na.rm = TRUE)
CvgStolThick <- CVg(varStolThick, meanStolThick)

# NUMBER OF FLOWERS
numFlwr <- lmer(Avg_num_flwrs ~ (1|Family_unique), 
                  data = commonGardenData, REML = TRUE)

summary(numFlwr)
pvalNumFlwr <- round(ranova(numFlwr)["Pr(>Chisq)"][2, 1], 3)
varNumFlwr <- as.data.frame(VarCorr(numFlwr), comp = "Variance")
meanNumFlwr <-  mean(commonGardenData$Avg_num_flwrs, na.rm = TRUE)
CvgNumFlwr <- CVg(varNumFlwr, meanNumFlwr)

# SEX/ASEX
sexAsex <- lmer((Reprod_biomass / Veget_biomass) ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(sexAsex)
pvalSexAsex <- round(ranova(sexAsex)["Pr(>Chisq)"][2, 1], 3)
varSexAsex <- as.data.frame(VarCorr(sexAsex), comp = "Variance")
meanSexAsex <-  mean(commonGardenData$Reprod_biomass / commonGardenData$Veget_biomass, 
                     na.rm = TRUE)
CvgSexAsex <- CVg(varSexAsex, meanSexAsex)

# HCN PRESENCE/ABSENCE
freqHCN <- lmer(HCN_Results ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(freqHCN)
pvalFreqHCN <- round(ranova(freqHCN)["Pr(>Chisq)"][2, 1], 3)
varFreqHCN <- as.data.frame(VarCorr(freqHCN), comp = "Variance")
meanFreqHCN <-  mean(commonGardenData$HCN_Results, 
                     na.rm = TRUE)
CvgFreqHCN <- CVg(varFreqHCN, meanFreqHCN)

# CLINEMAX

# Add clinemax to individuals using species scores from RDA
# Traits are first standardized by dividing by experiment-wide mean.
# Standardized trait values are multiplied by their cannonical regression
# coefficients on RDA. This is done for each trait and then summed. 
commonGardenData <- commonGardenData %>%
  mutate(clinemax = 
           (Time_to_germination / mean(Time_to_germination, na.rm = TRUE)) * species_scores["Germination"] +
           (Days_to_flower / mean(Days_to_flower, na.rm = TRUE)) * species_scores["Days_to_flower"] +
           (Num_Inf / mean(Num_Inf, na.rm = TRUE)) * species_scores["Num_Inf"] +
           (Reprod_biomass / mean(Reprod_biomass, na.rm = TRUE)) * species_scores["Reprod_biomass"] +
           (Veget_biomass / mean(Veget_biomass, na.rm = TRUE)) * species_scores["Veget_biomass"] +
           (Avg_bnr_wdth / mean(Avg_bnr_wdth, na.rm = TRUE)) * species_scores["Bnr_wdth"] +
           (Avg_bnr_lgth / mean(Avg_bnr_lgth, na.rm = TRUE)) * species_scores["Bnr_lgth"] +
           (Avg_petiole_lgth / mean(Avg_petiole_lgth, na.rm = TRUE)) * species_scores["Petiole_lgth"] +
           (Avg_peducle_lgth / mean(Avg_peducle_lgth, na.rm = TRUE)) * species_scores["Peducle_lgth"] +
           (Avg_num_flwrs / mean(Avg_num_flwrs, na.rm = TRUE)) * species_scores["Num_flwrs"] +
           (Avg_leaf_wdth / mean(Avg_leaf_wdth, na.rm = TRUE)) * species_scores["Leaf_wdth"] +
           (Avg_leaf_lgth / mean(Avg_leaf_lgth, na.rm = TRUE)) * species_scores["Leaf_lgth"] +
           (Avg_stolon_thick / mean(Avg_stolon_thick, na.rm = TRUE)) * species_scores["Stolon_thick"] +
           (HCN_Results / mean(HCN_Results, na.rm = TRUE)) * species_scores["FreqHCN"])

clineMax <- lmer(clinemax ~ (1|Family_unique), 
                data = commonGardenData, REML = TRUE)

summary(clineMax)
pvalClineMax <- round(ranova(clineMax)["Pr(>Chisq)"][2, 1], 3)
varClineMax <- as.data.frame(VarCorr(clineMax), comp = "Variance")
meanClineMax <-  mean(commonGardenData$clinemax, 
                     na.rm = TRUE)
# Absolute value since clinemax scores are negative
CvgClineMax <- abs(CVg(varClineMax, meanClineMax))


#### MODELS FOR TRAIT CHANGE WITH URBANIZATION: FAMILY MEANS ####

# Load in family mean data
familyMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")

## MODELS ##

# GERMINATION TIME
germTime_Fam <- lm(Num_Inf ~ Distance_to_core, 
                   data = familyMeans)

# Plot histogram of residuals and resudials vs. fitted.
# Residuals right skewed and resids vs. fitted shows fanning
germTime_resids <- residuals(germTime_Fam)
hist(germTime_resids)
plot(germTime_Fam)
plot(Num_Inf ~ Distance_to_core, data = familyMeans)
abline(germTime_Fam)
# Summarize. Distance to core significant
summary(germTime_Fam)

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
vegBio <- lmer(Veget_biomass ~ Distance_to_core + 
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

#### SUPPLEMENTARY ANALYSES ####

#### RDA AND CLINEMAX FOR FAMILY MEANS ####

#### MULTIVARIATE TRAIT CHANGE WITH URBANIZATION: POPULATION MEANS ####

# Load in family mean dataset
familyMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")

# Subset popMeans datafor use in RDA
familyMeans_forRDA <- familyMeans %>%
  na.omit() %>%
  select(-total_plants, -n_HCN, -Avg_seeds_per_flower) %>%
  select(Population, Distance_to_core, Time_to_germination:freqHCN)

# Perform RDA with multiple traits as response, distance as sole predictors
rdaFam <- rda(familyMeans_forRDA %>% 
                select(Time_to_germination:freqHCN) ~ 
                familyMeans_forRDA$Distance_to_core,
              scale = TRUE)
summary(rdaFam)

# Permutation based test of significance of distance term in RDA
set.seed(42)
anova.cca(rdaFam, by = "term", permutations = 1000)

# Extract canonical coefficients of traits onto RDA1 (i.e. 'species scores')
species_scoresFam <- scores(rdaFam)$species[,"RDA1"]

# Calculate cline_max according to Stock et al. 
# Traits are first standardized by dividing by experiment-wide mean.
# Standardized trait values are multiplied by their cannonical regression
# coefficients on RDA. This is done for each trait and then summed. 
allData <- allData %>%
  mutate(clinemax = 
           (Time_to_germination / mean(Time_to_germination, na.rm = TRUE)) * species_scoresFam["Time_to_germination"] +
           (Days_to_flower / mean(Days_to_flower, na.rm = TRUE)) * species_scoresFam["Days_to_flower"] +
           (Num_Inf / mean(Num_Inf, na.rm = TRUE)) * species_scoresFam["Num_Inf"] +
           (Reprod_biomass / mean(Reprod_biomass, na.rm = TRUE)) * species_scoresFam["Reprod_biomass"] +
           (Veget_biomass / mean(Veget_biomass, na.rm = TRUE)) * species_scoresFam["Veget_biomass"] +
           (Avg_bnr_wdth / mean(Avg_bnr_wdth, na.rm = TRUE)) * species_scoresFam["Avg_bnr_wdth"] +
           (Avg_bnr_lgth / mean(Avg_bnr_lgth, na.rm = TRUE)) * species_scoresFam["Avg_bnr_lgth"] +
           (Avg_petiole_lgth / mean(Avg_petiole_lgth, na.rm = TRUE)) * species_scoresFam["Avg_petiole_lgth"] +
           (Avg_peducle_lgth / mean(Avg_peducle_lgth, na.rm = TRUE)) * species_scoresFam["Avg_peducle_lgth"] +
           (Avg_num_flwrs / mean(Avg_num_flwrs, na.rm = TRUE)) * species_scoresFam["Avg_num_flwrs"] +
           (Avg_leaf_wdth / mean(Avg_leaf_wdth, na.rm = TRUE)) * species_scoresFam["Avg_leaf_wdth"] +
           (Avg_leaf_lgth / mean(Avg_leaf_lgth, na.rm = TRUE)) * species_scoresFam["Avg_leaf_lgth"] +
           (Avg_stolon_thick / mean(Avg_stolon_thick, na.rm = TRUE)) * species_scoresFam["Avg_stolon_thick"] +
           (freqHCN / mean(freqHCN, na.rm = TRUE)) * species_scoresFam["freqHCN"])

# Model testing for cline in cline_max
clineMax_modFam <- lm(clinemax ~ Distance_to_core, data = familyMeans)
summary(clineMax_modFam)
plot(clinemax ~ Distance_to_core, data = familyMeans)
abline(clineMax_modFam)


#Multiple regression analysis on abundance with distance, without controlling for inf #. 
Ab <- lm(Num.Ind ~ Distance_to_core*Morph,
         data = visitshare.polldata)
summary(Ab)
Anova(Ab, type = "III")

#Multiple regression analysis on abundance with distance, controlling for inf #. 
Ab <- lm(Num.Ind ~ Distance_to_core*Morph + Num_Inf,
         data = visitshare.polldata)
summary(Ab)
Anova(Ab, type = "III")

#### FIGURES AND TABLES: SUPPLEMENTAL ####

trait <- c("Banner length", "Banner width", "Clinemax", "Days to first flower",
           "Cyanogenesis", "Germination time", "Leaf length", "Leaf width",
           "Number of flowers", "Number of inflorescences", "Peduncle length",
           "Petiole length", "Reproductive biomass", "Sex/Asex ratio",
           "Stolon thickness", "Vegetative biomass")
means <- c(meanBnrLgth, meanBnrWdth, meanClineMax, meanFirstFlower,
           meanFreqHCN, meanGermTime, meanLeafLgth, meanLeafWdth,
           meanNumFlwr, meanNumInf, meanPedLgth, meanPetLgth, 
           meanRepBio, meanSexAsex, meanStolThick, meanVegBio)
CVGs <- c(CvgBnrLgth, CvgBnrWdth, CvgClineMax, CvgFirstFlower,
          CvgFreqHCN, CvgGermTime, CvgLeafLgth, CvgLeafWdth,
          CvgNumFlwr, CvgNumInf, CvgPedLgth, CvgPetLgth, 
          CvgRepBio, CvgSexAsex, CvgStolThick, CvgVegBio)
pvals <- c(pvalBnrLgth, pvalBnrWdth, pvalClineMax, pvalFirstFlower,
           pvalFreqHCN, pvalGermTime, pvalLeafLgth, pvalLeafWdth,
           pvalNumFlwr, pvalNumInf, pvalPedLgth, pvalPetLgth, 
           pvalRepBio, pvalSexAsex, pvalStolThick, pvalVegBio) / 2

dfCVG <- as.data.frame(cbind(trait, means, CVGs, pvals))

write_csv(dfCVG, "analysis/tables/geneticVar.csv")
