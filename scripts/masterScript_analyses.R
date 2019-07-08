

## Required by other scripts. Will eventually be migrated.
# install.packages("psych")
# install.packages("cowplot")
# install.packages("RcppEigen")
# install.packages("mnormt")
# install.packages("BH")
# install.packages("plogr")

# # Commands used to initialize packrat. 
# packrat::init(restart = TRUE, enter = FALSE, infer.dependencies = FALSE)
# packrat::on()
# 
# # install packages (these packages are loaded in this script)
# install.packages("lme4", dependencies = TRUE)
# install.packages("lmerTest", dependencies = TRUE)
# install.packages("tidyverse", dependencies = TRUE)
# install.packages("vegan", dependencies = TRUE)
# install.packages("car", dependencies = TRUE)
# 
# # .libPaths()
# packrat::status()
# packrat::snapshot(infer.dependencies = FALSE)

# Load require packages
library(lme4)
library(lmerTest)
library(tidyverse)
library(vegan)
library(car)

#### MULTIVARIATE TRAIT CHANGE WITH URBANIZATION: FAMILY MEANS ####

# Load in family mean dataset
famMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")

# Subset data for RDA
famMeans_forRDA <- famMeans %>%
  # select(-Seeds_per_flower, -Num_Cyano) %>%
  select(-total_plants, -n_HCN, -Avg_seeds_per_flower) %>%
  select(Population, Distance_to_core, Family_ID, Time_to_germination_C:Avg_stolon_thick_C, freqHCN) %>% 
  na.omit()

# Perform RDA with multiple traits as response, distance as sole predictors
rdaFam <- rda(famMeans_forRDA %>% 
                     select(Time_to_germination_C:freqHCN) ~ 
                famMeans_forRDA$Distance_to_core, 
              scale = TRUE, na.action = "na.omit")
summary(rdaFam)

# Permutation based test of significance of distance term in RDA
set.seed(42)
anova.cca(rdaFam, by = "term", permutations = 10000)

# Extract canonical coefficients of traits onto RDA1 (i.e. 'species scores')
species_scores <- scores(rdaFam)$species[,"RDA1"]

# Calculate cline_max according to Stock et al. 
# Traits are first standardized by dividing by experiment-wide mean.
# Standardized trait values are multiplied by their cannonical regression
# coefficients on RDA. This is done for each trait and then summed. 
indPlantData <- read_csv("data-clean/experimentalData_individualPlants.csv")
indPlantData <- indPlantData %>% 
  mutate(clinemax = 
           (Time_to_germination_C * species_scores["Time_to_germination_C"]) +
           (Days_to_flower_C * species_scores["Days_to_flower_C"]) +
           (Num_Inf_C * species_scores["Num_Inf_C"]) +
           (Reprod_biomass_C * species_scores["Reprod_biomass_C"]) +
           (Veget_biomass_C * species_scores["Veget_biomass_C"]) +
           (Avg_bnr_wdth_C * species_scores["Avg_bnr_wdth_C"]) +
           (Avg_bnr_lgth_C * species_scores["Avg_bnr_lgth_C"]) +
           (Avg_petiole_lgth_C * species_scores["Avg_petiole_lgth_C"]) +
           (Avg_peducle_lgth_C * species_scores["Avg_peducle_lgth_C"]) +
           (Avg_num_flwrs_C * species_scores["Avg_num_flwrs_C"]) +
           (Avg_leaf_wdth_C * species_scores["Avg_leaf_wdth_C"]) +
           (Time_to_germination_C * species_scores["Avg_leaf_lgth_C"]) +
           (Avg_stolon_thick_C * species_scores["Avg_stolon_thick_C"]) +
           (HCN_Results * species_scores["freqHCN"]))
  
famMeans_clineMax <- indPlantData %>% 
  group_by(Population, Family_ID, Distance_to_core) %>% 
  summarise(clinemax = mean(clinemax, na.rm = TRUE))

# Model testing for cline in cline_max
clineMax_mod <- lm(clinemax ~ Distance_to_core, data = famMeans_clineMax)
summary(clineMax_mod)

#### UNIVARIATE TRAIT CHANGE WITH URBANIZATION: FAMILY MEANS ####

## SIGNIFICANT ##

# Sig. 
germMod <- lm(Time_to_germination ~ Distance_to_core, data = famMeans)
summary(germMod)
plot(germMod)
hist(residuals(germMod)) # No transformation

# Sig.
ffMod <- lm(Days_to_flower ~ Distance_to_core, data = famMeans)
summary(ffMod)
plot(ffMod)
hist(residuals(ffMod)) # No transformation needed

# Sig.
vegMod <- lm(Veget_biomass ~ Distance_to_core, data = famMeans)
summary(vegMod)
plot(vegMod)
hist(residuals(vegMod)) # No transformation needed

# Sig.
blMod <- lm(Avg_bnr_lgth ~ Distance_to_core, data = famMeans)
summary(blMod)
plot(blMod)
hist(residuals(blMod)) # No transformation needed

# Sig
stMod <- lm(Avg_stolon_thick ~ Distance_to_core, data = famMeans)
summary(stMod)
plot(stMod)
hist(residuals(stMod)) # No transformation needed

# Marg
HCNMod <- lm(freqHCN ~ Distance_to_core, data = famMeans)
summary(HCNMod)
plot(HCNMod)
hist(residuals(HCNMod)) # No transformation

## NOT SIGNIFICANT ##

# NS
infMod <- lm(sqrt(Num_Inf) ~ Distance_to_core, data = famMeans)
summary(infMod)
plot(infMod)
hist(residuals(infMod)) # Square root transformation

# NS
repMod <- lm(Reprod_biomass ~ Distance_to_core, data = famMeans)
summary(repMod)
plot(repMod)
hist(residuals(repMod)) # No transformation needed

# Marg
bwMod <- lm(Avg_bnr_wdth ~ Distance_to_core, data = famMeans)
summary(bwMod)
plot(bwMod)
hist(residuals(bwMod)) #  No transformation needed

# NS
pedMod <- lm(Avg_peducle_lgth ~ Distance_to_core, data = famMeans)
summary(pedMod)
plot(pedMod)
hist(residuals(pedMod)) # No transformation needed

# NS
numFlwrsMod <- lm(Avg_num_flwrs ~ Distance_to_core, data = famMeans)
summary(numFlwrsMod)
plot(numFlwrsMod)
hist(residuals(numFlwrsMod)) # No transformation needed

# NS
leafWdthMod <- lm(Avg_leaf_wdth ~ Distance_to_core, data = famMeans)
summary(leafWdthMod)
plot(leafWdthMod)
hist(residuals(leafWdthMod)) # No transformation needed

# NS
leafLgthMod <- lm(Avg_leaf_lgth ~ Distance_to_core, data = famMeans)
summary(leafLgthMod)
plot(leafLgthMod)
hist(residuals(leafLgthMod)) # No transformation needed

# NS
petMod <- lm(Avg_petiole_lgth ~ Distance_to_core, data = famMeans)
summary(petMod)
plot(petMod)
hist(residuals(petMod)) # No transformation needed

# Marg
sexInvestMod <- lm(sex_asex ~ Distance_to_core, data = famMeans)
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

# Means and SDs for number of visits per inflorescence by morph
visitshare_polldata %>% 
  group_by(Morph) %>% 
  summarise(mean = mean(Visits_per_Inf),
            sd = sd(Visits_per_Inf))

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
popMeans <- read_csv("data-clean/experimentalData_popMeans.csv")
seedFlwrFieldData <- read_csv("data-clean/flwrSeedRatio_fieldPlants.csv") %>%
  select(-X1, -Comments) %>%
  group_by(Population, Distance_to_core) %>%
  summarize(Seeds_per_flower = mean(Seeds_per_flower)) 

# Retrieve seeds per flower from common garden plants
seedFlwrComGarData <- popMeans %>%
  select(Population, Distance_to_core, Avg_seeds_per_flower) %>% 
  rename("Seeds_per_flower" = "Avg_seeds_per_flower")

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

enviroData <- read_csv("data-clean/enviroData.csv") %>%
  left_join(., popMeans %>% select(Population, Distance_to_core))

cor(enviroData$Distance_to_core, enviroData$Imperv, use = "complete.obs")
cor(enviroData$Distance_to_core, enviroData$popDens, use = "complete.obs")

summary(lm(Imperv ~ Distance_to_core, data = enviroData))
summary(lm(popDens ~ Distance_to_core, data = enviroData))

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
df_sites  <- data.frame(scores(rdaFam, display = "sites", scaling = "sites")[,1:2]) %>%
  # rownames_to_column(var = "Population") %>%
  # mutate(Population = seq(1:27)) %>%
  cbind(., famMeans_forRDA %>% select(Population, Distance_to_core))

# Extract RDA1 and PC1 species scores
df2_species  <- data.frame(scores(rdaFam, display = "species", scaling = "species")[,1:2])     # loadings for PC1 and PC2
row.names(df_sites) <- seq(1:27)

# Plot site scores along first 2 axes. Colour points by distance
rda_plot <- ggplot(df_sites, aes(x = RDA1, y = PC1)) + 
  geom_point(size = 3, shape = 21, colour = "black", aes(fill = Distance_to_core)) +
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
  xlab("RDA1 (2.7%)") + ylab("PC1 (26%)")
rda_triplot

ggsave("analysis/figures/main-text/figure2A_RDA-triplot.pdf", 
       plot = rda_triplot, width = 8, height = 6, unit = "in", dpi = 600)

# figure 2B #

clineMax_plot <- ggplot(famMeans_clineMax, aes(x = Distance_to_core, y = clinemax)) +
  geom_point(size = 3, colour = "black") +
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
  ggplot(visitshare_polldata, aes(x = Distance_to_core, y = Visits_per_Inf)) +
  labs(x = "Source population distance\nto urban centre (km)",
       y = "Number of pollinator visits\n per inflorescence") +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Honey Bee"),
            stat = "smooth",
            method = "loess",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "HB",
                colour = "HB")) +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Bumble Bee"),
            stat = "smooth",
            method = "loess",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "BB",
                colour = "BB")) +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Sweat Bee"),
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

#### FIGURES: SUPPLEMENTARY MATERIALS ####

## FIGURE S1

# Figure S1A: Imperv. vs. distance
imperv_vs_distance <- ggplot(enviroData, aes(x = Distance_to_core, y = Imperv)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', colour = 'black', se = FALSE) + 
  ylab('Impervious surface (%)') + xlab('Source population distance from urban core (km)') +
  ng1
imperv_vs_distance

ggsave("analysis/figures/sup-mat/figureS1a_Imperv_vs.distance.pdf", 
       plot = imperv_vs_distance, width = 6, height = 6, unit = "in", dpi = 600)

# Figure S1B: Pop Dens vs. distance
popDens_vs_distance <- ggplot(enviroData, aes(x = Distance_to_core, y = popDens)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', colour = 'black', se = FALSE) + 
  ylab(expression(~Human~population~density~(per~km^2))) + xlab('Source population distance from urban core (km)') +
  ng1
popDens_vs_distance

ggsave("analysis/figures/sup-mat/figureS1b_popDens_vs.distance.pdf", 
       plot = popDens_vs_distance, width = 6, height = 6, unit = "in", dpi = 600)


## FIGURE S2

#' Generates biplot of with response variable against predictor variable
#'     Writes biplot to disk in outpath.
#'
#' @param df Dataframe containing variables that will be plotted as columns
#' @param response_var Variable to be plotted on y-axis
#' @param outpath Path to which plot will be written
#' @param figID Figure numbe and letter (e.g. Figure 2a)
#' 
#' @return None. Writes plot to disk.
create_Biplot <- function(df, response_var, outpath, figID){
  
  # print(path)
  response_vector <- df %>% pull(response_var)
  
  plot <- df %>%
    ggplot(., aes_string(x = "Distance_to_core", y = response_var)) +
    geom_point(colour = "black", size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", size = 1) + 
    ylab(response_var) + xlab("Source population distance from the urban core (km)") +
    ng1
    
    path <- paste0(outpath, figID, "_", response_var, "_by_distance", ".pdf")
    
    # Write dataframe
    ggsave(filename = path, plot = plot, device = "pdf",
           width = 6, height = 6, dpi = 300)
  
}


outpath <- "analysis/figures/sup-mat/"

# Figure S2a
create_Biplot(df = famMeans, response_var = "Time_to_germination", 
              outpath = outpath, figID = "figureS2a")
# Figure S2b
create_Biplot(df = famMeans, response_var = "Days_to_flower", 
              outpath = outpath, figID = "figureS2b")
# Figure S2c
create_Biplot(df = famMeans, response_var = "Veget_biomass", 
              outpath = outpath, figID = "figureS2c")
# Figure S2d
create_Biplot(df = famMeans, response_var = "Avg_bnr_lgth", 
              outpath = outpath, figID = "figureS2d")
# Figure S2e
create_Biplot(df = famMeans, response_var = "Avg_stolon_thick", 
              outpath = outpath, figID = "figureS2e")
# Figure S2f
create_Biplot(df = famMeans, response_var = "freqHCN", 
              outpath = outpath, figID = "figureS2f")

## FIGURE S3

# Figure S3a
create_Biplot(df = famMeans, response_var = "Num_Inf", 
              outpath = outpath, figID = "figureS3a")
# Figure S3b
create_Biplot(df = famMeans, response_var = "Reprod_biomass", 
              outpath = outpath, figID = "figureS3b")
# Figure S3c
create_Biplot(df = famMeans, response_var = "Avg_bnr_wdth", 
              outpath = outpath, figID = "figureS3c")
# Figure S3d
create_Biplot(df = famMeans, response_var = "Avg_peducle_lgth", 
              outpath = outpath, figID = "figureS3d")
# Figure S3e
create_Biplot(df = famMeans, response_var = "Avg_num_flwrs", 
              outpath = outpath, figID = "figureS3e")
# Figure S3f
create_Biplot(df = famMeans, response_var = "Avg_leaf_wdth", 
              outpath = outpath, figID = "figureS3f")
# Figure S3g
create_Biplot(df = famMeans, response_var = "Avg_leaf_lgth", 
              outpath = outpath, figID = "figureS3g")
# Figure S3h
create_Biplot(df = famMeans, response_var = "Avg_petiole_lgth", 
              outpath = outpath, figID = "figureS3h")
# Figure S3i
create_Biplot(df = famMeans, response_var = "sex_asex", 
              outpath = outpath, figID = "figureS3i")


## FIGURE S4

cols <- c("BB"="#FF0000","HB"="#F2AD00","SB"="#5BBCD6")
linetype <- c("BB"="dashed","HB"="dotted","SB"="dotdash")
shape <- c("BB"=21,"HB"=22,"SB"=24)
fill <- c("BB"="#FF0000","HB"="#F2AD00","SB"="#5BBCD6")
plotPoll_lin <-
  ggplot(visitshare_polldata, aes(x = Distance_to_core, y = Visits_per_Inf)) +
  labs(x = "Source population distance\nto urban centre (km)",
       y = "Number of pollinator visits\n per inflorescence") +
  geom_point(data = visitshare_polldata %>% filter(Morph == "Honey Bee"),
             size = 2.5, aes(colour = "HB", shape = "HB", fill = "HB")) +  
  geom_point(data = visitshare_polldata %>% filter(Morph == "Bumble Bee"),
             size = 2.5, aes(colour = "BB", shape = "BB", fill = "BB")) + 
  geom_point(data = visitshare_polldata %>% filter(Morph == "Sweat Bee"),
             size = 2.5, aes(colour = "SB", shape = "SB", fill = "SB")) +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Honey Bee"),
            stat = "smooth",
            method = "lm",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "HB",
                colour = "HB")) +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Bumble Bee"),
            stat = "smooth",
            method = "lm",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "BB",
                colour = "BB")) +
  geom_line(data = visitshare_polldata %>% filter(Morph == "Sweat Bee"),
            stat = "smooth",
            method = "lm",
            se = F, 
            # alpha = 0.7,
            size = 1.25,
            aes(linetype = "SB",
                colour = "SB")) +
  geom_smooth(method = "lm", size = 2, colour = "black", se = FALSE) +
  # coord_cartesian(xlim = c(2, 47), ylim = c(0, 2.5)) +
  # scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  scale_x_continuous(breaks = seq(0, 40, 10)) +
  scale_color_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                     values = cols) +
  scale_linetype_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                        values = linetype) +
  scale_fill_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                        values = fill) +
  scale_shape_manual(name = "", labels = c("Bumble bee", "Honey bee", "Sweat bee"),
                        values = shape) +
  ng1 + theme(legend.key.size = unit(1.5, "cm"),
              legend.position = "top",
              legend.direction = "horizontal")
plotPoll_lin

ggsave("analysis/figures/sup-mat/figureS4_visitsPerInf_by_Distance_linear.pdf", 
       plot = plotPoll_lin, width = 8, height = 6, unit = "in", dpi = 600)


#### TABLES: SUPPLEMENTARY MATERIALS ####

#' Creates vector with trait mean and output from standardized and 
#'     unstandardized regressions
#'
#' @param df Dataframe containing variables for regression
#' @param trait Trait to use as response variable in regression
#' 
#' @return model_vector. Trait meanm standardized and unstandardized 
#'     model betas, p-value, r squared
summarise_trait_regressions <- function(df, trait){
  
  # Retrieve response and predictors vars
  response_var <- df %>% pull(trait)
  predictor_var <- df %>% pull(Distance_to_core)
  
  # Calculate trait mean
  trait_mean <- round(df %>% pull(trait) %>% mean(., na.rm = TRUE), 4)
  
  # Run unstandardized model
  unstandard_mod <- lm(response_var ~ predictor_var)
  
  # Pull relevant coefficients
  unstand_beta <- round(summary(unstandard_mod)$coef[2, "Estimate"], 4)
  pval <- round(summary(unstandard_mod)$coef[2, "Pr(>|t|)"], 4)
  r_squared <- round(summary(unstandard_mod)$r.squared, 4)
  
  # Run regression with mean-standardized traits for all traits except HCN
  if(!trait == "freqHCN"){
    
    stand_trait = paste0(trait, "_C")
    stand_response_var <- df %>% pull(stand_trait)
    stand_mod <- lm(stand_response_var ~ predictor_var)
    stand_beta <- round(summary(stand_mod)$coef[2, "Estimate"], 4)

    model_vector <- c(trait, trait_mean, unstand_beta, stand_beta, pval, r_squared)
    
  }else{
    
    model_vector <- c(trait, trait_mean, unstand_beta, "NA", pval, r_squared)
    
  }
  
  return(model_vector)
  
}

germTimeVector <- summarise_trait_regressions(famMeans, "Time_to_germination")
FFVector <- summarise_trait_regressions(famMeans, "Days_to_flower")
vegBioVector <- summarise_trait_regressions(famMeans, "Veget_biomass")
bnrLVector <- summarise_trait_regressions(famMeans, "Avg_bnr_lgth")
stolVector <- summarise_trait_regressions(famMeans, "Avg_stolon_thick")
HCNVector <- summarise_trait_regressions(famMeans, "freqHCN")
numInfVector <- summarise_trait_regressions(famMeans, "Num_Inf")
repBioVector <- summarise_trait_regressions(famMeans, "Reprod_biomass")
bnrWVector <- summarise_trait_regressions(famMeans, "Avg_bnr_wdth")
pedVector <- summarise_trait_regressions(famMeans, "Avg_peducle_lgth")
numFlwrVector <- summarise_trait_regressions(famMeans, "Avg_num_flwrs")
leafWVector <- summarise_trait_regressions(famMeans, "Avg_leaf_wdth")
leafLVector <- summarise_trait_regressions(famMeans, "Avg_leaf_lgth")
petVector <- summarise_trait_regressions(famMeans, "Avg_petiole_lgth")

header <- c("Trait", "Mean", "Unstandarized beta", "Standardized beta", "P-value", "R squared")
tableS1 <- rbind(germTimeVector, FFVector, vegBioVector, bnrLVector, stolVector, HCNVector,
                 numInfVector, repBioVector, bnrWVector, pedVector, numFlwrVector, leafWVector,
                 leafLVector, petVector) %>% 
  as.data.frame() %>% 
  setNames(., header)

write_csv(tableS1, "analysis/tables/tableS1_traitRegSummary.csv", col_names = TRUE)


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
