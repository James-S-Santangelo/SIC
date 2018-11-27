# ken a thompson analysis

# ggplot themes


# load packages
library(tidyverse)

#load data
polldata <- read.csv('data-clean/pollinatorObservations.csv')
seeddata.field.popmean <- read.csv('data-clean/seedFlowerRatio_fieldData.csv')
indplant.cg <- read.csv('data-clean/experimentalData_individualPlants.csv')
popmean.cg <- read.csv('data-clean/experimentalData_popMeans.csv')

# plot themes

theme_KT <- theme(aspect.ratio=1.0,panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border=element_blank(),
                  axis.line = element_line(size=1),
                  axis.line.x = element_line(color="black", size = 1),
                  axis.line.y = element_line(color="black", size = 1),
                  axis.ticks=element_line(color="black"),
                  axis.text=element_text(color="black"),
                  axis.title=element_text(color="black"),
                  axis.title.y=element_text(size = 10, family = "Helvetica"),
                  axis.title.x=element_text(size = 10, family = "Helvetica"),
                  axis.text.x=element_text(size = 8, family = "Helvetica"),
                  axis.text.y=element_text(size = 8, family = "Helvetica"),
                  legend.position="none",
                  legend.title=element_blank(),
                  plot.title = element_text(hjust=0))


# analysis of pollinator observations
# where are the plants getting pollinated more?

popmeans.polldata <- polldata %>% 
  group_by(Population, Distance_to_core) %>% 
  summarize(MeanInf = mean(Num_Inf),
            MeanVisit = mean(Num_Visit)) %>% 
  mutate(Vis_per_Inf = MeanVisit/MeanInf)

visitshare.polldata <- polldata %>% 
  group_by(Population) %>% 
  #how many visits total per population
  mutate(SumVisit = sum(Num_Visit)) %>% 
  # group all by 'morph' of pollinator
  group_by(Morph, Population, SumVisit, Distance_to_core) %>% 
  summarize(VisitsPerMorph = sum(Num_Visit)) %>% 
  mutate(PropVisits = VisitsPerMorph / SumVisit)

indplant.cg2 <- indplant.cg %>% 
  mutate(did.flower = ifelse(Reprod_biomass > 0, 1, 0))


mylogit <- glmer(did.flower ~ Distance_to_core + (1|Population), data = indplant.cg2, family = "binomial")
summary(mylogit)



visitshare.polldata.bumble <- visitshare.polldata %>% 
  filter(Morph %in% c("Honey Bee", "Bumble Bee"))

merge(popmeans.polldata, visitshare.polldata.bumble, by = "Distance_to_core")

# analysis of field seed data

seeddata.field.2 <- seeddata.field %>% 
  mutate(SeedsPerFlower = Num.Seeds / Num.Flwrs)

# merge(seeddata.field.2, )

# add number of pollinator vists per infl to chart

seeddata.field.popmean$PolVis <- popmeans.polldata$MeanVisit

# general ggplot

Veget.Biom <- ggplot(indplant.cg, aes(x = Distance_to_core, y = Veget_biomass)) +
    geom_point(alpha = 1) +
    labs(x = "source population distance\nto urban centre (km)",
    y = "vegetative biomass (g)") +
    geom_smooth(method = lm, se = T) +
  theme_KT
Veget.Biom

Stolon.Thick <- ggplot(indplant.cg2, aes(x = Distance_to_core, y = did.flower)) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "reprod_bm") +
  geom_smooth(method = lm, se = T) +
  theme_KT
Stolon.Thick



summary(lm((Reprod_biomass) ~ Distance_to_core, data = popmean.cg))

# 1] "Population"           "Family_ID"            "Seed"                 "label"                "Time_to_germination"  "Row"                  "Column"               "Days_to_flower"      
# [9] "Num_Inf"              "HCN_Results"          "Reprod_biomass"       "Veget_biomass"        "Avg_bnr_wdth"         "Avg_bnr_lgth"         "Avg_petiole_lgth"     "Avg_peducle_lgth"    
# [17] "Avg_num_flwrs"        "Avg_leaf_wdth"        "Avg_leaf_lgth"        "Avg_stolon_thick"     "Avg_seeds_per_flower" "Latitude"             "Longitude"            "Distance_to_core"    
# [25] "Distance_to_cg"      
# > Veget.Biom <- ggplot(i

SeedFlwr.Ratio.Field <- ggplot(seeddata.field.popmean, aes(x = PolVis, y = Seeds_per_flower)) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "number of seeds per flower") +
  geom_smooth(method = lm, se = T) +
  theme_KT
SeedFlwr.Ratio.Field
  