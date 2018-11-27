# ken a thompson analysis
# sex in city

# load packages
library(psych)
library(tidyverse)
library(vegan)

#load data

seeddata.field.popmean <- read.csv('data-clean/seedFlwrRatio_fieldPlants_cleaned.csv')
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

# generate data frame of FAMILY_MEANS
familymeans.cg <- indplant.cg %>% 
  group_by(Population, Family_ID) %>% 
  summarise_all(funs(mean(., na.rm = TRUE)))

# multivariate cline analysis
# trying the CCA method stinch used

Data.MatrixTrait <- indplant.cg[,c(5, 12, 14)]

sc <- rda(Data.MatrixTrait ~ Distance_to_core, data = indplant.cg, scale = T)

# str(sc)
# 
# sc$terminfo


plot(sc)

# get phenotype variables
vars.for.cormat <- familymeans.cg %>% 
  select(Time_to_germination, Days_to_flower:Avg_stolon_thick)
    
# general PCA for dataframe
experimentPCA <- principal(vars.for.cormat, nfactors = 5)

# variance explained
experimentPCA$Vaccounted

#loadings
experimentPCA$loadings

PC.scores <- as.data.frame(experimentPCA$scores)

PC.scores$Distance_to_core <- familymeans.cg$Distance_to_core


# regression of PCs

expt.data.all$PC1 <- PC.scores$RC1
expt.data.all$PC3 <- PC.scores$RC3


# probability of flowering change along the gradient
mylogit <- glm(did.flower ~ Distance_to_core, data = indplant.cg2, family = "binomial")
summary(mylogit)


# pollinator observations
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

PC.gg <- ggplot(PC.scores, aes(x = Distance_to_core, y = RC4)) +
    geom_point(alpha = 1) +
    labs(x = "source population distance\nto urban centre (km)",
    y = "vegetative biomass (g)") +
    geom_smooth(method = lm, se = T) +
  theme_KT
PC.gg

Stolon.Thick <- ggplot(indplant.cg2, aes(x = Distance_to_core, y = did.flower)) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "reprod_bm") +
  geom_smooth(method = lm, se = T) +
  theme_KT
Stolon.Thick



summary(lm((Reprod_biomass) ~ Distance_to_core, data = popmean.cg))

SeedFlwr.Ratio.Field <- ggplot(seeddata.field.popmean, aes(x = PolVis, y = Seeds_per_flower)) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "number of seeds per flower") +
  geom_smooth(method = lm, se = T) +
  theme_KT
SeedFlwr.Ratio.Field
  