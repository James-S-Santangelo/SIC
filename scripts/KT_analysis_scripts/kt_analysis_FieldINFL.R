# ken a thompson analysis
# sex in city

# load packages
library(lme4)
library(lmerTest)
library(psych)
library(tidyverse)
library(vegan)

#load data

pop.info <- read.csv('data-clean/experimentalData_popMeans.csv')

seeddata.field <- read.csv('data-clean/seedFlwrRatio_fieldPlants_cleaned.csv')

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

# make a dataframe with population number and dist to core to join with other

# stats of individual plants
# linear model
seeds.per.flwr.lm <- lm(Seeds_per_flower ~ Distance_to_core, data = seeddata.field.popmeans)
summary(seeds.per.flwr.lm)

seeds.per.flwr.lmer <- lmer(sqrt(Seeds_per_flower) ~ Distance_to_core + (1|Population), data = seeddata.field)
anova(seeds.per.flwr.lmer, ddf = "Kenward-Roger")
# check for heterogeneity of variance
# plot(seeds.per.flwr.lmer)

# generate a 'pop means' dataset
seeddata.field.popmeans <- seeddata.field %>% 
  group_by(Population, Distance_to_core) %>% 
  summarise(Seeds_per_flower = mean(Seeds_per_flower))

# ggplot
Seeds.Per.Flower.Plot.Field <-
  ggplot(seeddata.field, aes(x = Distance_to_core, y = sqrt(Seeds_per_flower))) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "√(number of seeds per flower)") +
  geom_smooth(method = lm, se = T) +
  theme_KT
Seeds.Per.Flower.Plot.Field

ggsave(filename = "analysis/kt_figures/genetic_variance/seeds_per_flower_field.pdf", plot = Seeds.Per.Flower.Plot.Field, height = 3, width = 3)
