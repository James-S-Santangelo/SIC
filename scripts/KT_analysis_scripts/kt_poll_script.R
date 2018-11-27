# ken thompson SIC pollinator obs script

# load data
polldata <- read.csv('data-clean/pollinatorObservations.csv')

# load packages
library(tidyverse)



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
