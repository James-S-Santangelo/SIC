# ken thompson SIC pollinator obs script

# load data
polldata <- read.csv('data-clean/pollinatorObservations.csv')

polldata.popmean.dist <- polldata %>% 
  group_by(Population) %>% 
  summarise(Distance_to_core = mean(Distance_to_core))

# load packages
library(tidyverse)



# analysis of pollinator observations
# where are the plants getting pollinated more?
#try to impute zeros for unobserved polls
complete.polldata <- polldata %>% 
  group_by(Population,Morph, Distance_to_core) %>% 
  complete(Population, Morph, fill = list(Num_Visit = 0)) 
# summarise(count = sum(z))
# 
# popmeans.polldata <- polldata %>%
#   group_by(Population, Distance_to_core) %>%
#   summarize(MeanInf = mean(Num_Inf),
#             MeanVisit = mean(Num_Visit)) %>%
#   mutate(Vis_per_Inf = MeanVisit/MeanInf)

visitshare.polldata <- complete.polldata %>% 
  group_by(Population, Distance_to_core, Morph) %>% 
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>%  
  mutate(Visits_per_Inf = Num_Visit / Num_Inf) %>% 
  # replace NA in visits/inf with zero
  mutate(Visits_per_Inf = replace_na(Visits_per_Inf, 0))

popmeans.polldata <- polldata %>% 
  group_by(Population, Distance_to_core) %>% 
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>% 
  mutate(Visits_per_Inf = Num_Visit / Num_Inf)

visitshare.polldata.bumble <- visitshare.polldata %>% 
  filter(Morph == "Bumble Bee") 
visitshare.polldata.bumble$Distance_to_core <- polldata.popmean.dist$Distance_to_core

visitshare.polldata.honey <- visitshare.polldata %>% 
  filter(Morph == "Honey Bee") 
visitshare.polldata.honey$Distance_to_core <- polldata.popmean.dist$Distance_to_core

# some basic stats
test.int <- aov(Visits_per_Inf ~ Distance_to_core*Morph, data =visitshare.polldata)
car::Anova(test.int, type = 3)

summary(lm(Visits_per_Inf ~ Distance_to_core, data =popmeans.polldata))


#plots

plot.bumble <-
  ggplot(visitshare.polldata.bumble, aes(x = Distance_to_core, y = Num_Visit)) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "number of Bumblebee visits") +
  geom_smooth(method = loess, se = T) +
  # ylim(c(0, 80)) +
  theme_KT
plot.bumble



plot.allpoll <-
  ggplot(visitshare.polldata, aes(x = Distance_to_core, y = Num_Visit, colour = factor(Morph))) +
  geom_point(alpha = 1) +
  labs(x = "source population distance\nto urban centre (km)",
       y = "number of pollinator visits") +
  geom_smooth(method = loess, se = F) +
  # ylim(c(0, 80)) +
  theme_KT
plot.allpoll
