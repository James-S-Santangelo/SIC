setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/SIC - Sex in the city/*SIC")

# ken thompson SIC pollinator obs script 
# load data
polldata <- read.csv('data-clean/pollinatorObservations.csv')
# load packages
library(tidyverse)

# analysis of pollinator observations
#impute zeros for unobserved polls
complete.polldata <- polldata %>% 
  group_by(Population,Morph, Distance_to_core) %>% 
  complete(Population, Morph, fill = list(Num_Visit = 0)) 

#Summarize data by pollinator functional group
#Used for visit/infl analysis
visitshare.polldata <- complete.polldata %>% 
  group_by(Population, Distance_to_core, Morph) %>% 
  # calculate number of individuals observed per site per morph; number of inflorescence, number of visit per infl.
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>%  
  mutate(Visits_per_Inf = Num_Visit / Num_Inf) %>% 
  # replace NA in visits/inf with zero
  mutate(Visits_per_Inf = replace_na(Visits_per_Inf, 0))

#Summarize daya by population only (merge morphs)
#Used for total abundance analysis
popmeans.polldata <- polldata %>% 
  group_by(Population, Distance_to_core) %>% 
  summarise(Num.Ind = n(),
            Num_Visit = sum(Num_Visit),
            Num_Inf = mean(Num_Inf)) %>% 
  mutate(Visits_per_Inf = Num_Visit / Num_Inf)


#Stats
library(car)

#Visualize data: log transform vistis/inf to achieve normality
hist(log(visitshare.polldata$Visits_per_Inf))
log.vis.inf <- log(visitshare.polldata$Visits_per_Inf+1)
hist((popmeans.polldata$Num.Ind))

#Multiple regression analysis on log(visit/inf) with distance, morph, and their interaction. 
VI <- lm(log.vis.inf ~ Distance_to_core + Morph + Distance_to_core:Morph,
         data = visitshare.polldata)
summary(VI)
Anova(VI, type = "III") #Type 3 for interpreting interaction

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

#Plots
#Load theme
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
          legend.position = "right", legend.direction="vertical",
          legend.text=element_text(size=15), legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=15),legend.key.size = unit(1.0, "cm"))
#Filter by morph (bumblebee(b) and honeybee(h))
morph_b <- filter(visitshare.polldata, Morph == "Bumble Bee")
morph_h <- filter(visitshare.polldata, Morph == "Honey Bee")

#Main effect of distance on pollinator visitation rate
plot <-
  ggplot(visitshare.polldata, aes(x = Distance_to_core, y = Visits_per_Inf)) +
  geom_point(alpha = 1) +
  labs(x = "Source population distance\nto urban centre (km)",
       y = "Number of pollinator visits\n per inflorescence") +
  geom_smooth(method=lm, se=F, fullrange = T, size = 1.3, color = "black") +
  ng1
plot

#Plot with points colored by Morph
plot <-
  ggplot(subset(visitshare.polldata, !is.na(Morph)), aes(x = Distance_to_core, y = Visits_per_Inf, colour = Morph)) +
  geom_point(alpha = 1) +
  labs(x = "Source population distance\nto urban centre (km)",
       y = "Number of pollinator visits\n per inflorescence") +
  geom_smooth(method = "lm", se = F, linetype = "longdash") +
  geom_smooth(method = "lm", se = F, linetype = "dashed")+
  ng1
plot

