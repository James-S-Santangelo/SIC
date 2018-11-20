# ken a thompson analysis

# load packages
library(tidyverse)

polldata <- read.csv('data-clean/pollinatorObservations.csv')

# analysis of pollinator observations
# where are the plants getting pollinated more?

popmeans.polldata <- polldata %>% 