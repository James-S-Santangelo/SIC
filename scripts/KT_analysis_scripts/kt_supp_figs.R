# author: ken a thompson 
# title: clover multivariate phenomic cline supp mat figs
# date: 2018-12-11

# load packages

library(cowplot)
library(tidyverse)

# load data

familyMeans <- read_csv("data-clean/experimentalData_familyMeans.csv")

# ggplot theme
#Theme used for plots throughout script
ng1=theme(aspect.ratio=0.7,panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1),
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black",size=10),
          axis.title=element_text(color="black",size=1),
          axis.title.y=element_text(vjust=2,size=17),
          axis.title.x=element_text(vjust=0.1,size=10),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          strip.text.x = element_text(size = 7, colour = "black",face = "bold"),
          strip.background = element_rect(colour="black"),
          legend.position = "top", legend.direction="vertical",
          legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))


traitchange.function <- function(data, trait, trait_name){
  ggplot(data, aes_string(x = "Distance_to_core", y = trait)) +
    geom_point(size = 2, position = position_dodge(width = 0.1)) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", size = 2) +
    xlab("distance to core (km)") + ylab(trait_name) +
    ng1
}

# make the univariate plots
## hcn
UniVar.HCN <- traitchange.function(familyMeans, "freqHCN", "frequency of cyanogenesis")

## days to first flower & germination
UniVar.FirstFlower <- traitchange.function(familyMeans, "Days_to_flower", "# days to first flower")
UniVar.Germtime <- traitchange.function(familyMeans, "Time_to_germination", "# days to germination")

## banner petal size
UniVar.BnrL <- traitchange.function(familyMeans, "Avg_bnr_lgth", "banner petal length (mm)")
UniVar.BnrW <- traitchange.function(familyMeans, "Avg_bnr_wdth", "banner petal width (mm)")

## vegetative biomass
UniVar.VegBM <- traitchange.function(familyMeans, "Veget_biomass", "vegetative biomass (g)")

## stolon thickness
UniVar.Stolon <- traitchange.function(familyMeans, "Avg_stolon_thick", "stolon thickness (mm)")

Fig_S_X.UniVar <- plot_grid(UniVar.FirstFlower, UniVar.Germtime, UniVar.BnrL, UniVar.BnrW, UniVar.VegBM, UniVar.Stolon, UniVar.HCN, ncol = 3, labels = "AUTO")

