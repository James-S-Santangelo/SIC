## Sex and the City

### Authors: James S. Santangelo, Carole Advenard, L. Ruth Rivkin, and Ken A. Thompson

This repo contains code and data associated with a project testing the effects of urbanization on plant sexual and asexual reproduction in the model system white clover (_Trifolium repens_). The experiment has two questions:

1. Does urbanization alter relative investment in sexual vs. asexual reproduction?
2. Is urbanization associated with variation in plant phenological and reproductive traits?

#### Abstract

A growing body of evidence suggests that natural populations of dozens of species have undergone adaptive evolution in order to better tolerate the novel environmental conditions in urban areas. Invariably, studies of adaptive divergence in urban areas examine a single or few—often correlated—traits at a time from populations residing only at the most extreme urban and non-urban habitats, and do not control for environment-of-origin (e.g., maternal provisioning) effects. Thus, whether urbanization is driving divergence in many traits simultaneously in a manner that varies with the degree of local urbanization remains unclear. To examine whether urbanization drives clinal multivariate trait divergence, we generated seed families of white clover (Trifolium repens) collected from 27 populations along an urbanization gradient in Toronto, Canada, and used them to measure multiple phenotypic traits in a common garden. Our results show that families whose parents were collected from the most urban populations grew larger, had larger flower petals, experienced delayed germination and flowering, had thinner stolons, had reduced cyanogenesis, and were more attractive to pollinators. Each of these traits exhibited genetically-based changes that varied with the degree of urbanization of the source population. Field observations indicate that the pollinator community exhibits almost complete turnover between urban and nonurban sites, which potentially explains some of the observed divergence in reproductive traits. Our results suggest that urban populations are rapidly tuning their phenotypes to tolerate the local disturbances imposed by humans.


#### Using the code in this repository

To use the code in this repository to reproduce the manuscript's results, please follow the following steps:
1. `git clone` this repository or download it as a zip folder
2. Open `Rstudio`, go to `file > Open project` and open the `SIC.Rproj` Rproject associated with this repository
3. Once open, type `packrat::restore()` in the *R console* to download the packages required to perform the analyses.
