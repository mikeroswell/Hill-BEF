# Hill-BEF

Hill diversities other than richness might predict ecosystem function better --
which scaling exponent is best?


# Contributing authors

Michael Roswell1 mroswell@umd.edu
Rachael Winfree2 rwinfree@rutgers.edu
Mark A. Genung3 mark.genung@louisiana.edu

1. Dept. of Entomology, University of Maryland, College Park, MD 20742, USA
1. Department of Ecology, Evolution & Natural Resources, Rutgers University, NJ
08901, USA
1. Dept of Biology, University of Louisiana, Lafayette, LA 70503, USA


# Synopsis 

Hill diversities have been promoted as robust metrics that integrate relative
abundances more transparently into biodiversity measurement. Hill diversities
summarize the distribution of relative abundances in an assemblage, and are
governed by a control parameter that determines the relative emphasis on rare
versus common species. Little empirical or theoretical work guides how to select
the value of the control parameter, however, and most biodiversity studies use
species richness, emphasizing the contribution of the rarest species to
diversity. As an important justification for measuring and conserving
biodiversity is protecting of the functions and services conferred by diverse
ecosystems, here we ask how the relationship between biodiversity and ecosystem
function (BEF) changes with the specification of species diversity through the
Hill diversity control parameter. We tested which value of the control parameter
produces the strongest BEF relationship in two large-scale, observational
datasets, one on pollination by wild bees and the second on carbon storage by
tropical trees. We measured the strength of the correlation between diversity
and community-level as well as per-capita function. Unless richness is always
more predictive than other Hill diversities for both full-community and
per-capita function, adjusting the emphasis on rare versus common species may
reveal mechanisms underlying the BEF relationship, and provide more
readily-interpreted summaries of biodiversity change.

# Organization
```
Hill-BEF
| README.md
| .gitignore
|
----- code
|     | BEF_Hill_mg.R  # main analysis file
|     | simulate_correlations_basic.R # explore parameter space for underlying, abstract, non-biological relationships
|     | BEF_Hill_Forests_Feb22.R # explores diversity-ef relationships across ell-gradient and spatial scale for tropical forest datasets
|     | BEF_HIll_Forests2_Feb22.R # looks v similar to file without 2 in name. 
| 
----- data
|     | a.list.rds # for now data in .rds file
|     | z.list.rds # for now data in .rds file
|     | BC_biomass.csv # tropical forest data from Barro Colorado Island
|     | PA_biomass.csv # forest data Pasoh, Malaysia
|     | VB_biomass.csv # forest data Volcán Brava, Costa Rica
|     | YA_biomass.csv # forest data Yasuni, Ecuador
|
----- figures
      | basic_simulation_heatmaps.pdf # output from simluation_correlations_basic.R
```



