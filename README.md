# Hill-BEF

Hill diversities other than richness might predict ecosystem function better --
which scaling exponent is best?


# Contributing authors

Michael Roswell1 mroswell@umd.edu
Tina Harrison tinaharrison09@gmail.com
Mark A. Genung3,4 mark.genung@louisiana.edu

1. Dept. of Entomology, University of Maryland, College Park, MD 20742, USA
1. 
1. Dept of Biology, University of Louisiana, Lafayette, LA 70503, USA
1. Department of Ecology, Evolution & Natural Resources, Rutgers University, NJ
08901, USA


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
produces the strongest BEF relationship in several datasets, of real-world
ecosystem provision by wild, free-living organisms: pollination by wild bees,
above-ground carbon storage by tropical trees, and biomass provision from
reef-associated fish. We measured the strength of the correlation between
diversity and community-level as well as per-capita function. Unless richness is
always more predictive than other Hill diversities for both full-community and
per-capita function, adjusting the emphasis on rare versus common species may
reveal mechanisms underlying the BEF relationship, and provide more
readily-interpreted summaries of biodiversity change.

# Organization

Four directories. Anything used in the published manuscript will be at the top
level of one of these four directories, additional files from manuscript
preparation may be found in "archive" subdirectories.
- "code/" contains .R files for data wrangling, analysis, and producing figures. 
- "data/" contains files of various sorts, including source and derived
binaries, produced and manipulated along the workflow from downloading from a
repository to analyses
- "figures/" may be empty on git, but is the destination directory of any graph
produced by the code
- "Rmds" contains .rmd documents comprised of text and code developed for this
project. 

It should be possible to run all scripts in "code/" from the parent directory to
recreate all analyses and figures in the published MS. For archived files,
mileage will vary; we recommend reaching out to the authsors directly with any
questions about this content.

We expect that the files from the top-level of the "code" and "data" directories 
at a minimum will be archived according to journal guidelines by the time of 
publication. 

```
Hill-BEF
| README.md
| .gitignore
|
----- code
|     | Download_BEF_data.R # code that pulls datasets from online repositories. 
|     | format_BEF_data.R # format all datasets in similar way
|     | BEF_correlation_analysis.R # compute correlations and make graphs
|     | cartoon_fig.R # generate fake data to illustrate methods
|
----------- archive
|     |     | towards_general_partition.R # run Tina's partitions on any data (now lefcheck)
|     |     | simulate_correlations_basic.R # explore parameter space for underlying, abstract, non-biological relationships
|     |     | BEF_Hill_Forests_Feb22.R # explores diversity-ef relationships across ell-gradient and spatial scale for tropical forest datasets
|     |     | BEF_HIll_Forests2_Feb22.R # looks v similar to file without 2 in name. 
|     |     | Tina_paritition.R # minor changes from Tina's code
| 
----- data 
|     | lecheck.csv # key data from Lefcheck figshare
|     | bee_data.rds # data download from Genung et al. 2022 dryad
|     | bci.rdata # data download from Condit et al. 2019 dryad
|     | PA_biomass.csv # forest data Pasoh, Malaysia
|     | VB_biomass.csv # forest data Volcán Brava, Costa Rica
|     | YA_biomass.csv # forest data Yasuni, Ecuador
|     | fish_for_analysis.csv # formatted reef fish data
|     | bees_for_analysis.csv # formatted pollination data
|     | bci_for_analysis.csv # formatted BCI data
|     | team_forests_for_analysis.csv # formatted TEAM network forest data
|     | bef_data_for_analyses.csv # harmonized data for analyses
|
----- figures
|      | basic_simulation_heatmaps.pdf # output from simluation_correlations_basic.R
|
----- Rmds
|     | Hill_Postulae.Rmd # some general properties of Hill numbers

```

# Use and licensing: 
Condit et al. 2019 data are not subject to copyright and are available under a
CC 1.0 license from [https://doi.org/10.15146/5xcp-0d46](Dryad)

Lefcheck et al. 2021 data are used under a CC BY 4.0 license and are available
from
[https://figshare.com/articles/dataset/Species_richness_and_identity_both_determine_the_biomass_of_global_reef_fish_communities/16847029?file=31149808](figshare)

Genung et al. 2022 data are not subject to copyright and are available under a
CC 1.0 license from
[https://datadryad.org/stash/dataset/doi:10.5061/dryad.qnk98sfkc](Dryad)

The data from the now-defunct TEAM network ...

Other materials (e.g. code, figures) present here are freely available under a
CC BY 4.0 license. You can cite this dataset by ...






