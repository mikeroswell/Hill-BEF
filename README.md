# Biodiversity-ecosystem function relationships change in sign and magnitude across the Hill diversity spectrum


## Contributing authors

Michael Roswell1 mroswell@umd.edu
Tina Harrison2 tinaharrison09@gmail.com
Mark A. Genung2  mark.genung@louisiana.edu

1. Dept. of Entomology, University of Maryland, College Park, MD 20742, USA
1. Dept of Biology, University of Louisiana, Lafayette, LA 70503, USA


## Synopsis 

This repository contains data and code for our submission to the Phil Trans RSB
special issue on Biodiversity monitoring, "Biodiversity-ecosystem function
relationships change in sign and magnitude across the Hill diversity spectrum".
The repo also contains lots of notes and the starts of ideas that are not in the
MS. Upon article acceptance, the relevant components were migrated to 
[Zenodo](https://zenodo.org/doi/10.5281/zenodo.7846422)

## Keywords 
abundance, biodiversity, diversity profile, ecosystem function, Hill number,
rarity

## Summary 
Motivated by accelerating anthropogenic extinctions, decades of
biodiversity-ecosystem function (BEF) experiments show that ecosystem function
declines with species loss from local communities. Yet, at the local scale,
changes in species’ total and relative abundances are more common than species
loss. The consensus best biodiversity measures are Hill numbers, which use a
scaling parameter, ℓ, to emphasize rarer versus more common species. Shifting
that emphasis captures distinct, function-relevant biodiversity gradients beyond
species richness. Here, we surmised that Hill numbers that emphasize rare
species more than richness may distinguish large, complex, and presumably
higher-functioning assemblages from smaller and simpler ones. In this study, we
tested which values of ℓ produce the strongest BEF relationships in community
datasets of ecosystem functions provided by wild, free-living organisms. We
found that ℓ values that emphasized rare species more than richness most often
correlated most strongly with ecosystem functions. As emphasis shifted to more
common species, BEF correlations were often weak and/or negative. We argue that
unconventional Hill diversities that shift emphasis towards rarer species may be
useful for describing biodiversity change, and that employing a wide spectrum of
Hill numbers can clarify mechanisms underlying BEF relationships.


## Repo Organization

Four directories. Anything used in the published manuscript is contained in
"code", "data", or possibly "figures". "Rmds" contains notes and ideas not
included directly in the MS. Within the three directories with
publication-relevant content, the content for the published MS will be at the
top level of the directories, additional files from manuscript preparation may
be found in "archive" subdirectories.
- "code/" contains .R files for data wrangling, analysis, and producing figures. 
- "data/" contains files of various sorts, including source and derived
binaries, produced and manipulated along the workflow from downloading from a
repository to analyses
- "figures/" may be empty on git, but is the destination directory of any graph
produced by the code
- "Rmds" contains .rmd documents comprised of text and code developed for this
project. 

It should be possible to run all scripts in the top level of "code/" from the
parent directory to recreate all analyses and figures in the published MS. For
archived files, mileage will vary; we recommend reaching out to the authors
directly with any questions about this content.

## File organization

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
|      | # this directory contains many .pdf outputs that can be recreated with code.
|
----- Rmds
|     | Hill_Postulae.Rmd # some general properties of Hill numbers
|     | Jensens_inequality_with_weights.Rmd # sketch of Jensen's inequality
|     | empirical_ellBEF.Rmd # notes and snippets used in MS prep
|     | lefcheck_data_notes.Rmd # super early-stage project notes

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

The data from Cavanaugh et al. 2014 GEB was originally available from the
now-defunct TEAM network. We are republishing these data to make them available
again. We reached out to data holders to discuss this.

Other materials (e.g. code, figures) present here are freely available under a
CC BY 4.0 license. You can download and cite the version of record from Zenodo 
at https://zenodo.org/doi/10.5281/zenodo.7846422






