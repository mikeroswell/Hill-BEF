# format datasets
library(tidyverse)
# the bee data consist of a bunch of lists


# format the bee data
a.list = readRDS("data/a.list.rds")
z.list = readRDS("data/z.list.rds")
a.data <- a.list %>% lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% as.matrix}) 
z.t.data <- z.list %>% lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% as.matrix}) 
z.data = sapply(1:length(z.t.data),
                function(se) {z = z.t.data[[se]]/a.data[[se]]
                z[is.nan(z)] = 0; z})

# Might go a bit flatter for the lefcheck data
lefcheck<-read.csv("data/lefcheck.csv")

lefcheck_by_site<-lefcheck %>% group_by(Province
                      , Ecoregion
                      , Realm
                      , SiteCode
                      , SiteLat
                      , SiteLong
                      , SPECIES_NAME) %>% 
  summarize(Biomass = sum(Biomass)
            , Abundance = sum(Abundance)
            , depths = n_distinct(Depth)
            , survs = n_distinct(SurveyID))