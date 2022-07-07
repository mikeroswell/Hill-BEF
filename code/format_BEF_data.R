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

lefcheck_by_site<-lefcheck %>% 
  group_by(Province) %>% 
  mutate(Province_sites = n_distinct(SiteCode)) %>% 
  filter(Province_sites >10) %>% 
  select(-Province_sites) %>% 
  group_by(Province
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
# look at bci
load("data/bci_tree8.rdata")
head(bci.tree8)

bci.tree8 %>% 
  group_by(treeID) %>% 
  summarize(agbs = n_distinct(agb), stems = n_distinct(stemID)) %>% 
  arrange(desc(agbs))

bci_sum <- bci.tree8 %>% group_by(sp) %>% 
  summarize(Biomass = sum(agb)
            , Abundance = n()
            , quadrats =n_distinct(quadrat))





# View(lefcheck %>% group_by(Province) %>% 
#   summarize(abundance = sum(Abundance)
#             , sites = n_distinct(SiteCode)
#   ))

# read in the tree data from TEAM
locs<- c("PA", "VB", "YA")
sixSiteForests<-map_dfr(locs, function(loc){
    dat <- read.csv(paste0("data/", loc, "_biomass.csv"))
    data.frame(loc, dat)})

ssf<-sixSiteForests %>% group_by(loc, Gen_sp, Year, X1ha.Plot.Number ) %>% 
  summarize(Abundance = n(), totBiomass = sum(biomass), totCarbon = sum(carbon))

ssf
