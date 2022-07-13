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



range(bci.tree8$gx, na.rm=T) # 0 to 1000 is the range of x (10 ha)
range(bci.tree8$gy, na.rm=T) # 0 to 500 is the range of y (5 ha)

row.divs <- 11
col.divs <- 6

x.seq = seq(from=0, to=1000, length.out =row.divs) # 10 divisions on x, length 11 because we need endpoint
y.seq = seq(from=0, to= 500, length.out = col.divs) #  5 divisions on y, length  6 because we need endpoint

counter = 0 # this will help us combine unique i and j into one variable for the output list
bci.50 = vector('list', 50) # empty list, one for each of the 1-ha plots we want to make

for(i in 1:(row.divs-1)){ # x axis
  for(j in 1:(col.divs-1)){ # y axis
    counter = counter+1
    # "clunky but it works" way to select the right area
    bci.50[[counter]] = bci.tree8[which(bci.tree8$gx >= x.seq[i] & bci.tree8$gx < x.seq[i+1] &
                                          bci.tree8$gy >= y.seq[j] & bci.tree8$gy < y.seq[j+1]),]
    bci.50[[counter]]$subplot = counter
  }
}

# checking the number of rows in each of the 50 new  1-ha plots, looks reasonable
plot(unlist(lapply(bci.50, nrow)), ylim=c(0,12000))

# checking what the r-binded df looks like, looks good
bci_50subplots = do.call(rbind, bci.50)
str(bci_50subplots)

# but need to make sure all the unique subplots showed up
unique(bci_50subplots$subplot)

# looks good






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
