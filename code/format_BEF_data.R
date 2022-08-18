# format datasets
library(tidyverse)


# Bee data from Genung et al. 2022 Realistic Loss Order:
# Visitation and estimated pollen deposition from Bees on 3 plant species at 
# 25 sites set in a grid in Central New Jersey Sampled in 3 years
# Although each plant species was sampled in the same 25 sites, we treat the 
# samples from a given plant species as a distinct "metacommunity", as 
# pollination of *that plant* is the ecosystem function in this context. 
# For source data we used .rds files provided by M Genung. 


# format the bee data
a.list = readRDS("data/a.list.rds")
z.list = readRDS("data/z.list.rds")
# name the plants
#################################################
### NEED TO CHECK THE ORDER, THIS IS MADE UP!  ##
##################################################
pnames <- c( "Polemonium_reptans","Phacelia_tanacetifolia", "Monarda_fistulosa")
names(a.list) <- pnames
names(z.list) <- pnames

a.list
a.data <- a.list %>% 
  lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% 
    as.matrix}) 
z.t.data <- z.list %>% lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% as.matrix}) 
z.data = sapply(1:length(z.t.data),
                function(se) {z = z.t.data[[se]]/a.data[[se]]
                z[is.nan(z)] = 0; z})

z.data
a.data
names(z.data) <- pnames

fix_bee <- function(x, value_name){pivot_longer(data.frame(x) %>% 
                                      rownames_to_column("bee_species")
                                    , cols = 2:last_col()
                                    , names_to = c("site", "year")
                                    , names_sep = "\\."
                                    , values_to = value_name) 
  }
bee_abund <- purrr::map_dfr(a.data, function(x){fix_bee(x, value_name = "visit_rate")}, .id = "plant")
bee_pcf <- purrr::map_dfr(z.data, function(x){fix_bee(x, value_name = "pollen_dep")}, .id = "plant")

bees_1819 <- left_join(bee_abund, bee_pcf) %>% 
  group_by(plant, bee_species, site) %>% 
  summarize(visit_rate = mean(visit_rate), pollen_dep = mean(pollen_dep)) %>% 
  ungroup()

bees_1819
# Fish data from Lefcheck et al. 2021 Nature Communications "Species Richness
# and Identity..."
# downloaded with script "code/Download_BEF_data.R"
# source("code/Download_BEF_data.R")

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
  summarize(Abundance = n(), totBiomass = sum(biomass), totCarbon = sum(carbon)) %>% filter(Year == 2012)

# look at CV abundance
ssf %>% group_by(loc) %>% summarize(sd(Abundance)/mean(Abundance))

# and pcf relationships
# ssf 

pdf("figures/pcf_in_TEAMS.pdf")
ssf %>% 
  mutate(siteCode = str_replace(.data$X1ha.Plot.Number, "VG-.*-", "")) %>% 
  ggplot(aes(Abundance, totBiomass/Abundance, color = siteCode)) +
  geom_smooth(se = FALSE)+
  geom_point(alpha =0.5)+
  facet_wrap(~loc)+
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()
dev.off()
               