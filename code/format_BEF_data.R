# format datasets
library(tidyverse)


# Bee data are from same dataset as Genung et al. 2022 Rare&Declining:
# Visitation and estimated pollen deposition from Bees on 3 plant species at 
# 25 sites set in a grid in Central New Jersey Sampled in 2 years. However,
# visits and samples from the two years have been combined in this analysis and,
# because "abundance" is measured as a per-minute visitation rate and sample
# times weren't always identical across sites and years, the numbers are a bit
# different from those in the Dryad repo https://doi.org/10.5061/dryad.qnk98sfkc

# Although each plant species was sampled in the same 25 sites, we treat the 
# samples from a given plant species as a distinct "metacommunity", as 
# pollination of *that plant* is the ecosystem function in this context. 

# MG provided data source files. Raw data housed in Winfree Lab database at
# Rutgers U.

fix_bee <- function(x, value_name){pivot_longer(data.frame(x) %>% 
                                                  rename(bee_species = X)
                                                , cols = 2:last_col()
                                                , names_to = "site"
                                                , values_to = value_name) 
}

beef <- list.files(path = "data", pattern = ".*_Apr22\\.csv")[1]
pnames <- c( "Polemonium_reptans","Phacelia_tanacetifolia", "Monarda_fistulosa")

bee_dup <- map_dfr(list.files(path = "data", pattern = ".*_Apr22\\.csv"), function(beef){
  dat <- read.csv(paste0("data/",beef))
  value_name <- gsub( "_.*", "", beef)
  plant_code <- gsub("(^.*_output_)(.*)(_Ap.*$)", "\\2", beef) 
  plant_ID <- pnames[gsub("(^[a-z]{1})([a-z]*_)([a-z]{1})(.*$)"
                          , "\\1\\3"
                          , tolower(pnames))==plant_code]
  return(data.frame(fix_bee(dat, value_name), plant_ID))
})


bee_dat <- left_join(bee_dup %>% filter(is.na(func)) %>% select(-func)
                     , bee_dup %>% filter(is.na(abund)) %>% select(-abund)
                     , by = c("bee_species", "site", "plant_ID") ) %>% 
  select(plant_ID, site, bee_species, abund, ef = func) %>% 
  arrange(plant_ID, site, bee_species)

write.csv(bee_dat, "data/bef-scale_bees_taxonomy_issues.csv", row.names = FALSE)

unique(bee_dat$bee_species) %>% sort

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
               