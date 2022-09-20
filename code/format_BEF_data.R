# format datasets
library(tidyverse)


# Bee data are from same dataset as Genung et al. 2022 Rare&Declining:
# Visitation and estimated pollen deposition from Bees on 3 plant species at 
# 25 sites set in a grid in Central New Jersey Sampled in 2 years. We use only 
# 2017 data. Dryad repo https://doi.org/10.5061/dryad.qnk98sfkc

# Although each plant species was sampled in the same 25 sites, we treat the 
# samples from a given plant species as a distinct "metacommunity", as 
# pollination of *that plant* is the ecosystem function in this context. 



fix_bee <- function(x, value_name){pivot_longer(data.frame(x) 
                                                , cols = 3:last_col()
                                                , names_to = c("site", "year")
                                                , names_sep = "\\."
                                                , values_to = value_name) 
}


pnames <- c( "Polemonium_reptans","Phacelia_tanacetifolia", "Monarda_fistulosa")

bee_data <- readRDS("data/bee_data.rds")

# filter to just the sites of interest 
bee_abundances <- map_dfr(1:3, function(x){
  df <- data.frame(plantID = pnames[[x]], bee_data[[x]] %>% 
                     rownames_to_column("bee_gs"))
  return(df)
  })

bee_functions <- map_dfr(10:12, function(x){
  df <- data.frame(plantID = pnames[[x-9]], bee_data[[x]]%>% 
                     rownames_to_column("bee_gs"))
  return(df)
})


bee_abundances
bee_functions


bee_ab <- fix_bee(bee_abundances, "abund")
bee_ef <- fix_bee(bee_functions, "ef")

bee_complete <- left_join(bee_ab, bee_ef, by = c("plantID", "site", "year", "bee_gs"))

bee_complete %>% group_by(year, plantID) %>% summarize(sum(abund>0))
# 2017 looks a bit better

bee_2017 <- bee_complete %>% filter(year == 17)

# check out bee names

unique(bee_2017$bee_gs) %>% sort
# sort out some taxonomic issues (overall, pretty minor)

# make a couple species complexes

bee_2017$bee_gs[grepl("Hylaeus_aff", bee_2017$bee_gs)|grepl("Hylaeus_mod", bee_2017$bee_gs)] <- "Hylaeus_affinis_modestus"
bee_2017$bee_gs[grepl("hitch", bee_2017$bee_gs)|grepl("wee", bee_2017$bee_gs)] <- "Lasioglossum_hitchensi_weemsi"
# first, drop weird bombus 
bee_2017 <- bee_2017 %>% filter(!grepl("seeTN", bee_2017$bee_gs) & !grepl("Bombyl", bee_2017$bee_gs))
# need to group the complexes data now
bee_fixed <- bee_2017 %>% 
  mutate(study = "genung") %>% 
  group_by(syst = plantID, site, gen_sp = bee_gs, study) %>% 
  summarize(abund = sum(abund), ef = sum(ef)) %>% 
  filter(abund > 0) # keeping zeroes is unhelpful 

# head(bee_fixed)

# save the flat bee data to file 
write.csv(bee_fixed, "data/bees_for_analysis.csv", row.names = FALSE)


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

lefcheck_by_site %>% group_by(Realm) %>% summarize(n_distinct(Province))
lefcheck_by_site %>% 
  mutate(temperate = abs(SiteLat) > 23.44 ) %>% 
  group_by(temperate) %>% 
  summarize(n_distinct(Province))

reef_fish <- lefcheck_by_site %>% 
  ungroup() %>% 
  mutate(temperate = abs(SiteLat) > 23.44 ) %>% 
  group_by(Province) %>% 
  mutate(tempSite = ifelse(temperate, SiteCode, NA)) %>% 
  mutate(tempProv = (n_distinct(tempSite)/n_distinct(SiteCode))>0.5) %>% 
  mutate(syst = paste(Province, tempProv, sep = "_")
         , study = "lefcheck") %>% 
  select(study
         , syst
         , site = SiteCode
         , gen_sp = SPECIES_NAME
         , abund = Abundance
         , ef = Biomass)




# remove highly depauperate sites/systems with problematic sites
reef_fish <- reef_fish %>% 
  group_by(syst, site) %>% 
  mutate(smallSite = n()<5 & sum(abund)<30) %>% 
  ungroup() %>% 
  group_by(syst) %>% 
  mutate(small_syst = n_distinct(site)<10) %>% 
  filter(!smallSite & !small_syst) 
  
  
write.csv(reef_fish, "data/fish_for_analysis.csv")


# look at bci
load("data/bci.rdata")
head(bci.tree8)

# this is a bit slow here and unecessary 
# bci.tree8 %>% 
#   group_by(treeID) %>% 
#   summarize(agbs = n_distinct(agb), stems = n_distinct(stemID)) %>% 
#   arrange(desc(agbs))

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
    counter <- counter+1
    # "clunky but it works" way to select the right area
    bci.50[[counter]] <- bci.tree8[which(bci.tree8$gx >= x.seq[i] & bci.tree8$gx < x.seq[i+1] &
                                          bci.tree8$gy >= y.seq[j] & bci.tree8$gy < y.seq[j+1]),]
    bci.50[[counter]]$subplot <- paste0("HaSub_", counter)
  }
}

# checking the number of rows in each of the 50 new  1-ha plots, looks reasonable
# plot(unlist(lapply(bci.50, nrow)), ylim=c(0,12000))

# checking what the r-binded df looks like, looks good
bci_50subplots = do.call(rbind, bci.50)
# str(bci_50subplots)
# 
# # but need to make sure all the unique subplots showed up
# unique(bci_50subplots$subplot)

# looks good

bci_carbon <- bci_50subplots %>% 
  group_by(site = as.character(subplot), gen_sp = sp) %>% 
  summarize( abund = n(), ef = 0.5*sum(agb) ) %>% # carbon is 50% woody biomass
  mutate(syst = "BCI_sub_agc", study = "condit" )


bci_carbon
write.csv(bci_carbon, "data/bci_for_analysis",  row.names = FALSE)


# View(lefcheck %>% group_by(Province) %>% 
#   summarize(abundance = sum(Abundance)
#             , sites = n_distinct(SiteCode)
#   ))

# read in the tree data from TEAM
locs<- c("PA", "VB", "YA")
sixSiteForests<-map_dfr(locs, function(loc){
    dat <- read.csv(paste0("data/", loc, "_biomass.csv"))
    data.frame(loc, dat)})

# get areal extents

sixSiteForests %>% 
  filter(Year == 2012) %>% 
  group_by(loc, X1ha.Plot.Number) %>% 
  summarize(lon = mean(X1ha.Plot.X.Coordinate, na.rm= TRUE)
            , lat = mean(X1ha.Plot.Y.Coordinate, na.rm= TRUE)) %>% 
  drop_na() %>% # wd be good to have VB coords, currently missing
  group_by(loc) %>% 
  summarize(mdist = max(fossil::earth.dist(data.frame(.data$lon, .data$lat))))

TEAM_Forests<-sixSiteForests %>% 
  filter(Year == 2012) %>% 
  group_by(syst = loc
           , gen_sp = Gen_sp
           , site = X1ha.Plot.Number ) %>% 
  summarize(abund = n(), ef = sum(carbon)) %>% 
  mutate(study = "TEAM")


write.csv(TEAM_Forests, "data/team_forests_for_analysis.csv", row.names = FALSE)


bef_data <- bind_rows(
  TEAM_Forests
  , bci_carbon
  , reef_fish
  , bee_fixed
)

bef_data <- bef_data %>% 
  mutate(study_plus = if_else(study == "lefcheck"
                              , paste("reef_fish", c("temperate", "tropical")[1+as.numeric(grepl("TRUE", syst))], sep = "_") # TRUE
                              , if_else(study %in% c("condit", "TEAM") #FALSE
                                        , "tree_carbon" # TRUE
                                        , if_else(study == "genung" # fALSE
                                                  , "bee_pollination" # TRUE
                                                  , study))) )# FALSE




  
# summary(bef_data)
# head(bef_data)
# 
# unique(bef_data$study_plus)

write.csv(bef_data, "data/bef_data_for_analyses.csv", row.names = FALSE)

# # look at CV abundance
# ssf %>% group_by(loc) %>% summarize(sd(Abundance)/mean(Abundance))

# and pcf relationships
# ssf 
# 
# pdf("figures/pcf_in_TEAMS.pdf")
# ssf %>% 
#   mutate(siteCode = str_replace(.data$X1ha.Plot.Number, "VG-.*-", "")) %>% 
#   ggplot(aes(Abundance, totBiomass/Abundance, color = siteCode)) +
#   geom_smooth(se = FALSE)+
#   geom_point(alpha =0.5)+
#   facet_wrap(~loc)+
#   theme_classic() +
#   scale_x_log10()+
#   scale_y_log10()
# dev.off()
#                