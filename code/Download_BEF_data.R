## download data for Bees and Fish

# load libraries
library(RMariaDB)
library(tidyverse)
library(lubridate)

# files to access Winfree Lab SQL database (do not share!)
source("psw_wlab.R") #simply save a .R file as psw="PASSWORD"
source("user_wlab.R")

# helpers for data cleanups
lower <- function (x) {tolower(names(x))}
wlab_host<-"benedick.sebs.rutgers.edu"

conn <- dbConnect(RMariaDB::MariaDB(), user=user_wlab, password=psw_wlab,                  
                  dbname="bef_scale", port = 3306 , host= wlab_host)

dbListTables(conn)

vis_spec <- dbReadTable(conn, "nsf1718_spec")
# if we want to filter to only bees (but why?)
vis_spec <- vis_spec %>% filter(basicID=="bee")
# initial data cleanup
vis_spec$species[vis_spec$genus=="Lasioglossum" & vis_spec$species== "nr. obscurum" ] <- "obscurum"
#I think this is a fine choice; I examined one of the two specimens myself and buy it.  MR 20190218
vis_spec$species[vis_spec$genus=="Sphecodes" & vis_spec$species== "smilacinae?" ] <- "smilacinae"
# vis_spec$species[vis_spec$genus=="Ceratina" & vis_spec$species== "calcarata cf C. dupla" ] <- "calcarata_dupla_mikmaqi"
vis_spec$species[vis_spec$genus=="Lasioglossum" & vis_spec$species== "hitchensi/weemsi" ] <- "hitchensi_weemsi"

names(vis_spec)
vis_spec <- vis_spec %>%
  select(uniqueID
         , bee_genus = genus
         , bee_species = species
         , sex
         , fdate
         , year
         , plant_genus
         , plant_species
         , latitude
         , longitude
         , site
         )

species_deposition <- dbReadTable(conn, "single_visits")
names(species_deposition)

# summarize as per pollinator_id means (at species level?)

sampling_times <- dbReadTable(conn, "visit_rate_obs")
names(sampling_times)


# remember to disconnect 
dbDisconnect(conn)
# fish data from Lefcheck et al 2021
download.file("https://smithsonian.figshare.com/ndownloader/files/31149808"
              , destfile = "data/lefcheck.xlsx")
lefcheck<-readxl::read_excel("data/lefcheck.xlsx", sheet = 1)

write.csv(lefcheck, "data/lefcheck.csv", row.names =F)
file.remove("data/lefcheck.xlsx")

# not downloading Genung et al. 2022 data from dryad b/c those data didn't have years combined. 
# 
# # explore data 
# names(lefcheck)
# # is site the thing?
# lefcheck %>% 
#   group_by(SiteCode) %>% 
#   summarize(survs = n_distinct(SurveyID), depths = n_distinct(Depth)) %>% 
#   arrange(desc(survs))
# 
# # there are some sites with a ton of depths which makes me interested in that,
# # but also many sites less Could look at the biomass/depths correlation but for 
# # now maybe easier to just go for a higher group
# 
# lefcheck %>% summarize(ERs = n_distinct(Ecoregion)
#                         , Ps = n_distinct(Province)
#                         , Rs = n_distinct(Realm))
# # are these nested?
# lefcheck %>% group_by(Ecoregion) %>% summarize(Ps = n_distinct(Province)) %>% arrange(desc(Ps))
# 
# # look at sites per ER
# 
# lefcheck %>% 
#   group_by(Ecoregion) %>% 
#   summarize(sites = n_distinct(SiteCode)) %>% 
#   ggplot(aes(sites))+
#   geom_histogram()+
#   scale_x_log10()+
#   theme_classic()
# 
# 
# # messy, lots have <10
# 
# # try sites per Province?
# lefcheck %>% 
#   group_by(Province) %>% 
#   summarize(sites = n_distinct(SiteCode)) %>% 
#   ggplot(aes(sites))+
#   geom_histogram()+
#   scale_x_log10()+
#   theme_classic()
# 
# # Realm doesn't help much, better to keep replicates at Province level I think. 
# lefcheck %>% 
#   group_by(Realm) %>% 
#   summarize(sites = n_distinct(SiteCode)) %>% 
#   ggplot(aes(sites))+
#   geom_histogram()+
#   scale_x_log10()+
#   theme_classic()
# 
