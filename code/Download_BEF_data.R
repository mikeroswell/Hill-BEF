## download data

# load libraries
library(tidyverse)

# bee data from Genung et al. 2022
download.file("https://datadryad.org/stash/downloads/file_stream/1678058"
              , destfile = "data/bee_data.rds")


# bee_data <- readRDS("data/bee_data.rds")
# head(bee_data)
# will clean up in format_data.R

# fish data from Lefcheck et al 2021
download.file("https://smithsonian.figshare.com/ndownloader/files/31149808"
              , destfile = "data/lefcheck.xlsx")
lefcheck<-readxl::read_excel("data/lefcheck.xlsx", sheet = 1)

write.csv(lefcheck, "data/lefcheck.csv", row.names =F)
file.remove("data/lefcheck.xlsx")

# BCI tree data from Condit et al. 2019

# get the tree zip file
download.file("https://datadryad.org/stash/downloads/file_stream/148942"
              , destfile = "data/BCI_from_Dryad.zip")

#unzip locally 
unzip("data/BCI_from_Dryad.zip", exdir = "data")
# pull just the file we need
file.copy("data/home/fullplotdata/tree/bci.tree8.rdata", "data/bci.rdata")
# remove other stuff
unlink("data/home", recursive = TRUE)
file.remove("data/BCI_from_Dryad.zip")
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
