## download data from the internet
library(tidyverse)

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
