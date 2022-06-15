# Code to fit neg. binomials and reclaim the simulation parameters
library(MASS)
library(tidyverse)
library(furrr)
library(tictoc)
library(MeanRarity)

source("code/format_BEF_data.R")

# learn how to pull out neg. binomial params from a distribution
?fitdistr
afit<-fitdistr(sample_infinite(fit_SAD()[[3]], 100)
               , densfun = "negative binomial")

glm(family ="gamma")

# at least for now, this function only takes a df as its argument. 
# this means that columns have to be named according to the lefcheck data
# `Province` names a collection of sites
# `Abundance` is abundance, `Biomass` is species-level EF
# `SiteCode` names an individual site and `SPECIES_NAME` an individual sp. 
par_getter <- function(dat){
  if(length(dat$Abundance) < 10){return(data.frame(Province = dat$Province
                                                   , meanSiteAb = sum(dat$Abundance)/length(dat$Abundance)
                                                   , gammaRichObs = length(unique(dat$SPECIES_NAME))
                                                   , n_sites = length(unique(dat$SiteCode))
                                                   , size_est = NA
                                                   , size_sd = NA
                                                   , mu_est = NA
                                                   , mu_sd = NA
                                                   , EF_disp = NA
                                                   , EF_slope = NA))
  }
  
  nbpar <- MASS::fitdistr(dat$Abundance, densfun = "Negative Binomial"
                          , lower=c(1e-5, 1e-5))
  lmDat<-dat %>% group_by(SiteCode) %>% summarize(Ab = sum(Abundance), EF = sum(Biomass))
  glmPar <- glm(lmDat$EF ~ log(lmDat$Ab), family = Gamma(link = "log")) 
  glmSum <- summary(glmPar)
  return(data.frame(Province = dat$Province
                    , meanSiteAb = sum(dat$Abundance)/length(dat$Abundance)
                    , gammaRichObs = length(unique(dat$SPECIES_NAME))
                    , n_sites = length(unique(dat$SiteCode))
                    , size_est = nbpar$estimate[1]
                    , size_sd = nbpar$sd[1]
                    , mu_est = nbpar$estimate[2]
                    , mu_sd = nbpar$sd[2]
                    , EF_disp = glmSum$dispersion
                    , EF_slope = glmSum$coefficients[[2]]))
  }

lefchecked <-  lefcheck %>% 
  group_by(Province) %>% 
  filter(n()>10) 
province_params <- map_dfr(unique(lefchecked$Province), function(prov){
  dat <- lefchecked %>% 
    filter(Province == prov & Biomass>0) %>% 
   
    ungroup() %>% 
    drop_na()
  print(length(dat$Abundance))
  par_getter(dat)
})
 

pdf("figures/simulation_params_from_empirical_data.pdf")
pairs(province_params[,c(2,3,4,5, 7, 9, 10 )]
      , gap = 1/6
      , upper.panel = NULL)
dev.off()

