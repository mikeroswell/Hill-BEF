# Code to fit neg. binomials and reclaim the simulation parameters
library(MASS)
library(tidyverse)
library(furrr)
library(tictoc)
library(MeanRarity)

source("code/format_BEF_data.R")

# learn how to pull out neg. binomial params from a distribution
# ?fitdistr
# afit<-fitdistr(sample_infinite(fit_SAD()[[3]], 100)
#                , densfun = "negative binomial")
# 
# glm(family ="gamma")

library(ggplot2)
# function to fit neg binomial to abundances of species at the per-site level

nbpar <- function(ab){
  MASS::fitdistr(ab, densfun = "Negative Binomial"
               , lower=c(1e-9, 1e-9))}

trunc<-function(x){x[x>0]}

# simulate an abundance vector
set.seed(100)

trials<-667
size = 0.4
mu = 30

site_abundance <-rnbinom(n = trials, size = size, mu = mu)



# fit the distribution and get simulation parameters back out
nbpar(site_abundance) # returns something very close to simulated parameters

# fit again with zeros omitted
x<-nbpar(site_abundance[site_abundance>0]) # different parameters

drift <- data.frame()
for(driftSteps in c(1:40)){
  mypar <- nbpar(trunc(site_abundance))
  size <- mypar$estimate[[1]]
  mu <- mypar$estimate[[2]]
  site_abundance <- rnbinom(n = trials, size = size, mu = mu) 
  drift[driftSteps,"driftSteps"]<- driftSteps
  drift[driftSteps,"size"]<- size
  drift[driftSteps,"mu"]<- mu
  
}



drift %>% ggplot(aes(driftSteps, mu)) +
  geom_point() +
  theme_classic()

nbpar(trunc(rnbinom(667, size = 0.29, mu = 365)))
# function to fit gamma distributions to the nb parameters. Right now using a
# very klugey appraoch to making this a one-parameter problem, there is
# certainly a better one (uniroot?)
gpar <- function(nbparam){
  MASS::fitdistr(nbparam, densfun = "gamma"
                 , lower=c(1e-5, 0.999999999), upper = c(1e14, 1.0000000001) )}

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
  

  
  getNB<-map_dfr(unique(dat$SiteCode), function(site){
    print(site)
   nb <- dat %>% filter(SiteCode == site ) %>% pull(Abundance) %>% nbpar()
   return(data.frame(SiteCode = site
    , NBsize_est = nb$estimate[1]
     , NBsize_sd = nb$sd[1]
     , NBmu_est = nb$estimate[2]
     , NBmu_sd = nb$sd[2]))
   })
   getMeanGamma <- gpar(getNB$NBmu_est)$estimate[[1]]
   getSizeGamma <- gpar(getNB$NBsize_est)$estimate[[1]]
  
  lmDat<-dat %>% group_by(SiteCode) %>% summarize(Ab = sum(Abundance), EF = sum(Biomass))
  glmPar <- glm(lmDat$EF ~ log(lmDat$Ab), family = Gamma(link = "log")) 
  glmSum <- summary(glmPar)
  return(data.frame(Province = dat$Province
                    , meanSiteAb = sum(dat$Abundance)/length(dat$Abundance)
                    , gammaRichObs = length(unique(dat$SPECIES_NAME))
                    , n_sites = length(unique(dat$SiteCode))
                    , gammaShape_mu = getMeanGamma
                    , gammaShape_size = getSizeGamma
                    , EF_disp = glmSum$dispersion
                    , EF_slope = glmSum$coefficients[[2]]))
}

lefchecked <-  lefcheck %>% 
  group_by(Province) %>% 
  filter(n()>10) 

nc<- 7
plan(strategy = multiprocess, workers = nc)
province_params <- future_map_dfr(unique(lefchecked$Province), function(prov){
    dat <- lefchecked %>% 
      filter(Province == prov & Biomass>0 & !is.na(Abundance)) %>% 
     
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

