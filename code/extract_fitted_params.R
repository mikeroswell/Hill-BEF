# Code to fit neg. binomials and reclaim the simulation parameters
library(MASS)
library(tidyverse)
library(furrr)
library(tictoc)
library(MeanRarity)

# learn how to pull out neg. binomial params from a distribution
?fitdistr
afit<-fitdistr(sample_infinite(fit_SAD()[[3]], 100)
               , densfun = "negative binomial")

glm(family ="gamma")

province_params <- map_dfr(unique(lefcheck$Province), function(prov){
  dat <- lefcheck %>% 
    filter(Province == prov & Biomass>0) %>% 
    drop_na()
  print(length(dat$Abundance))
  if(length(dat$Abundance) < 10){return(data.frame(Province = prov
                                                  , meanSiteAb = sum(dat$Abundance)/length(dat$Abundance)
                                                  , obsRich = length(unique(dat$SPECIES_NAME))
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
  return(data.frame(Province = prov
                    , meanSiteAb = sum(dat$Abundance)/length(dat$Abundance)
                    , obsRich = length(unique(dat$SPECIES_NAME))
                    , size_est = nbpar$estimate[1]
                    , size_sd = nbpar$sd[1]
                    , mu_est = nbpar$estimate[2]
                    , mu_sd = nbpar$sd[2]
                    , EF_disp = glmSum$dispersion
                    , EF_slope = glmSum$coefficients[[2]]))
})

pdf("figures/simulation_params_from_empirical_data.pdf")
pairs(province_params[,c(2,3,4, 6, 8, 9 )]
      , gap = 1/1000
      , upper.panel = NULL)
dev.off()

dprovince_params %>% pivot_longer(totAb, size_est, size_sd, mu_est, mu_sd, EF_disp, EF_slope
                                 , names_to = "param"
                                 , values_to = "val"
                                 ) %>% 
  ggplot(aes())
  


