# analyze data with Tina partition
library(tidyverse)
library(MeanRarity)
# loop over ell

ell_vec<-seq(-10, 10, 0.2)
res <- map_dfr(ell_vec, function(ell){
  lefcheck_by_site %>% summarize(ell = ell
                                 , D = rarity(Abundance, ell)
                                 , ab = sum(Abundance)
                                 , EFt = sum(Biomass)
                                 , pcf = EFt/ab)
} )

first_out<-map_dfr(unique(res$Province), function(Province){
    map_dfr(ell_vec, function(ell){
      dat = res %>% filter(Province == Province, ell==ell) 
  
      abMod = lm(log(dat$ab) ~ log(dat$D))
      pcfMod = lm(log(dat$pfc) ~ log(dat$D))
      ab.CI = confint(abMod, 2)
      pcf.CI = confint(pcfMod, 2)
      EF.CI = confint(lm(log(dat$EFt)~log(dat$D)), 2)
      return(data_frame(
        pcf.slope = coef(pcfMod)[[2]]
        , ab.slope = coef(abMod)[[2]]
        , EF.slope = pcf.slope + ab.slope
        , EF.cor = cor(log(dat$EFt)
                    , log(dat$D))
        , Ab.cor = ab.slope* sd(log(dat$D))/(sd(
          mean(log(dat$pcf)) * log(dat$ab)
        ))
        , pcf.cor = pcf.slope * sd(log(dat$D))/(sd(
          mean(log(dat$ab)) * log(dat$pcf)
        ))
      
      ))
    })
  })
