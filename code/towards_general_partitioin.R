# analyze data with Tina partition
library(tidyverse)
library(MeanRarity)
library(furrr)
# loop over ell

ell_vec<-seq(-10, 10, 0.2)

# set number of "cores" to use for parallelization
nc<-7

plan(strategy = "multiprocess", workers = nc)

res <- furrr::future_map_dfr(ell_vec, function(ell){
  lefcheck_by_site %>% summarize(ell = ell
                                 , D = rarity(Abundance, ell)
                                 , ab = sum(Abundance)
                                 , EFt = sum(Biomass)
                                 , pcf = EFt/ab)
} )

# # check my work
# res %>% ggplot(aes(ell, D, color = Province))+
#   geom_point(alpha = 0.2)+
#   theme_classic()


first_out <- map_dfr(unique(res$Province), function(Province){
    future_map_dfr(ell_vec, function(ell){
      dat = res %>% filter(Province == !!Province, ell == !!ell) 
  # was struggling with the quoting stuff, the !! fixes
      abMod = lm(Ln(dat$ab) ~ Ln(dat$D))
      pcfMod = lm(Ln(dat$pcf) ~ Ln(dat$D))
      ab.CI = confint(abMod, 2)
      pcf.CI = confint(pcfMod, 2)
      EF.CI = confint(lm(Ln(dat$EFt)~Ln(dat$D)), 2)
    
      pcf.slope = coef(pcfMod)[[2]]
      ab.slope = coef(abMod)[[2]]
      EF.slope = pcf.slope + ab.slope
      EF.cor = cor(Ln(dat$EFt)
                  , Ln(dat$D))
      ab.cor = ab.slope* sd(Ln(dat$D))/(sd(
        mean(Ln(dat$pcf)) * Ln(dat$ab)
      ))
      pcf.cor = pcf.slope * sd(Ln(dat$D))/(sd(
        mean(Ln(dat$ab)) * Ln(dat$pcf)
      ))
      return(data.frame(pcf.slope
                        , ab.slope
                        , EF.slope
                        , EF.cor
                        , ab.cor
                        , pcf.cor
                        , ell
                        , Province
                        , dataset = "lefcheck"
      ))
    })
  })

pdf("figures/lefcheck_partitions.pdf")
map(unique(first_out$Province), function(province){
  first_out %>% 
    filter(Province == !!province) %>%
    pivot_longer(cols = ends_with(".slope")
                 , values_to = "slope"
                 , names_to ="component") %>% 
    ggplot(aes(ell, slope, color = component))+
    geom_line()+
    theme_classic()+
    geom_hline(yintercept = 0, size = 0.2)  +
    geom_vline(xintercept = 1, size = 0.2) +
    scale_color_manual(values = c("red", "black", "blue")) +
    ggtitle(province)
    
})

dev.off()
# for now save data so I can restart computer
write.csv(first_out, "data/lefcheck_results.csv", row.names = F)
first_out<-read.csv("data/lefcheck_results.csv")
