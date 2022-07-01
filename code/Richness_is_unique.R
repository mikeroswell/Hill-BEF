# what is going on with richness and why is it different?
# one feature of richness is that it is more discrete than other diversities
# it is also the flip-point between emphasizing species progressively more the 
# rarer they are vs. the common they are. 
#
# The goal of this script is to drop a lot of the complexity of previous 
# simulations and see what there is to see about richness. 

library(MeanRarity)
library(tidyverse)
library(patchwork)

seeds <- 1:20

pdf("figures/dumb_simulation_unevenness.pdf")
map(seeds, function(seed){
  set.seed(seed)
  # set number of sites
  nSites <- 20
  # make a large variation in richness
  richVec <- rnbinom(nSites, size = 3, mu = 90)
  # hist(richVec)
  # set evennesses
  eVec <- runif(nSites, min = 0.02, max = 0.3)
  # set site abundances (independent from richness for now)
  # abVec <- rnbinom(nSites, size = 2, mu = 3000)
  abVec <- rnbinom(nSites, size = 2, mu = 1000 + 15 * richVec )
  # also (risky) set site-level EF as a noisy function of abundance
  # EFVec <- rlnorm(nSites, meanlog = 2 + 0.5 *log(abVec), sdlog = 0.7) 
  EFVec <- rnorm(nSites, mean = 5+20/(2+exp(-1/700*(abVec-2500))), sd = 3)
  
  ellVec <- seq(-5, 5, 0.2)
  # generate SADs (do lognormal for now)
  BEFdat <- map_dfr(1:nSites, function(site){
    rich <- richVec[site]
    # print(rich)
    evenness <- eVec[site]
    # print(evenness)
    rel_ab <- fit_SAD(rich = rich
                  , simpson = evenness*(rich-1)+1
                  , distr = "lnorm")[[3]]
    ab <- sample_infinite(rel_ab, size = abVec[site] )
    print(paste("rich_sim =", rich, "\nrich obs =", sum(ab>0)))
    map_dfr(ellVec, function(ell){
      D <- rarity(ab, l = ell)
      return(data.frame(site
                        , evenness
                        , totAb = abVec[site]
                        , ell
                        , D
                        , EF = EFVec[site]))
    })
  })
  
  BEF_cors <- BEFdat %>% group_by(ell) %>% summarize(BEF = cor(log(D), log(EF))
                                                     , ell_ab = cor(log(D)
                                                                    , log(totAb)
                                                                    ))
  
  ab_EF <- BEFdat %>% ggplot(aes(totAb, EF)) +
    geom_point() +
    theme_classic()
  
  ab_ell <- BEF_cors %>% ggplot(aes(ell, ell_ab)) + 
    geom_point() + 
    ylim(c(-1,1)) + 
    labs(title = seed, y = "correlation D vs. totAb (log-log)") +
    theme_classic() 
  
  
  
  BEF <- BEF_cors %>% ggplot(aes(ell, BEF)) + 
    geom_point() + 
    ylim(c(-1,1)) + 
    labs(title = seed) +
    theme_classic() 
  
  (ab_EF + ab_ell)/BEF
  
})

dev.off()
