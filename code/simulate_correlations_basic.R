# script by michael roswell to simulate parameter space
library(tidyverse)
library(MeanRarity)
library(tictoc)
library(furrr)
# function to simulate correlated vector 
# from https://stackoverflow.com/a/52880462/8400969

simcor <- function (x, ymean=0, ysd=1, correlation=0) {
  n <- length(x)
  
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * 
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}



n_sp <- 350 # number of possible species
n_comm <- 20 # number of communities
ab_mean <- 35 # mean species abundance per community
ab_sd <- 2 # sd of mean species abundance across communities
ab_disp <- 0.2 # dispersion of negative binomial
cor_vec<-seq(-0.8,0.8,0.1) # correlation between species abundance and per-capita ecosystem function
ef_mean<-25 # constant, mean per-capita function
ef_sd<- 5 # conditional SD of per-capita function
ell_vec<-seq(-5,15,0.25)

to_rare<-function(x){
  sapply(x, function(y){
    ifelse(y>0, sum(x)/y, y)
  })
  }

Ln = function(x) {
  ifelse((x != 0 & !is.na(x)), log(x), 0)}
Exp = function(x, y) {ifelse(x != 0, x^y, 0)}

nonzero <-function(x, y){x[(x*y)!=0]}
#ignore non-functionals (absent or 0 function)
cor_rm<-function(x, y){
  cor(nonzero(x,y), nonzero(y,x))
}
# test it 
# x <- c(0, 0, 0, 5, 6, 7, 8)
# y <- c(0, 1, 0, 7, 3, 5, 4)
# cor(x,y)
# cor_rm(x,y)
# cor_rm(Ln(x), Ln(y))


reps<-99

# testing values for simulation
iter<-1
rare_ef_cor<- 0.6
ell<-1


nc <- 7 # cores for parallelization
plan(strategy = multiprocess, workers = nc)

tic()
BEF_sim_simple<-future_map_dfr(1:reps, function(iter){
  map_dfr(cor_vec, function(rare_ef_cor){
    # abundance: species in rows, communities in columns
    # ab = sapply(1:n_comm, function(x){rnbinom(n_sp, size = ab_disp, mu = ab_mean)})
    # in v2, add variation in community size
    ab = sapply(1:n_comm, function(x){rnbinom(n_sp
                                              , size = ab_disp
                                              , mu = rnorm(1
                                                           , mean = ab_mean
                                                           , sd = ab_sd))})
    rarities = apply(ab, 2, to_rare)
    # for all these communities in this iteration, set per-capita function based on species abundance
    ef = apply(ab, 2, function(x){simcor(to_rare(x)
                                         , correlation = rare_ef_cor
                                         , ymean = ef_mean
                                         , ysd = ef_sd)})
    # summarize ef
    species_ef = ab*ef
    per_capita_ef = apply(species_ef, 2, sum)/ apply(ab, 2, sum)
    total_ef = apply(species_ef, 2, sum)
    map_dfr(ell_vec, function(ell){
      D = apply(ab, 2, function(x){rarity(x, l = ell)})
      per_capita_BEF = cor_rm(sapply(D, Ln), sapply(per_capita_ef, Ln))^2
      per_community_BEF = cor_rm(sapply(D, Ln), sapply(total_ef, Ln))^2
      return(data.frame(rare_ef_cor
                        , ell
                        , v_rich = var(apply(ab, 2, function(x){sum(x>0)}))
                        , per_capita_BEF
                        , per_community_BEF))
    })
  })
})
toc() # 5.5 minutes on MR macbook pro
# check variation in richness as driver

BEF_summary<-BEF_sim_simple %>% 
  group_by(rare_ef_cor, ell) %>% 
  summarize_all(.funs = "mean")

 
# BEF_summary %>% ggplot(aes(v_rich, per_community_BEF))+
#   geom_point() + # will funnel out if this is the thing (variation in bef increasing with variation in richness)
#   theme_classic() +
#   labs(x = ("variance in richness")) 

pdf("figures/basic_simulation_heatmaps.pdf")
# make heatmap for per-capita function
BEF_summary %>% ggplot(aes(rare_ef_cor, ell, fill = per_capita_BEF))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis_c() +
  labs(x = "cor(rarity, per-capita ef)"
         , y = "Hill scaling exponent ell"
         , fill = "BEF (per-capita) r2 \nin log-log space")

# make heatmap for total function

BEF_summary %>% ggplot(aes(rare_ef_cor, ell, fill = per_community_BEF))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis_c() +
  labs(x = "cor(rarity, per-capita ef)"
       , y = "Hill scaling exponent ell"
       , fill = "BEF (per-community) r2 \nin log-log space")

dev.off()


