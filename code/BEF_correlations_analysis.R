# analyses of correlations between log diversity, log function, and also some 
# relationships between those two variables and abundance and per-capita 
# function.

# load packages. 
# devtools::install_github("mikeroswell/MeanRarity")
library(MeanRarity) 
library(tidyverse)
library(furrr)
library(tictoc)

# helper functions
# cubeRt <- function(x){sign(x)*abs(x)^(1/3)}

# load data
bef_data <- read.csv("data/bef_data_for_analyses.csv")
bef_data <- bef_data %>% filter(abund >0)

# vector of Hill diversity scaling factors to consider
ell_vec<-seq(-10, 10, 0.05)



# tempting to look at greater density of points near ell = 0
# ell_vec <- c(exp(seq(-12,3, 0.1)), -(exp(seq(-12,3, 0.1))))
# ell_vec <- cubeRt(c(seq(-80, 80, 0.5)
#                     , seq(-0.5, 0.5, 0.01)))

# set number of "cores" to use for parallelization
# on intel mac at least, hyperthreading is used so the maximum is double the 
# number of physical cores
nc <- 7



# this could use some work to become more general but works fine for lefcheck
# remember that this works on groups only if they exist (e.g. site)
# this was actually slower
# sum_by_ell <- function(dat
#                        , ell_vec
#                        ,  ...){
#   dat2 <- dat %>% mutate(ab = sum(abund)
#                          , EFt = sum(ef)
#                          , pcf = EFt/ab)
#   res <- furrr::future_map_dfr(ell_vec, function(ell){
#     dat2 %>% 
#       group_by(ab, EFt, pcf, .add = TRUE) %>% 
#       summarize(ell = ell
#                 , D = rarity(abund, ell)
#                 
#                 # , meanlat = mean(abs(SiteLat))
#       )})
#   return(res)
# }

# housekeeping: make a table for coloring all graphs
set.seed(42)
clr_fac <- sample(letters)
clr_tab <- bef_data %>% 
  group_by(study_plus, syst) %>%
  summarize(Nsites = n_distinct(site)) %>% 
  mutate(clr = clr_fac[row_number(syst)])



sum_by_ell <- function(dat
                       , ell_vec
                       ,  ...){
  furrr::future_map_dfr(ell_vec, function(ell){
    dat %>% summarize(ell = ell
                       , D = rarity(abund, ell)
                       , ab = sum(abund)
                       , EFt = sum(ef)
                       , pcf = EFt/ab
                     )})
}


# times the crunching, if you want
tictoc::tic()
future::plan(strategy = "multiprocess", workers = nc)
bef_by_ell <- sum_by_ell(bef_data %>% group_by(site, syst, study, study_plus), ell_vec = ell_vec )

tictoc::toc()
# under a minute MBP

bef_by_ell

fit_lms <- function(sub){
  furrr::future_map_dfr(ell_vec, function(ell){
    dat = sub %>% filter(ell == !!ell) 
    # was struggling with the quoting stuff, the !! fixes
    abMod = lm(Ln(dat$ab) ~ Ln(dat$D))
    pcfMod = lm(Ln(dat$pcf) ~ Ln(dat$D))
    # ab.CI = confint(abMod, 2)
    # pcf.CI = confint(pcfMod, 2)
    # EF.CI = confint(lm(Ln(dat$EFt)~Ln(dat$D)), 2)
    
    pcf.slope = coef(pcfMod)[[2]]
    ab.slope = coef(abMod)[[2]]
    EF.slope = pcf.slope + ab.slope
    EF.cor = cor(Ln(dat$EFt)
                 , Ln(dat$D))
    ab.cor = cor(Ln(dat$ab), Ln(dat$D))
    pcf.cor = cor(Ln(dat$D), Ln(dat$pcf))
    EF.spe = cor(Ln(dat$EFt)
                 , Ln(dat$D)
                 , method = "spearman")
    ab.spe = cor(Ln(dat$ab), Ln(dat$D), method = "spearman")
    pcf.spe = cor(Ln(dat$D), Ln(dat$pcf), method = "spearman")
    return(data.frame(pcf.slope
                      , ab.slope
                      , EF.slope
                      , EF.cor
                      # , ab.cor.part
                      # , pcf.cor.part
                      , EF.cor.spe = EF.spe
                      , ab.cor.spe = ab.spe
                      , pcf.cor.spe = pcf.spe
                      , ab.cor
                      , pcf.cor
                      , ell
                      
    ))
  })
}



tic()
bef_cors <- map_dfr(unique(bef_by_ell$syst), function(syst){
  dat <- bef_by_ell %>% 
    filter( syst == !!syst)
      return(data.frame(fit_lms(dat)
                      , syst = syst
                      , study = dat$study[1]
                      , study_plus = dat$study_plus[1]
   ))
})
toc()

# just over a minute on MBP


# D-EF Correlation graph

pdf("figures/D-EF_cor.pdf", width = 6.5, height = 2.5)
bef_cors %>% 
  group_by(syst) %>% 
  left_join(clr_tab) %>% 
  mutate(study_plus = factor(study_plus
                             , levels = c("tree_carbon"
                                          , "bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"))) %>% 
  ggplot(aes(ell, EF.cor, color = clr))+
  # get ref lines on first so they don't overlap anything
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 0.2) +
  geom_line(size = 0.7) +
  facet_grid(~study_plus) +
  theme_classic()+
  scale_color_viridis_d(option = "plasma") + 
  theme(legend.position = "none"
         , panel.spacing = unit(1, "lines"))+ 
  xlim(-5,5) +
  labs(y = "correlation between \nlog(D) and log(EF)")
dev.off()

# D-ab and D-pcf graph
pdf("figures/EF_cor_partions.pdf", width = 6.5, height = 4.5) 
# note that the height is klugey, should be fixed. 

bef_cors %>% 
  group_by(syst) %>% 
  mutate(best_ell = ell[which.max(abs(.data$EF.cor))]) %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="slope_component") %>%
  pivot_longer(cols = contains("cor")
               , values_to = "correlation"
               , names_to ="correlation_component") %>%
  filter(!grepl("spe", correlation_component ) & !grepl("EF", correlation_component )) %>% 
  left_join(clr_tab) %>% 
  mutate(study_plus = factor(study_plus
                             , levels = c("tree_carbon"
                                          , "bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"))) %>% 
  ggplot(aes(ell, correlation, color = clr))+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 0.2) +
  geom_line(size = 0.7) +
  facet_grid(factor(correlation_component, levels = c("ab.cor", "pcf.cor"))~study_plus) +
  theme_classic()+
  scale_color_viridis_d(option = "plasma") + 
  geom_vline(xintercept = 1, size = 0.2) +
  theme( legend.position = "none"
         , panel.spacing = unit(1.2, "lines"))+ 
  xlim(-5,5)
dev.off()


# just diversity
# D-ab and D-pcf graph
pdf("figures/diversity_profiles.pdf", width = 4, height = 4) 

bef_by_ell %>% 
  group_by(study_plus) %>% 
  mutate(colcol = c(2:4, (1:3)+0.2, (1:3)+0.4)[dense_rank(as.numeric(as.factor(syst)))]) %>% 
  group_by(site, syst, colcol) %>% 
  mutate(ss = paste(site, syst)
         , study_plus = factor(study_plus
                             , levels = c("tree_carbon"
                                          , "bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"))
         ) %>% 
  ggplot(aes(ell, D ))+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 0.2) +
  geom_line(size = 0.2, aes(group =ss, color = colcol)) +
  facet_wrap(~study_plus, scales = "free") +
  theme_classic()+
  scale_color_viridis_c(option = "viridis") + 
  geom_vline(xintercept = 1, size = 0.2) +
  theme( legend.position = "none"
         , panel.spacing = unit(1, "lines"))+ 
  xlim(-5,5) +
  scale_y_log10()
dev.off()





pdf("figures/EF_slope_starting_place.pdf")
# challenge of finding a diverging palette that allows discrimination near middle
bef_cors %>% 
  group_by(syst) %>% 
  mutate(best_ell = ell[which.max(abs(.data$EF.cor))]) %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="slope_component") %>%
  pivot_longer(cols = contains("cor")
               , values_to = "correlation"
               , names_to ="correlation_component") %>%
  filter(!grepl("spe", correlation_component )) %>% 
  ggplot(aes(ell, slope, color = best_ell))+
  geom_point(size = 0.2) +
  facet_grid(study_plus~factor(slope_component, levels = c("EF.slope", "ab.slope", "pcf.slope")), scales = "free") +
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  # labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_gradient2(low = "blue", high = "darkorange2", mid = "grey30") +
  geom_vline(xintercept = 1, size = 0.2) +
  theme( panel.spacing = unit(1.3, "lines")) 
dev.off()

  
# slopes and counterfactuals in blue, black, and red
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

