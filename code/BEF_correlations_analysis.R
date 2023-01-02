# analyses of correlations between log diversity, log function, and also some 
# relationships between those two variables and abundance and per-capita 
# function.

# load packages. 
# devtools::install_github("mikeroswell/MeanRarity")
library(MeanRarity) 
library(tidyverse)
library(furrr)
library(tictoc)
library(egg)

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
nc <- 9



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

write.csv(bef_by_ell, "data/bef_by_ell.csv", row.names = FALSE)

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

bef_by_ell <- read.csv("data/bef_by_ell.csv")

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


write.csv(bef_cors, "data/bef_correlations.csv", row.names = FALSE)
bef_cors <- read.csv("data/bef_correlations.csv")
# D-EF Correlation graph

lab_clean <- function(x){str_replace(x, "_", " ") %>% 
                                  str_replace("_", " \\*")}

quartz(file = "figures/Fig3_D-EF_cor.pdf"
       , width = 6.5
       , height = 2.5
       , type = "pdf")
f3 <- bef_cors %>% 
  group_by(syst) %>% 
  left_join(clr_tab) %>% 
  mutate(study_plus = factor(study_plus
                             , levels = c("bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon"))) %>%  
  ggplot(aes(ell, EF.cor, color = clr))+
  # get ref lines on first so they don't overlap anything
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 1) +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  geom_line(size = 0.7) +
  facet_grid(~study_plus
             , labeller = labeller(study_plus = lab_clean)) +
  theme_classic()+
  # scale_color_manual(values = hcl.colors(16, "hawaii")) +
  scale_color_viridis_d(option = "plasma") +
  theme(legend.position = "none"
        , panel.spacing = unit(1.1, "lines")
        , strip.background = element_blank()
        , strip.placement = "outside"
        # , strip.text.x = element_blank() #element_text(vjust = 1, hjust = 0.3 )
        )+ 
  xlim(-5,5) +
  labs(x = "Hill diversity scaling factor, \"\U2113\""
       , y = "\nBEF correlation") 

tag_facet(f3, open = NULL, close = NULL) + 
  theme(strip.text.x = element_text(size = rel(0.8)
                                    , margin = margin(b = 5)))

dev.off()

quartz(file = "figures/Fig3_D-EF_cor_SPEARMAN.pdf", width = 6.5, height = 2.5, type = "pdf")
f3 <- bef_cors %>% 
  group_by(syst) %>% 
  left_join(clr_tab) %>% 
  mutate(study_plus = factor(study_plus
                             , levels = c("bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon"))) %>%  
  ggplot(aes(ell, EF.cor.spe, color = clr))+
  # get ref lines on first so they don't overlap anything
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 1) +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  geom_line(size = 0.7) +
  facet_grid(~study_plus
             , labeller = labeller(study_plus = lab_clean)) +
  theme_classic()+
  # scale_color_manual(values = hcl.colors(16, "hawaii")) +
  scale_color_viridis_d(option = "plasma") +
  theme(legend.position = "none"
        , panel.spacing = unit(1.1, "lines")
        , strip.background = element_blank()
        # , strip.text.x = element_blank() #element_text(vjust = 1, hjust = 0.3 )
  )+ 
  xlim(-5,5) +
  labs(x = "Hill diversity scaling factor \"\U2113\""
       , y = "\nBEF correlation") 

tag_facet(f3, open = NULL, close = NULL) + 
  theme(strip.text.x = element_text(size = rel(0.8)
                                    , margin = margin(b = 5)))

dev.off()


# D-ab and D-pcf graph
quartz(file = "figures/Fig4_EF_cor_partions.pdf"
       , width = 6.5
       , height = 4.8
       , type = "pdf") 
# note that the height is klugey, should be fixed. 

f4 <- bef_cors %>% 
  filter(ell>-5 & ell <5) %>% 
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
                             , levels = c("bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon"))
         , correlation_component = factor(correlation_component
                                          , levels = c("ab.cor", "pcf.cor")
                                          , labels = c("abundance \ncorrelation", "per-capita function \ncorrelation"))) %>%  

  ggplot(aes(ell, correlation, color = clr))+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 1) +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  geom_line(size = 0.7) +
  facet_grid(correlation_component ~ study_plus
             , labeller = labeller(study_plus = lab_clean)
             , switch = "y") +
  theme_classic()+
  scale_color_viridis_d(option = "plasma") + 
  theme( legend.position = "none"
         , panel.spacing = unit(1.1, "lines")
         , strip.text.y = element_text(size = rel(1))
         , strip.placement = "outside"
         , strip.background = element_blank()
         , axis.title.y = element_blank()
         ) + 
  labs(x = "Hill diversity scaling factor, \"\U2113\"") +
  xlim(-5,5)

tag_facet(f4, open = NULL, close = NULL) + 
  theme(strip.text.x = element_text(size = rel(0.8)
                                    , margin = margin(b = 5)))
dev.off()


# bef_cors %>% 
#   group_by(syst) %>% 
#   mutate(best_ell = ell[which.max(abs(.data$EF.cor))]) %>% 
#   filter(best_ell < -4) %>% 
#   pull(syst)


bef_data %>% 
  filter(syst == "Cold Temperate Northeast Pacific_TRUE") %>% 
  group_by(site) %>% 
  summarize(rich =n(), totAb = sum(abund)) %>% 
  arrange(totAb)
  
  # ggplot(aes(totAb, rich))+
  # geom_point()+
  # theme_classic()
  # 
bef_by_ell %>% 
  filter(syst == "Cold Temperate Northeast Pacific_TRUE") %>% 
  filter(ell %in% c(-10, -1, 0, 1, 1.5, 2, 5, 10)) %>% 
  ggplot(aes(ab, D))+
  geom_point()+
  facet_wrap(~ell, scales = "free")+
  theme_classic()

bef_by_ell %>% 
  filter(syst == "Cold Temperate Northeast Pacific_TRUE") %>% 
  filter(ell %in% c(-10, -1, 0, 1, 1.5, 2, 5, 10)) %>% 
  ggplot(aes(log(D), log(EFt)))+
  geom_smooth(method = "lm")+
  geom_point()+
  facet_wrap(~ell, scales = "free")+
  theme_classic()

bef_cors %>% 
  filter(syst == "Cold Temperate Northeast Pacific_TRUE") %>%
  ggplot(aes(ell, EF.cor))+
  geom_point()+
  theme_classic()

bef_cors %>% 
  group_by(syst, study_plus) %>% 
  summarize(best_ell = ell[which.max(abs(.data$EF.cor))]) %>% 
  filter(abs(best_ell)<2) %>% 
  ungroup() %>% 
  summarize(mean(best_ell))

pdf("figures/Fig2_variant_best_ell_histogram_with_r2_shading.pdf"
    , width = 4.5, height = 3.5)

bef_cors %>%
  mutate(study_plus = factor(study_plus 
                           , levels = c("bee_pollination"
                                         , "reef_fish_temperate"
                                         , "reef_fish_tropical"
                                         , "tree_carbon")
                           , labels = lab_clean(c("bee_pollination"
                                                 , "reef_fish_temperate"
                                                 , "reef_fish_tropical"
                                                 , "tree_carbon")))) %>% 
  group_by(syst, study_plus) %>% 
  summarize(best_ell = ell[which.max(abs(.data$EF.cor))]
            , best_R2 = (EF.cor[which.max(abs(.data$EF.cor))])^2) %>% 
  ggplot(aes(x = best_ell
             , fill = study_plus
             , color = study_plus
             , alpha = best_R2
             , group = interaction(best_R2, study_plus))) +
  geom_histogram()+
  geom_hline(yintercept = 0, size = 0.4, linetype = "solid")+
  geom_vline(xintercept = 1, size = 1, linetype = "solid") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  theme_classic() +
  labs(x = "Hill diversity scaling factor \U2113 \nwith highest BEF correlation"
       , y = "community datasets", fill = "system", alpha = "best_R2" ) +
  guides(color = "none")
dev.off()

quartz(file = "figures/Fig2_solid_best_ell_histogram.pdf"
       , width = 4.5
       , height = 3.5
       , type = "pdf")
bef_cors %>%
  mutate(study_plus = factor(study_plus 
                             , levels = c("bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon")
                             , labels = lab_clean(c("bee_pollination"
                                                    , "reef_fish_temperate"
                                                    , "reef_fish_tropical"
                                                    , "tree_carbon")))) %>% 
  group_by(syst, study_plus) %>% 
  summarize(best_ell = ell[which.max(abs(.data$EF.cor))]
            , best_R2 = (EF.cor[which.max(abs(.data$EF.cor))])^2) %>% 
  ggplot(aes(x = best_ell
             , fill = study_plus
             # , alpha = best_R2
             , group = interaction(best_R2, study_plus))
         ) +
  geom_histogram(color = "black")+
  # geom_hline(yintercept = 0, size = 0.2, linetype = "solid")+
  geom_vline(xintercept = 1, size = 1, linetype = "solid") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  theme_classic() +
  labs(x = "Hill diversity scaling factor, \"\U2113\" \nwith highest BEF correlation"
       , y = "community datasets", fill = "system"
       # , alpha = "best_R2" 
       ) +
  guides(color = "none")
dev.off()


pdf("figures/Fig2_variant_best_ell_biplot.pdf", width = 6, height = 3.5)
bef_cors %>% 
  group_by(syst, study_plus) %>% 
  mutate(study_plus = factor(study_plus 
                             , levels = c("bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon")
                             , labels = lab_clean(c("bee_pollination"
                                                    , "reef_fish_temperate"
                                                    , "reef_fish_tropical"
                                                    , "tree_carbon")))) %>% 
  summarize(best_ell = ell[which.max(abs(.data$EF.cor))]
            , best_R2 = (EF.cor[which.max(abs(.data$EF.cor))])^2) %>% 
  ggplot(aes(x = best_ell
             , y = best_R2
             , fill = study_plus
            
             , shape = study_plus
             # , alpha = best_R2
             # , group = interaction(best_R2, study_plus)
             ) , color = "black") +
  geom_point(size = 2.4) +
  geom_vline(xintercept = 1, size = 1, linetype = "solid") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  scale_shape_manual(values = 21:24) + 
  theme_bw() +
  theme(panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , strip.background = element_blank()) +
  labs(x = "Hill diversity scaling factor, \"\U2113\" \nwith highest BEF correlation"
       , y = bquote('BEF '*r^2*''), shape = "system", fill = "system" ) 
dev.off()
# just diversity
# D-ab and D-pcf graph
quartz(file = "figures/S1_diversity_profiles.pdf"
       , width = 4.5
       , height = 4
       , type = "pdf") 

bef_by_ell %>% 
  group_by(study_plus) %>% 
  mutate(colcol = c(2:4, (1:3)+0.2, (1:3)+0.4)[dense_rank(as.numeric(as.factor(syst)))]) %>% 
  group_by(site, syst, colcol) %>% 
  mutate(ss = paste(site, syst)
         , study_plus = factor(study_plus
                             , levels = c( "bee_pollination"
                                          , "reef_fish_temperate"
                                          , "reef_fish_tropical"
                                          , "tree_carbon"))
         ) %>% 
  ggplot(aes(ell, D ))+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_line(size = 0.2, aes(group =ss, color = colcol)) +
  facet_wrap(~study_plus, scales = "free") +
  theme_classic()+
  scale_color_viridis_c(option = "viridis") + 
  geom_vline(xintercept = 1, size = 1, linetype = "solid") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, size = 1, linetype = "dotted") +
  theme( legend.position = "none"
         , panel.spacing = unit(1, "lines"))+ 
  xlim(-5,5) +
  scale_y_log10() +
  labs(x = "Hill diversity scaling factor, \"\U2113\"", y = "Hill diversity")
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

# quick data explorations to support some results text
# what is the range of correlations for ell>1 for bee data?

# here, looking at range of R2 for D-EF relnshp for ell>1
bef_cors %>% 
  filter( ell >1) %>% 
  group_by(study_plus, syst) %>% 
  summarize(R2.max = max(EF.cor^2)
            , R2.mean = mean(EF.cor^2)
            , R2.min = min(EF.cor^2)
  ) %>% 
  pivot_longer(cols = starts_with("R2")
               , names_prefix = "R2\\."
               , names_to = "quant"
               , values_to = "val") %>% 
  ggplot(aes(val, fill = study_plus))+
  geom_histogram() +
  facet_wrap(~quant) +
  theme_classic()


# EF-D for negative ell (raw correlation)
bef_cors %>% 
  filter( ell <0) %>% 
  group_by(study_plus, syst) %>% 
  summarize(R.max = max(EF.cor)
            , R.mean = mean(EF.cor)
            , R.min = min(EF.cor)
  ) %>% 
  pivot_longer(cols = starts_with("R")
               , names_prefix = "R\\."
               , names_to = "quant"
               , values_to = "val") %>% 
  ggplot(aes(val, fill = study_plus))+
  geom_histogram() +
  facet_wrap(~quant) +
  theme_classic()


# ab-D R2
bef_cors %>% 
  ggplot(aes(ell, ab.cor^2, color = study_plus))+ 
  geom_point() +
  theme_classic() +


# T.test stuff
bef_cors %>% filter(ell %in% c(-10, -1, 0, 1, 1.5, 10)) %>% 
  group_by(ell) %>% 
  summarize(mean_R2 = mean(EF.cor^2)
            , lcl = t.test(EF.cor)$conf.int[1]
            , ucl = t.test(EF.cor)$conf.int[2]
            , two_sided_p = t.test(EF.cor)$p.value)

