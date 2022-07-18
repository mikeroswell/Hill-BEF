# analyze data with Tina partition
source("code/format_BEF_data.R")

library(MeanRarity)
library(furrr)
# loop over ell

cubeRt <- function(x){sign(x)*abs(x)^(1/3)}

ell_vec<-seq(-10, 10, 0.2)
ell_vec <- c(exp(seq(-12,3, 0.1)), -(exp(seq(-12,3, 0.1))))
ell_vec <- cubeRt(c(seq(-80, 80, 0.5)
                    , seq(-0.5, 0.5, 0.01)))

# set number of "cores" to use for parallelization
# on intel mac at least, hyperthreading is used so the maximum is double the 
# number of physical cores
nc<-7

plan(strategy = "multiprocess", workers = nc)

# this could use some work to become more general but works fine for lefcheck
# remember that this works on groups only if they exist (e.g. site)
sum_by_ell <- function(dat
                       , ell_vec
                       , abundance = "Abundance"
                       , EF = "Biomass"
                       ,  ...){
  furrr::future_map_dfr(ell_vec, function(ell){
    dat %>% summarize(ell =ell
                       , D = rarity(eval(as.name(abundance)), ell)
                       , ab = sum(eval(as.name(abundance)))
                       , EFt = sum(eval(as.name(EF)))
                       , pcf = EFt/ab
                       # , meanlat = mean(abs(SiteLat))
                     )})
}

# res2<-sum_by_ell(ab_df %>% group_by(site) %>%  mutate(SiteLat = as.numeric(site)), ell_vec = ell_vec)
fores <- sum_by_ell(ssf %>% group_by(loc, X1ha.Plot.Number), ell_vec = ell_vec, EF = "totCarbon")
# fores %>% 
#   ggplot(aes(ell, D, color = X1ha.Plot.Number)) +
#   geom_point(size = 0.2) +
#   theme_classic()+
#   scale_color_viridis_d()+
#   scale_y_log10() +
#   xlim(c(-10, 10))

fores %>% filter(ell == -10) %>% ggplot(aes(ab, pcf, color = loc)) + geom_point()+theme_classic() + scale_x_log10()+scale_y_log10()
fores %>% filter(ell == -10) %>% ggplot(aes(ab, D, color = loc)) + geom_point()+theme_classic()+ scale_x_log10()+scale_y_log10()+facet_wrap(~loc, scales = "free")

ssf %>% group_by()


res<-sum_by_ell(lefcheck_by_site, ell_vec = ell_vec) 

# res %>% filter(meanlat > 45) %>% 
#   ggplot(aes(ell, D, color = meanlat)) +
#   geom_point(size = 0.2) +
#   theme_classic()+
#   scale_color_viridis_c()+
#   scale_y_log10() +
#   xlim(c(-10, 2))


# # check my work
# res %>% ggplot(aes(ell, D, color = Province))+
#   geom_point(alpha = 0.2)+
#   theme_classic()


fit_lms <- function(sub, dataset = "lefcheck"){
  future_map_dfr(ell_vec, function(ell){
    dat = sub %>% filter(ell == !!ell) 
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
    # ab.cor.part = ab.slope* sd(Ln(dat$D))/(sd(
    #   mean(Ln(dat$pcf)) * Ln(dat$ab)
    # ))
    # pcf.cor.part = pcf.slope * sd(Ln(dat$D))/(sd(
    #   mean(Ln(dat$ab)) * Ln(dat$pcf)
    # ))
    # 
    ab.cor.raw = cor(Ln(dat$ab), Ln(dat$D))
    pcf.cor.raw = cor(Ln(dat$D), Ln(dat$pcf))
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
                      , ab.cor.raw
                      , pcf.cor.raw
                      , ell
                      # , Latitude = dat$meanlat
                      , dataset = dataset
    ))
  })
}


forest_out <- map_dfr(unique(locs), function(loc){
  sub <- fores %>% filter(loc== !!loc)
  data.frame(fit_lms(sub = sub, dataset = "TEAM"), loc = loc)
})

forest_out


pdf("figures/quickDirty")
forest_out %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="component") %>%
  pivot_longer(cols = contains("cor")
               , values_to = "correlation"
               , names_to ="correlation_component") %>%
  filter(!grepl("spe", correlation_component )) %>% 
  group_by(loc) %>% 
  ggplot(aes(ell, correlation, color = loc))+
  geom_point(size = 0.2) +
  facet_wrap(~correlation_component) +
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  # labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_d() +
  geom_vline(xintercept = 1, size = 0.2) # +
  # scale_x_continuous(trans = scales::trans_new("cuberoot"
  #                                              , function(x){sign(x)*(abs(x))^(1/3)}
  #                                              , inverse = function(x){x^3}
  #                                              , breaks = function(x){round(x,2)})
  #                    , breaks = c(-5, -3, -1, -0.5, -0.1, 0.1,  0.5, 1,  3,5 )) +
  # theme(axis.text.x = element_text(angle = 90))
dev.off()
# sim_out <-fit_lms(res2, dataset = "sim1")

first_out <- map_dfr(unique(res$Province), function(Province){
  sub <- res %>% filter(Province == !!Province)
  data.frame(fit_lms(sub = sub), Province = Province)
  })

# get highest correlations

highells<-res %>% 
  mutate(highell = Province %in% 
           c(first_out %>% 
               group_by(Province) %>% 
               filter(EF.cor == max(EF.cor)) %>% 
               ungroup() %>% 
               filter(ell>2) %>% 
               pull(Province))
  )


CVab<- highells %>% 
  filter(ell ==1) %>% 
  group_by(highell, Province) %>% 
  summarize(cvAbund = sd(ab)/mean(ab)
           , cvRich = sd(D)/mean(D) 
           , cvPCF = sd(pcf)/mean(pcf)) 



CVab %>% ggplot(aes(cvRich, cvPCF, color = highell))+
  # geom_smooth() +
  geom_point()+
  theme_classic()



pdf("figures/highell_from_richness_and_abundance_variation.pdf")

CVab %>% ggplot(aes(cvAbund, cvRich, color = highell))+
  # geom_smooth() +
  geom_point()+
  theme_classic()

dev.off()

pdf("figures/highell_bc_high_cv_abund.pdf")

CVab %>% ggplot(aes( highell, cvAbund))+
  geom_boxplot()+
  geom_beeswarm(size = 2)+
  theme_classic() +
  labs(x = "high ell most predictive")

dev.off()

highells_raw <- lefcheck_by_site %>% 
  mutate(highell = Province %in% 
           c(first_out %>% 
               group_by(Province) %>% 
               filter(EF.cor == max(EF.cor)) %>% 
               ungroup() %>% 
               filter(ell>2) %>% 
               pull(Province)))

pdf("figures/look_for_humps.pdf")
highells_raw %>% group_by(SiteCode) %>% 
  ggplot(aes(log(Abundance), log(Biomass/Abundance), color =highell, shape = SiteCode)) +
  geom_smooth()+
  facet_wrap(~Province)+
  theme_classic()+
  theme(legend.position = "none") 
    
dev.off()
# for now save data so I can restart computer
# write.csv(first_out, "data/lefcheck_results.csv", row.names = F)
# first_out<-read.csv("data/lefcheck_results.csv")

pdf("figures/rank_correlations.pdf")
pdf("figures/cuberoot_correlations.pdf")
first_out %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="component") %>%
  pivot_longer(cols = contains("cor")
               , values_to = "correlation"
               , names_to ="correlation_component") %>%
  filter(!grepl("spe", correlation_component )) %>% 
  group_by(Province, Latitude) %>% 
  ggplot(aes(ell, correlation, color = Latitude))+
  geom_point(size = 0.2) +
  facet_wrap(~correlation_component) +
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_c(option = 
                          "plasma") +
  geom_vline(xintercept = 1, size = 0.2) +
  scale_x_continuous(trans = scales::trans_new("cuberoot"
                                       , function(x){sign(x)*(abs(x))^(1/3)}
                                       , inverse = function(x){x^3}
                                       , breaks = function(x){round(x,2)})
                     , breaks = c(-5, -3, -1, -0.5, -0.1, 0.1,  0.5, 1,  3,5 )) +
  theme(axis.text.x = element_text(angle = 90))
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

sim_out %>% pivot_longer(cols = ends_with(".slope")
                         , values_to = "slope"
                         , names_to ="component") %>% 
  ggplot(aes(ell, slope, color = component))+
  geom_line()+
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 0.2) +
  scale_color_manual(values = c("red", "black", "blue")) 

# slopes and correlation coefficients from all "Provinces" together
pdf("figures/lefcheck_together.pdf")
# absolute value of slope
first_out %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="component") %>% 
  group_by(Province, Latitude) %>% 
  ggplot(aes(ell, abs(slope), color = Latitude, shape = Province))+
  geom_line(size = 1, alpha = 0.7)+
  # geom_smooth()+
  facet_wrap(~component) +
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_c(option = 
                        "plasma") +
  geom_vline(xintercept = 1, size = 0.2)

# raw slope
first_out %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="component") %>% 
  group_by(Province, Latitude) %>% 
  ggplot(aes(ell, slope, color = Latitude, shape = Province))+
  geom_line(size = 1, alpha = 0.7)+
  # geom_smooth()+
  facet_wrap(~component) +
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_c(option = 
                          "plasma") +
  geom_vline(xintercept = 1, size = 0.2)

# absolute value of r
first_out %>% group_by(Province, Latitude) %>% 
  ggplot(aes(ell, abs(EF.cor), color = Latitude, shape = Province)) +
  geom_line(size = 1, alpha = 0.7)+
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_c(option = 
                          "plasma") +
  geom_vline(xintercept = 1, size = 0.2)
# raw correlation
first_out %>% group_by(Province, Latitude) %>% 
  ggplot(aes(ell, EF.cor, color = Latitude, shape = Province)) +
  geom_line(size = 1, alpha = 0.7)+
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  labs(color = "mean absolute \nlatitude for \nprovince")+
  scale_color_viridis_c(option = 
                          "plasma") +
  geom_vline(xintercept = 1, size = 0.2)
dev.off()  

# play with BCI, this won't make sense until subplots are created though. 
bci_out<-bci_sum %>% sum_by_ell(ell_vec = ell_vec)


bci_slopes <- future_map_dfr(ell_vec, function(ell){
  dat = bci_out %>% filter(ell == !!ell)
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
                    , dataset = "bci"
  ))
})

pdf("figures/bci_partitions_lumped.pdf")
bci_slopes %>% 
  pivot_longer(cols = ends_with(".slope")
               , values_to = "slope"
               , names_to ="component") %>% 
  ggplot(aes(ell, slope, color = component))+
  geom_line()+
  theme_classic()+
  geom_hline(yintercept = 0, size = 0.2)  +
  geom_vline(xintercept = 1, size = 0.2) +
  scale_color_manual(values = c("red", "black", "blue")) +
  ggtitle("BCI")

dev.off()