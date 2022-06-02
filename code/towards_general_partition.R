# analyze data with Tina partition
source("code/format_BEF_data.R")

library(MeanRarity)
library(furrr)
# loop over ell

ell_vec<-seq(-10, 10, 0.2)

# set number of "cores" to use for parallelization
# on intel mac at least, hyperthreading is used so the maximum is double the 
# number of physical cores
nc<-7

plan(strategy = "multiprocess", workers = nc)

# this could use some work to become more general but works fine for lefcheck
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
                       , meanlat = mean(abs(SiteLat))
                     )})
}

# res2<-sum_by_ell(ab_df %>% group_by(site) %>%  mutate(SiteLat = as.numeric(site)), ell_vec = ell_vec)

res<-sum_by_ell(lefcheck_by_site, ell_vec = ell_vec) 


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
                      , Latitude = dat$meanlat
                      , dataset = dataset
    ))
  })
}


# sim_out <-fit_lms(res2, dataset = "sim1")

first_out <- map_dfr(unique(res$Province), function(Province){
  sub <- res %>% filter(Province == !!Province)
  data.frame(fit_lms(sub = sub), Province = Province)
  })

# for now save data so I can restart computer
# write.csv(first_out, "data/lefcheck_results.csv", row.names = F)
# first_out<-read.csv("data/lefcheck_results.csv")

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