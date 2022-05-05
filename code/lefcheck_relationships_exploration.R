lsum<-lefcheck %>% 
  group_by(Province, SiteCode) %>% 
  summarize(amu = mean(Abundance), avar = sd(Abundance), rich = n_distinct(SPECIES_NAME), fn = sum(Biomass))


pdf(file = "figures/abundance_richness_lefcheck.pdf", width = 10, height = 7.5)
lsum %>% ggplot(aes(amu, rich)) + 
  geom_point() +
  facet_wrap(~Province, scales = "free") +
  labs(x = "mean species abundance", y = "site richness") +
  theme_classic()
dev.off()

pdf(file = "figures/sd_abundance_richness_lefcheck.pdf", width = 10, height = 7.5)
lsum %>% ggplot(aes(avar, rich)) + 
  geom_point() +
  facet_wrap(~Province, scales = "free") +
  labs(x = "sd abundance", y = "site richness") +
  theme_classic()
dev.off()

pdf(file = "figures/abundance_meanvar_lefcheck.pdf", width = 10, height = 7.5)
lsum %>% ggplot(aes(amu, avar)) + 
  geom_point() +
  facet_wrap(~Province, scales = "free") +
  labs(x = "mean species abundance", y = "sd abundance") +
  theme_classic()
dev.off()

pdf(file = "figures/abundance_EF_lefcheck.pdf", width = 10, height = 7.5)
lsum %>% ggplot(aes(amu, fn)) + 
  geom_point() +
  facet_wrap(~Province, scales = "free") +
  labs(x = "mean species abundance", y = "EF") +
  theme_classic()
dev.off()
  
