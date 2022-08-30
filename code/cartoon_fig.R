# Script to generate cartoon figure of Hill number profiles highlighting 
# - use of ell
# - ell > 1

# load packages
library(MeanRarity)
library(tidyverse)
library(patchwork)

# example abundance vectors. One has high abundance and rel. low dominance
com1 = c(80, 50, 20, 5, 3, NA, NA)
# other has higher richness, but also higher dominance and lower abundance
com2 = c(20, 1, 1, 1, 1, 1, 1)

# vector of hill scaling param
ell_vals = seq(-2, 5, by=0.05)

# compute diversity profiles with MeanRarity
D1 = MeanRarity::divpro(com1, -3, 5, 0.05)
D2 = MeanRarity::divpro(com2, -3, 5, 0.05)

D.dat <- bind_rows(D1, D2, .id = "com")
a.dat <- data.frame(com1, com2, sp = 1:length(com1)) %>% 
  pivot_longer(cols = c("com1", "com2")
               , names_to = "com"
               , names_prefix = "com"
               , values_to = "abund")


# barplot for RADs 
bp <- a.dat %>% ggplot(aes(sp, abund, fill = as.factor(com))) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(x = "species rank", y = "species abundance" ) + 
  guides(fill = "none") 


# diversity profiles on same axes
dp <- D.dat %>% 
  ggplot(aes(ell, d, color = com)) +
  geom_vline(xintercept = c(-1:1), linetype = c(3:1)) +
  geom_line(size = 1.2) +
  labs(y = "Hill diversity", color = "community") +
  theme_classic()

# print to .pdf for full-page spread using patchwork
pdf("figures/cartoon_fig_horizontal.pdf", height = 2.5, width = 6.5)
bp + dp
dev.off()
