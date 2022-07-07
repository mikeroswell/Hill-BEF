# Play with simple data to get a picture of BEF magic near richness
library(tidyverse)
library(MeanRarity)
# aba <- c(1, 1,1, 1, 1,1,1,1,1, 19)
# # aba <- c(4,4,4,4,4,4,4,4,4,4,10)
# # abb <- c( 1,2, 2, 2, 2, 2)
# 
# abb <- c( 2,2, 2, 2, 2, 3,3,3,3,3,3)
# abc <- c(1,1,1,1,1, 2, 3, 4, 6)
# 
# 
# aba <- c(4,4,4,4,4,5)
# abb <- c(3,3,3,3,4,8)
# abc <- c( 1, 2, 2, 2, 11)


aba <- c(rep(c(2,5,25), 18),1)
abb <- c(rep(70, 37),5, 222, rep(14, 5))
abc <- c( 1, 2, 2, 2, 11)
sum(aba)
sum(abb)/5

577/25


length(abb)
length(aba)
comms <- list(aba, abb, abc)



Dpro<-map_dfr(seq(-10, 100, 0.1), function(ell){
  map_dfr(1:3, function(comm){
    data.frame(D = rarity(comms[[comm]], ell)
               , ell
               , comm
               , EF = c(1,2,3)[comm])
  })
})

Dpro %>% 
  select(D, ell, comm) %>%  
  pivot_wider(names_from = comm, values_from= D) %>% 
  mutate(comm_1_more_diverse = `1` > `2`) %>% 
  ggplot(aes(ell, as.numeric(comm_1_more_diverse)) )+
  geom_point() +
  theme_classic()

Dpro %>% 
  ggplot(aes(ell, D, color = comm)) + 
  geom_point(size = 0.2)+
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = 0.5, color = "blue") +
  geom_vline(xintercept = 0, color = "red") +
  theme_classic() +
  scale_y_log10(limits= c(10, 680)) 
 


Dpro %>% group_by(ell) %>% 
  summarize(BEF = cor(D, EF)) %>% 
  ggplot(aes(ell, BEF)) + 
  geom_point() +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = 0.5, color = "blue") +
  geom_vline(xintercept = 0, color = "red") +
  theme_classic()
