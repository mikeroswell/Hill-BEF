library(tidyverse)
library(MeanRarity)




# For these 7 SADs, Diversity is always the same at ell = 1 (5) and abundance is always 13.
# for ell -> infinity, diversity is same for more, flattened, and squash, dom, and sub. It is
# also identical for ev, norare

# for ell -> -Inf, diversity is same for more, norare, and sub; for flattened, ev, and
# squash is intermediate and sub and dom are way low.

# Shannon and Simpson are both good at differentiating these SADS, and no rank changes in there. 

# ell = 1.5 doesn't distinguish more from squash (nor any ell>1), or flattened from norare (these have different max rarities tho)
more<-c(5,3,3,1,1)
flattened <- c(3,3,3,3,1)
norare <- c(5,2,2,2,2)
ev <- c(3,3,3,2,2)
squash<-c(4,4,3,1,1)
dom<-c(9,1,1,1,1)
sub<-c(5,5,1,1,1)
# rarity(more, l = 1.5)
# rarity(flattened, l = 1.5)
# rarity(squash, l =1.5)
# rarity(more, l = 0.5)
# rarity(flattened, l = 0.5)
# rarity(squash, l =0.5
SADS <- c("more"
          , "flattened"
          , "squash"
          , "ev"
          , "norare"
          , "dom"
          , "sub"
)

# SADS <- c("rich", "poor", "bigpoor")

dpros<- map_dfr(SADS, function(SAD){
  df<-divpro(get(SAD), ell_low = -250, ell_hi = 250, by = 0.01)
  return(data.frame(df, SAD))
})

dpros %>% ggplot(aes(ell, d, color = SAD))+
  geom_point(size = 0.1)+
  theme_classic() +
  geom_vline(xintercept = 0, lty = 5) +
  geom_vline(xintercept = -1, lty = 3) +
  geom_vline(xintercept = 1, lty = 1)  +
  geom_vline(xintercept = 1.5, lty = 2)  +
  xlim(c(-2, 2))
  

dpros %>% ggplot(aes(ell, d, color = SAD))+
  geom_point(size = 0.1)+
  theme_classic() +
  # geom_vline(xintercept = 0, lty = 5) +
  # geom_vline(xintercept = -1, lty = 3) +
  # geom_vline(xintercept = 1, lty = 1)  +
  geom_vline(xintercept = 1.5, lty = 2) # +
  # xlim(c(-2, 2))


# make some very simple examples to build some intuitions
# different richnesses

rich<- c(rep(5, 7),1)
rarity(rich, 1.5)
poor <- c(31, rep(1,4))
rarity(poor, 1.5)

# big vs small
big <- c(40, 30, 25, 20, 18, 16, 15, 1)
rarity(big, 1.5)
bigpoor <- c(60,50, 30, 20, 1)
rarity(bigpoor, 1.5)

# I have this intuition that skew is what matters, but not comfortable enough with this to pin it down
rarity(c(more, 1), 1.5)
rarity(c(more, 3), 1.5)
rarity(c(more, 5), 1.5)
# adding one hyperdominant (which makes existing rare ultrarare)
# is like adding to more rare spp
rarity(c(more, 19), 1.5)
rarity(c(more, 1,1), 1.5)
