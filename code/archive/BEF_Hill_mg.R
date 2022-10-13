# code by Mark Genung, Tina Harrison, Michael Roswell to test correlations
# between Hill diversity and ecosystem function across different Hill diversity
# scaling exponent values

library(MeanRarity)
library(ggplot2)

n.datasets = 3 # how many separate datasets to analyze

# abundace data
a.list = readRDS("data/a.list.rds") # readRDS("/Users/Mark/Downloads/a.list.rds")
# function data
z.list = readRDS("data/z.list.rds") # readRDS("/Users/Mark/Downloads/z.list.rds")

bef.r2 = vector("list", n.datasets)
pcf.r2 = vector("list", n.datasets)
D = vector("list", n.datasets)
t.out = vector("list", n.datasets)
c.out = vector("list", n.datasets)

for(i in 1:n.datasets){
  
  # (1) load and format data
  #####
  
  # load data (species (rows) by site (cols) of either abun or function data)
  input.a = a.list[[i]] # abundance 
  input.z = z.list[[i]] # function
  
  # convert column to row names (won't always be needed)
  a = input.a
  rownames(a) = input.a[,1]
  a = a[,-1] # abundance
  
  z = input.z
  rownames(z) = input.z[,1]
  z = z[,-1] # function
  
  # z is already total function in this case, so we need to do this (won't always be needed)
  z = z/a
  z = as.matrix(z)
  z[is.nan(z)] = 0
  
  # TH code
  require(magrittr)
  dim(a)
  dim(z)
  
  S = nrow(a) # gamma
  N = ncol(a) # total site number
  A = colSums(a) # site abundances
  t = colSums(a*z) # site total function
  p = apply(a, 2, function(x) x/sum(x)) # relative abundance matrix
  cwm.func = colSums(p*z) # site cwm function # average per-capita function
  #####
  
  # (2) measure q-diversity
  #####
  
  Ln = function(x) {ifelse(x != 0, log(x), 0)}
  Exp = function(x, y) {ifelse(x != 0, x^y, 0)}
  # use MeanRarity::rarity for this
  # D.nab = function(z, q) {ifelse(q != 1,
  #                                sum(Exp(z/sum(z), q))^(1/(1 - q)),
  #                                exp(-sum((z/sum(z))*Ln(z/sum(z)))) )}
  # example: here is how you get simpson diversity for each site
  # apply(p, 2, D.nab, q = 2) 
  # here is how you get a range of q-diversities for each site
  # first define the range of ell you are interested (change as desired)
  ell_vals = seq(-2, 1, by=0.5) 
  D[[i]] = sapply(ell_vals, function(ell) {
    apply(a, 2, rarity, l = ell)
    }) %>% t
  # D # sites in columns, orders of q in rows
  #####
  
  # for figures
  #####
  
  # save total function and cwm function
  t.out[[i]] = t
  c.out[[i]] = cwm.func
  
  # r2
  bef.r2[[i]] = apply(D[[i]], 1, function(x) {cor(log(x), log(t))^2})
  pcf.r2[[i]] = apply(D[[i]], 1, function(x) {cor(log(x), log(cwm.func))^2})
  #####
}

# main figs
# sorry for using so many lines to do these figs
# select species 1 to n.datasets
# 1= P reptans, 2 = P tanacetifolia, 3 = M fistulosa
se = 3

q.cols =colorspace::heat_hcl(length(ell_vals)) # change if you want different colors

# set theme for ggplots

mytheme<-theme_classic()

# BEF
# refactor this so its programatic. Seems like it should be easy with ggplot. Added advantage of legends
ggplot() +
  stat_smooth(aes(x = D[[se]][1, ]
                  , y = log(t.out[[se]]))
              , method = "lm"
              , se = F
              , col = q.cols[1]
              , size = 2) +
  stat_smooth(aes(x = D[[se]][2, ]
                  , y = log(t.out[[se]]))
              , method = "lm"
              , se = F
              , col = q.cols[2]
              , size = 2) +
  stat_smooth(aes(x = D[[se]][3, ]
                  , y = log(t.out[[se]]))
              , method = "lm"
              , se = F
              , col = q.cols[3]
              , size = 2) +
  stat_smooth(aes(x = D[[se]][4, ]
                  , y = log(t.out[[se]]))
              , method = "lm"
              , se = F
              , col = q.cols[4]
              , size = 2) +
  stat_smooth(aes(x = D[[se]][5, ]
                  , y = log(t.out[[se]]))
              , method = "lm"
              , se = F
              , col = q.cols[5]
              , size=2) +
  stat_smooth(aes(x=D[[se]][6,],y=log(t.out[[se]])),method="lm",se=F,col=q.cols[6],size=2) +
  stat_smooth(aes(x=D[[se]][7,],y=log(t.out[[se]])),method="lm",se=F,col=q.cols[7],size=2) +
  
  geom_point(aes(x=D[[se]][1,], y=log(t.out[[se]])), fill=q.cols[1], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][2,], y=log(t.out[[se]])), fill=q.cols[2], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][3,], y=log(t.out[[se]])), fill=q.cols[3], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][4,], y=log(t.out[[se]])), fill=q.cols[4], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][5,], y=log(t.out[[se]])), fill=q.cols[5], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][6,], y=log(t.out[[se]])), fill=q.cols[6], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  geom_point(aes(x=D[[se]][7,], y=log(t.out[[se]])), fill=q.cols[7], pch=21, 
             stroke=0.2, size=4.5, alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Effective # of Species") +
  scale_y_continuous("Total Function")

# PCF
ggplot() +
  
  stat_smooth(aes(x=D[[se]][1,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[1],size=2) +
  stat_smooth(aes(x=D[[se]][2,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[2],size=2) +
  stat_smooth(aes(x=D[[se]][3,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[3],size=2) +
  stat_smooth(aes(x=D[[se]][4,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[4],size=2) +
  stat_smooth(aes(x=D[[se]][5,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[5],size=2) +
  stat_smooth(aes(x=D[[se]][6,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[6],size=2) +
  stat_smooth(aes(x=D[[se]][7,], y=log(c.out[[se]])),method="lm",se=F,col=q.cols[7],size=2) +
  
  geom_point(aes(x=D[[se]][1,], y=log(c.out[[se]])), fill=q.cols[1], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][2,], y=log(c.out[[se]])), fill=q.cols[2], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][3,], y=log(c.out[[se]])), fill=q.cols[3], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][4,], y=log(c.out[[se]])), fill=q.cols[4], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][5,], y=log(c.out[[se]])), fill=q.cols[5], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][6,], y=log(c.out[[se]])), fill=q.cols[6], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  geom_point(aes(x=D[[se]][7,], y=log(c.out[[se]])), fill=q.cols[7], pch=21,
             stroke=0.2,size=4.5,alpha=0.5) +
  
  mytheme +
  scale_x_continuous("Effective # of Species") +
  scale_y_continuous("Per-Capita Function")

# these plots are for all three species so you don't need to pick one

# r2 plot - BEF
ggplot() +
  
  geom_line(aes(x=ell_vals, y=bef.r2[[1]]), size=2) +
  geom_point(aes(x=ell_vals, y=bef.r2[[1]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  geom_line(aes(x=ell_vals, y=bef.r2[[2]]), size=2) +
  geom_point(aes(x=ell_vals, y=bef.r2[[2]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  geom_line(aes(x=ell_vals, y=bef.r2[[3]]), size=2) +
  geom_point(aes(x=ell_vals, y=bef.r2[[3]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  mytheme +
  scale_x_continuous("ell") +
  scale_y_continuous("BEF r2")

# r2 plot - PCF
ggplot() +
  
  geom_line(aes(x=ell_vals, y=pcf.r2[[1]]), size=2) +
  geom_point(aes(x=ell_vals, y=pcf.r2[[1]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  geom_line(aes(x=ell_vals, y=pcf.r2[[2]]), size=2) +
  geom_point(aes(x=ell_vals, y=pcf.r2[[2]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  geom_line(aes(x=ell_vals, y=pcf.r2[[3]]), size=2) +
  geom_point(aes(x=ell_vals, y=pcf.r2[[3]]), fill=q.cols, pch=21, stroke=1, size=8) +
  
  mytheme +
  scale_x_continuous("ell") +
  scale_y_continuous("PCF r2")

