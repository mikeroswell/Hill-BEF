# analysis file for Hill-BEF relationships with forests and bees
# Mark Genung, Michael Roswell

library(gratia) 
library(tidyverse)
library(MeanRarity)
library(mgcv)
# remotes::install_github("TimTeaFan/dplyover")
# library(dplyr)
# choose site and q values method
#####
site              = "YA"                               # c("PA", "BC", "YA")
ell.seq             = seq(from= -5.0, to= 5.0, by=0.20)  # range of q (Hill numbers)
area.interval.vec = c(25,33,50,100)                 # don't change
# how does area interval work? e.g. area interval 20 gives a 5x5 grid (100/20 = 5)
#####

# big loop, each time subdivides the plots to different sizes

# create output matrices and lists
ef.cor.out  = matrix(NA, length(area.interval.vec), length(ell.seq))
pcf.cor.out = matrix(NA, length(area.interval.vec), length(ell.seq))
div.out     = vector("list", length(area.interval.vec))
c.ef.out    = vector("list", length(area.interval.vec))
c.pcf.out   = vector("list", length(area.interval.vec))

for(are.int in 1:length(area.interval.vec)){
  
  # data reading and processing at the beginning of the big loop
  # probably would be more efficient outside
  start = read.csv(paste("data/", site, "_biomass.csv",
                         sep=""))
  
  if(site=="PA"){start = start[which(start$Year==2012),]}
  start$abun = 1
  colnames(start)[which(colnames(start)=="X1ha.Plot.Number")] = "site"
  
  colnames(start)[which(colnames(start)=="X1ha.Plot.X.Coordinate")] = "x"
  colnames(start)[which(colnames(start)=="X1ha.Plot.Y.Coordinate")] = "y"
  
  area.interval = area.interval.vec[are.int]
  
  # set up for area subdivision loop
  #####
  start$small.site = NA
  
  area.vec = seq(from=0, to=100, by=area.interval)
  area.vec = area.vec[-length(area.vec)] # we don't want the last value
  area.vec
  
  area.grid = expand.grid(area.vec, area.vec)
  #####
  
  # subdivide plots by area
  #####
  for(i in 1:nrow(area.grid)){
    start$small.site[which(start$x>area.grid[i,1] & start$x<(area.grid[i,1]+area.interval) &
                             start$y>area.grid[i,2] & start$y<(area.grid[i,2]+area.interval))] = 
      paste(start$site[which(start$x>area.grid[i,1] & start$x<(area.grid[i,1]+area.interval) &
                               start$y>area.grid[i,2] & start$y<(area.grid[i,2]+area.interval))],
            "smallsite", i)
  }
  
  if(any(is.na(start$small.site))){
    start = start[-which(is.na(start$small.site)),]
  }
  
  start$site = paste(start$site, start$small.site)
  #####

  ##########################
  # above is all about formatting the data. It should be kept separate.
  
  # Now, a more general function to compute metrics/ fit gams etc.
  
  hill_ef<-function(dataset
                    , site.col = "site"
                    , tax.col = "Gen_sp"
                    , fun.col = "carbon"
                    , ell.vec = seq(-5, 5, 0.5)){
  # set up abundance and function matrices
  #####
  # not sure we need these but here. 
  nsites = dataset %>% summarize(n_distinct({{site.col}}))
  nsp = dataset %>% summarize(n_distinct({{tax.col}}))
  
  # s.ab = matrix(0, nsites, nsp)
  # s.fn = matrix(0, nsites, nsp)
  site_d = dataset %>% 
    group_by(!!rlang::sym(site.col), !!rlang::sym(tax.col)) %>% 
    summarize(sp_ab = n()) %>% 
    ungroup() %>% 
    group_by({{site.col}}) %>% 
    summarize(dplyover::over(.x = ell.vec
                             , ~rarity(.data$sp_ab, l = .x)
                             , .names - "div_l_{x}")) 
      
  site_sums = dataset %>% 
    group_by({{site.col}}) %>% 
    summarize(ab = n()
              , ef = sum(!!rlang::sym(fun.col))
              ) %>% 
    mutate(pcf = ef/ab)
  
  site_D = left_join(site_sums, site_d )
  return(site_D)
  }

  
  map(ell.vec, function(ell)
  map(unique(dataset[ , site.col])){
    map_dfr()
    for(j in 1:nsp){
      s.ab[i,j] = sum(start$abun[which(start$site==unique(start$site)[i] &
                                         start$Gen_sp==unique(start$Gen_sp)[j])])
      s.fn[i,j] = sum(start$carbon[which(start$site==unique(start$site)[i] &
                                           start$Gen_sp==unique(start$Gen_sp)[j])])
    }
    print(paste("making matrices",i,"of",nsites))
  }
  
  c.ef = rowSums(s.fn)
  c.pcf = rowSums(s.fn)/rowSums(s.ab)
  
  rownames(s.ab) = unique(start$site)
  colnames(s.ab) = unique(start$Gen_sp)
  
  rownames(s.fn) = unique(start$site)
  colnames(s.fn) = unique(start$Gen_sp)
  #####
  
  # get the ef~div correlation
  #####

  div = matrix(NA, length(ell.seq), nsites)
  # for(z in 1:length(ell.seq)){
  #   
  #   q.loop = ell.seq[z]
  #   
  #   # get hill diversity
  #   if(ell.loop != 1){
  #     for(j in 1:nrow(s.ab)){
  #       q.hold = s.ab[j,]
  #       q.hold = q.hold[-which(q.hold==0)]
  #       div[z,j] = sum((q.hold/sum(q.hold))^ell.seq[z])^(1/(1-ell.seq[z]))
  #     }
  #   }
  #   
    
  map(1:length(ell.seq), function(ell){
    map(1:nrow(s.ab), function(x){
      div[ell, x] <<- rarity(s.ab[x, ], ell.seq[ell])
          })
    
    # if(ell.loop == 1){
    #   for(j in 1:nrow(s.ab)){
    #     q.hold = s.ab[j,]
    #     q.hold = q.hold[-which(q.hold==0)]
    #     div[z,j] = 2.71828^(-sum((q.hold/sum(q.hold))*log(q.hold/sum(q.hold))))
    #   }
    # }
    
    # models for diversity and function
    # the gam version takes longer
    ef.cor.out[are.int, ell]  <<- summary(gam(c.ef  ~ s(div[ell, ], k=3)))$r.sq
    pcf.cor.out[are.int, ell] <<- summary(gam(c.pcf ~ s(div[ell,], k=3)))$r.sq
    #ef.cor.out[q,z]  = cor(c.ef , div[z,], method="spearman")
    #pcf.cor.out[q,z] = cor(c.pcf, div[z,], method="spearman")
    
    # track progress
    print(paste("different ell values", ell,"of",length(ell.seq)))

  
  # save the diversity matrices, we will need one of them later
  div.out[[ell]] <<- div
  
  # save community-level EF, we will need it later
  c.ef.out[[ell]]  <<- c.ef
  c.pcf.out[[ell]] <<- c.pcf
  
  print(paste("different areas", are.int, "of",length(area.interval.vec)))
  })
}

# if you want to require that ef.cor.out and pcf.cor.out are not < 0
#ef.cor.out[  ef.cor.out<0] = 0
#pcf.cor.out[pcf.cor.out<0] = 0

# 2nd order question - what is the shape of the most predictive relationship
# save models for this

# first find the most predictive (both in terms of area and q)
pre.ef.best.z  = which( ef.cor.out==max( ef.cor.out), arr.ind = TRUE)
pre.pcf.best.z = which(pcf.cor.out==max(pcf.cor.out), arr.ind = TRUE)

# then get the diversity values associated with the above
ef.best.z  = div.out[[ pre.ef.best.z[1]]][ pre.ef.best.z[2],]
pcf.best.z = div.out[[pre.pcf.best.z[1]]][pre.pcf.best.z[2],]

ef.best.model   = gam(c.ef.out [[ pre.ef.best.z[1]]] ~ s( ef.best.z, k=3))
pcf.best.model  = gam(c.pcf.out[[pre.pcf.best.z[1]]] ~ s(pcf.best.z, k=3))
#####

# notes
#####
# as q-> inf, only the proportional abundance of the most abundant species matters
# div(q) approaches (1/(max(prop. abun)))
# so when a species increasingly dominates the community
# that denominator will increase, but that lowers div(q). so div(high q) measures evenness
# so if div(high q) has positive r, it means increases in evenness increase function
# and and negative r means increases in evenness decrease function

# as q-> neg inf, only the proportional abundance of the LEAST abundant species matters
# i don't have the exact math but roughly it means less even communities are more diverse
# i think the limit, for any site, is sum(abun)/min(abun)
# so if div(high -q) is most predictive with a positive r, it means low evenness increases
# function, and negative r means low evenness increases function
#####

# save output
#####
saveRDS( ef.cor.out, paste(site,  "ef.cor.GAM.k3.rds", sep=""))
saveRDS(pcf.cor.out, paste(site, "pcf.cor.GAM.k3.rds", sep=""))

saveRDS( ef.best.model, paste(site,  "ef.best.GAM.k3.rds", sep=""))
saveRDS(pcf.best.model, paste(site, "pcf.best.GAM.k3.rds", sep=""))
#####

# Yasuni (Ecuador) figures
#####
YAef  = readRDS("YAef.cor.GAM.k3.rds" ); YA25ef  = pmax(0, YA25ef )

YApcf = readRDS("YApcf.cor.GAM.k3.rds"); YA25pcf = pmax(0, YA25pcf)

ggplot()+
  
  geom_point(aes(x=ell.seq,y=YAef[1,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=ell.seq,y=YAef[2,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=ell.seq,y=YAef[3,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=ell.seq,y=YAef[4,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, Yasuni") +
  mytheme

ggplot()+
  
  geom_point(aes(x=ell.seq, y=YApcf[1,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=ell.seq, y=YApcf[2,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=ell.seq, y=YApcf[3,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=ell.seq, y=YApcf[4,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, Yasuni") +
  mytheme

# what does the most predictive relationship look like
YAef.best   = readRDS( "YAef.best.GAM.k3.rds")
YApcf.best  = readRDS("YApcf.best.GAM.k3.rds")

plot( YAef.best) # negative linear
plot(YApcf.best) # concave up
#####

# BCI Power Line Road (Panama) figures
#####
BCef  = readRDS("BCef.cor.GAM.k3.rds" ); BC25ef  = pmax(0, BC25ef )

BCpcf = readRDS("BCpcf.cor.GAM.k3.rds"); BC25pcf = pmax(0, BC25pcf)

ggplot()+

  geom_point(aes(x=q.seq,y=BCef[1,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=q.seq,y=BCef[2,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=q.seq,y=BCef[3,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=q.seq,y=BCef[4,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, BCI Road") +
  mytheme

ggplot()+
 
  geom_point(aes(x=q.seq, y=BCpcf[1,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=q.seq, y=BCpcf[2,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=q.seq, y=BCpcf[3,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=q.seq, y=BCpcf[4,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, BCI Road") +
  mytheme

# what does the most predictive relationship look like
BCef.best   = readRDS( "BCef.best.GAM.k3.rds")
BCpcf.best  = readRDS("BCpcf.best.GAM.k3.rds")

plot( BCef.best) # positive linear
plot(BCpcf.best) # concave down
#####

# Pasoh (Malaysia) figures
#####
PAef  = readRDS("PAef.cor.GAM.k3.rds" ); PA25ef  = pmax(0, PA25ef )

PApcf = readRDS("PApcf.cor.GAM.k3.rds"); PA25pcf = pmax(0, PA25pcf)

ggplot()+
  
  geom_point(aes(x=q.seq,y=PAef[1,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=q.seq,y=PAef[2,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=q.seq,y=PAef[3,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=q.seq,y=PAef[4,]), size=4, stroke=0.3, pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, Pasoh") +
  mytheme

ggplot()+
  
  geom_point(aes(x=q.seq, y=PApcf[1,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[1]) +
  geom_point(aes(x=q.seq, y=PApcf[2,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[2]) +
  geom_point(aes(x=q.seq, y=PApcf[3,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[3]) +
  geom_point(aes(x=q.seq, y=PApcf[4,]), size=4,stroke=0.3,pch=21, alpha=0.5, fill=scol[4]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, Pasoh") +
  mytheme

# what does the most predictive relationship look like
PAef.best   = readRDS( "PAef.best.GAM.k3.rds")
PApcf.best  = readRDS("PApcf.best.GAM.k3.rds")

plot( PAef.best) # positive linear
plot(PApcf.best) # concave up
#####

