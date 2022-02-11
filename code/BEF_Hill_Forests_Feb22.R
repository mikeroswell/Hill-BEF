# choose site and area subdivision method
#####
site          = "YA"                               # c("PA", "BC", "YA")
area.interval = 100                                # c(20, 25, 33, 50, 100)
q.seq         = seq(from= -5.0, to= 5.0, by=0.20)  # range of q (Hill numbers)

# to get the figures to work,
# you need to run all the listed sites and area.intervals separately (up to line 149)
# how does area interval work? e.g. area interval 20 gives a 5x5 grid (100/20 = 5)
#####

# read and process data, set up for area subdivision loop
#####
start = read.csv(paste("/Users/Mark/Desktop/OldDesktop/GlobalDom_final/", site, "_biomass.csv",
                 sep=""))

if(site=="PA"){start = start[which(start$Year==2012),]}
start$abun = 1
colnames(start)[which(colnames(start)=="X1ha.Plot.Number")] = "site"

colnames(start)[which(colnames(start)=="X1ha.Plot.X.Coordinate")] = "x"
colnames(start)[which(colnames(start)=="X1ha.Plot.Y.Coordinate")] = "y"

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


# set up abundance and function matrices
#####

nsites = length(unique(start$site))
nsp = length(unique(start$Gen_sp))

s.ab = matrix(0, nsites, nsp)
s.fn = matrix(0, nsites, nsp)

for(i in 1:nsites){
  for(j in 1:nsp){
    s.ab[i,j] = sum(start$abun[which(start$site==unique(start$site)[i] &
                                       start$Gen_sp==unique(start$Gen_sp)[j])])
    s.fn[i,j] = sum(start$carbon[which(start$site==unique(start$site)[i] &
                                         start$Gen_sp==unique(start$Gen_sp)[j])])
  }
  print(i)
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
ef.cor.out = c()
pcf.cor.out = c()

div.z = vector("list", length(q.seq))

div = matrix(NA, length(q.seq), nsites)
for(z in 1:length(q.seq)){
  
  q.loop = q.seq[z]
  
  # get hill diversity
  if(q.loop != 1){
    for(j in 1:nrow(s.ab)){
      q.hold = s.ab[j,]
      q.hold = q.hold[-which(q.hold==0)]
      div[z,j] = sum((q.hold/sum(q.hold))^q.seq[z])^(1/(1-q.seq[z]))
    }
  }
  
  if(q.loop == 1){
    for(j in 1:nrow(s.ab)){
      q.hold = s.ab[j,]
      q.hold = q.hold[-which(q.hold==0)]
      div[z,j] = 2.71828^(-sum((q.hold/sum(q.hold))*log(q.hold/sum(q.hold))))
    }
  }
  
  # models for diversity and function
  # the gam version takes longer
  ef.cor.out[z]  = summary(gam(c.ef  ~ s(div[z,], k=3)))$r.sq
  pcf.cor.out[z] = summary(gam(c.pcf ~ s(div[z,], k=3)))$r.sq
  #ef.cor.out[z]  = cor(c.ef , div[z,], method="spearman")
  #pcf.cor.out[z] = cor(c.pcf, div[z,], method="spearman")
  
  # track progress
  print(z)
}

# on hold for now
# 2nd order question - what is the shape of the most predictive relationship
# save models for this
#ef.best.z = which(ef.cor.out==max(ef.cor.out))
#ef.best.z = which(ef.cor.out==max(ef.cor.out))

#ef.best.model  = gam(c.ef ~ s(div[ef.best.z, ], k=3))
#pcf.best.model = gam(c.ef ~ s(div[pcf.best.z,], k=3))
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
saveRDS(ef.cor.out, paste(site, area.interval, "ef.cor.GAM.k3.rds", sep=""))
saveRDS(pcf.cor.out, paste(site, area.interval, "pcf.cor.GAM.k3.rds", sep=""))
#####

# Yasuni (Ecuador) figures
#####
YA20ef   = readRDS("YA20ef.cor.GAM.k3.rds"  );   YA20ef = pmax(0, YA20ef  )
YA25ef   = readRDS("YA25ef.cor.GAM.k3.rds"  );   YA25ef = pmax(0, YA25ef  )
YA33ef   = readRDS("YA33ef.cor.GAM.k3.rds"  );   YA33ef = pmax(0, YA33ef  )
YA50ef   = readRDS("YA50ef.cor.GAM.k3.rds"  );   YA50ef = pmax(0, YA50ef  )
YA100ef  = readRDS("YA100ef.cor.GAM.k3.rds" );  YA100ef = pmax(0, YA100ef )

YA20pcf  = readRDS("YA20pcf.cor.GAM.k3.rds" ); YA20pcf  = pmax(0, YA20pcf )
YA25pcf  = readRDS("YA25pcf.cor.GAM.k3.rds" ); YA25pcf  = pmax(0, YA25pcf )
YA33pcf  = readRDS("YA33pcf.cor.GAM.k3.rds" ); YA33pcf  = pmax(0, YA33pcf )
YA50pcf  = readRDS("YA50pcf.cor.GAM.k3.rds" ); YA50pcf  = pmax(0, YA50pcf )
YA100pcf = readRDS("YA100pcf.cor.GAM.k3.rds"); YA100pcf = pmax(0, YA100pcf)

# on hold for now
#YA.combined = rbind(YA20ef, YA25ef, YA33ef, YA50ef, YA100ef)
#YA.ef.best = which(YA.combined == max(YA.combined), arr.ind = TRUE)

scol = viridis(5)

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=YA20ef), fill=scol[1],  col=scol[1],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA25ef), fill=scol[2],  col=scol[2],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA33ef), fill=scol[3],  col=scol[3],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA50ef), fill=scol[4],  col=scol[4],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA100ef), fill=scol[5], col=scol[5],  alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=YA20ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=YA25ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=YA33ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=YA50ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=YA100ef),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, Ecuador") +
  mytheme

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=YA20pcf), fill=scol[1], col=scol[1], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA25pcf), fill=scol[2], col=scol[2], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA33pcf), fill=scol[3], col=scol[3], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA50pcf), fill=scol[4], col=scol[4], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=YA100pcf),fill=scol[5], col=scol[5], alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=YA20pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=YA25pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=YA33pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=YA50pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=YA100pcf),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, Ecuador") +
  mytheme

# needs more specific setup to work
#draw(gam(c.ef   ~ s(div[34,], k=3))) # negative linear
#draw(gam(c.pcf  ~ s(div[50,], k=3))) # concave up

#####

# BCI Power Line Road (Panama) figures
#####
BC20ef   = readRDS("BC20ef.cor.GAM.k3.rds"  );   BC20ef = pmax(0, BC20ef  )
BC25ef   = readRDS("BC25ef.cor.GAM.k3.rds"  );   BC25ef = pmax(0, BC25ef  )
BC33ef   = readRDS("BC33ef.cor.GAM.k3.rds"  );   BC33ef = pmax(0, BC33ef  )
BC50ef   = readRDS("BC50ef.cor.GAM.k3.rds"  );   BC50ef = pmax(0, BC50ef  )
BC100ef  = readRDS("BC100ef.cor.GAM.k3.rds" );  BC100ef = pmax(0, BC100ef )

BC20pcf  = readRDS("BC20pcf.cor.GAM.k3.rds" ); BC20pcf  = pmax(0, BC20pcf )
BC25pcf  = readRDS("BC25pcf.cor.GAM.k3.rds" ); BC25pcf  = pmax(0, BC25pcf )
BC33pcf  = readRDS("BC33pcf.cor.GAM.k3.rds" ); BC33pcf  = pmax(0, BC33pcf )
BC50pcf  = readRDS("BC50pcf.cor.GAM.k3.rds" ); BC50pcf  = pmax(0, BC50pcf )
BC100pcf = readRDS("BC100pcf.cor.GAM.k3.rds"); BC100pcf = pmax(0, BC100pcf)

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=BC20ef), fill=scol[1],  col=scol[1],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC25ef), fill=scol[2],  col=scol[2],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC33ef), fill=scol[3],  col=scol[3],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC50ef), fill=scol[4],  col=scol[4],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC100ef), fill=scol[5], col=scol[5],  alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=BC20ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=BC25ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=BC33ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=BC50ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=BC100ef),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, Panama") +
  mytheme

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=BC20pcf), fill=scol[1], col=scol[1], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC25pcf), fill=scol[2], col=scol[2], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC33pcf), fill=scol[3], col=scol[3], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC50pcf), fill=scol[4], col=scol[4], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=BC100pcf),fill=scol[5], col=scol[5], alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=BC20pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=BC25pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=BC33pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=BC50pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=BC100pcf),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, Panama") +
  mytheme

# needs more specific setup to work
#draw(gam(c.ef   ~ s(div[22,], k=3))) # positive linear
#draw(gam(c.pcf  ~ s(div[45,], k=3))) # concave down

#####

# Pasoh (Malaysia) figures
#####
PA20ef   = readRDS("PA20ef.cor.GAM.k3.rds"  );   PA20ef = pmax(0, PA20ef  )
PA25ef   = readRDS("PA25ef.cor.GAM.k3.rds"  );   PA25ef = pmax(0, PA25ef  )
PA33ef   = readRDS("PA33ef.cor.GAM.k3.rds"  );   PA33ef = pmax(0, PA33ef  )
PA50ef   = readRDS("PA50ef.cor.GAM.k3.rds"  );   PA50ef = pmax(0, PA50ef  )
PA100ef  = readRDS("PA100ef.cor.GAM.k3.rds" );  PA100ef = pmax(0, PA100ef )

PA20pcf  = readRDS("PA20pcf.cor.GAM.k3.rds" ); PA20pcf  = pmax(0, PA20pcf )
PA25pcf  = readRDS("PA25pcf.cor.GAM.k3.rds" ); PA25pcf  = pmax(0, PA25pcf )
PA33pcf  = readRDS("PA33pcf.cor.GAM.k3.rds" ); PA33pcf  = pmax(0, PA33pcf )
PA50pcf  = readRDS("PA50pcf.cor.GAM.k3.rds" ); PA50pcf  = pmax(0, PA50pcf )
PA100pcf = readRDS("PA100pcf.cor.GAM.k3.rds"); PA100pcf = pmax(0, PA100pcf)

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=PA20ef), fill=scol[1],  col=scol[1],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA25ef), fill=scol[2],  col=scol[2],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA33ef), fill=scol[3],  col=scol[3],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA50ef), fill=scol[4],  col=scol[4],  alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA100ef), fill=scol[5], col=scol[5],  alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=PA20ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=PA25ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=PA33ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=PA50ef), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=PA100ef),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("EF~d, Malaysia") +
  mytheme

ggplot()+
  
  #stat_smooth(aes(x=q.seq, y=PA20pcf), fill=scol[1], col=scol[1], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA25pcf), fill=scol[2], col=scol[2], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA33pcf), fill=scol[3], col=scol[3], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA50pcf), fill=scol[4], col=scol[4], alpha=0.25, span=0.1, size=2) +
  #stat_smooth(aes(x=q.seq, y=PA100pcf),fill=scol[5], col=scol[5], alpha=0.25, span=0.1, size=2) +
  
  geom_point(aes(x=q.seq,y=PA20pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[1]) +
  geom_point(aes(x=q.seq,y=PA25pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[2]) +
  geom_point(aes(x=q.seq,y=PA33pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[3]) +
  geom_point(aes(x=q.seq,y=PA50pcf), size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[4]) +
  geom_point(aes(x=q.seq,y=PA100pcf),size=4,stroke=0.3,pch=21,col="black",alpha=0.5,fill=scol[5]) +
  
  geom_hline(aes(yintercept=0)) +
  
  scale_x_continuous("q") +
  scale_y_continuous("PCF~d, Malaysia") +
  mytheme

# needs more specific setup to work
#draw(gam(c.ef  ~ s(div[33,], k=3))) # positive linear
#draw(gam(c.pcf ~ s(div[51,], k=3))) # concave up
#####

