# Tina sim without shiny
library(tidyverse)
library(furrr)
library(tictoc)




# helpers
Ln = function(x) {ifelse(x != 0, log(x), 0)}
Exp = function(x, y) {ifelse(x != 0, x^y, 0)}
get.Dq = function(vec, q) {
  if(sum(vec) == 0) {return(0)}
  ifelse(q != 1,
         sum(Exp(vec/sum(vec), q))^(1/(1 - q)),
         exp(-sum((vec/sum(vec))*Ln(vec/sum(vec))))) }

# # Tina's code to pull out "families" of curves
# 
# 
# ell.cor %>% summarise(ell.pos.cor = ell[which.max(Df.cor)],
#                       ell.neg.cor = ell[which.min(Df.cor)],
#                       ell.min.cor = ell[which.min(abs(Df.cor))])
# 
# studies = add_column(studies, group = do.call(cbind, Splits) %>% apply(1, which))
# 
# Splits = with(studies, 
#               list(ell.pos.cor >= 2 & ell.neg.cor <=-2,
#                    ell.pos.cor >= 2 & ell.neg.cor > -2,
#                    ell.pos.cor < ell.neg.cor,
#                    ell.pos.cor > -2 & ell.pos.cor < 2 & 
#                      ell.neg.cor > -2 & ell.neg.cor < 2 &
#                      ell.pos.cor > ell.neg.cor,
#                    ell.pos.cor > -2 & ell.pos.cor < 2 & 
#                      ell.neg.cor <= -2))
# 
# do.call(cbind, Splits) %>% rowSums


# Fixed parameters
sites <- 20
species_number <- 200

# variable ones
abundance_mean_v <- floor(10^seq(1,2.8, 0.4)) #150
abundance_dispersion_v <- seq(0.05, 2.15, 0.35) #1.5
function_slope_v <- seq(-2, 2, 0.5) # -1.1
function_error_v <- 10^seq(-1, 2, 0.5) #50


# abundance_mean  <- 150
# abundance_dispersion <- 1.5
# function_slope <- -1.1
# function_error <- 50

abundance_mean <- 10
abundance_dispersion <- 0.4
function_slope <- -0.5
function_error <- 0.31

tic()
pdf("figures/try_param_combos_for_peak.pdf")
nc <- 7
plan(strategy = multiprocess, workers = nc)
future_map(abundance_mean_v, function(abundance_mean){
  map(abundance_dispersion_v, function(abundance_dispersion){
  map(function_slope_v, function(function_slope){
  map(function_error_v, function(function_error){

   

        
# set site SAD params
site.abund_mu <- rgamma(sites, shape = abundance_mean)
site.abund_size = rgamma(sites, shape = abundance_dispersion)
  
a <- sapply(1:as.numeric(species_number), function(i) {
    sapply(1:as.numeric(sites), function(k) {
      rnbinom(1,  
              mu = site.abund_mu[k], 
              size = site.abund_size[k])})
  }) %>% t

a[which.max(rowSums(a)), which(colSums(a) == 0)] <- 1

ab_df<-data.frame(t(a)) %>% 
  rownames_to_column("site") %>% 
  pivot_longer(cols = 2:species_number+1, names_to= "sp", values_to = "Abundance")

ab_df<-ab_df %>% 
  mutate(Biomass = rgamma(1, rate = function_error/exp(5+function_slope*Ln(Abundance)), shape = function_error))

# 
amu <- apply(a, 2, mean)
avar <-apply(a, 2, sd)
rich <- apply(a, 2, function(x)sum(x >0))
#


z <- sapply(1:as.numeric(sites), function(i) {
  y_true = exp(5 + function_slope*Ln(a[,i])) # exp link glm
  shape = function_error # larger shape values = lower error
  rgamma(species_number, rate = shape/y_true, shape = shape)
}) # close sapply

# plot(a, z)
# hist(rich)
# plot(amu, avar)

q <- seq(-3, 5, by = 0.2)
# derived function and diversity
f <- colSums(a*z)
Dq <- sapply(q, function(q) apply(a, 2, get.Dq, q))
# sites in rows, q' in columns
# x.ax.val = seq(min(Dq), max(Dq), by = 1)


  # plot(a, z, # col=alpha("gray", 0.2),
  #      main = "site selection function")
  # for(i in 1:as.numeric(sites)) {
  #   X = a[,i]
  #   #this is fiddly
  #   mod = glm(z[,i] ~ X, family=Gamma(link="log"))
  #   x = seq(1, max(a, by=1))
  #   y = predict(mod, newdata=data.frame(X=x), type="response")
  #   points(x, y, type="l")
  # }
  # 
  # 
  # 
  # plot(range(q), range(Dq), type="n",
  #      xlab = "q", ylab = "Dq",
  #      main = "site diversity profiles")
  # apply(Dq, 1, function(x) points(q, x,
  #                                   type="l", lwd=2))
  # abline(v=0, lty=2)
  # evenness = log(Dq[,q==2])/log(Dq[,q==0])
  # mtext(paste("mean even =", round(mean(evenness), 2)), side=3, line=-2)
  # mtext(paste("sd even =", round(sd(evenness), 2)), side=3, line=-3)
  # 
  # 
  # 
  # plot(f ~ Dq[,1], ylim = range(f), xlim=range(Dq),
  #      type = "n", xlab = "Dq",
  #      main = "q-BEFs")
  # Col=rainbow(ncol(Dq))
  # sapply(1:length(q), function(x) points(Dq[, x], f, col=Col[x]))
  # ### linear mod
  # # plot(f ~ Dq[,1], ylim = range(f), xlim=range(Dq), type = "n")
  # sapply(1:length(q), function(x) abline(lm(f ~ Dq[, x]), col=Col[x], lwd=2))
  BEF <-  tryCatch(apply(Dq, 2, function(x) { cor(log(f), log(x)) })
                   , error = function(e) NA
                   , warning = function(w) NA )

  
  

  befDat<-tryCatch(data.frame(ell =1-q, BEF, rb = as.factor(1:ncol(Dq))), error = function(e){e})  
  if(all(is.na(BEF)) |!is.data.frame(befDat)){
    print(data.frame(abMu = abundance_mean
                     , abDisp = abundance_dispersion
                     , abEF_slope = function_slope
                     , abEF_err = function_error))
    return(NULL)}
  if(abs(range(BEF)[[2]]-range(BEF)[[1]])<0.4 & 
    ( !(13<which.max(BEF) & which.max(BEF)<26) | !(13<which.min(BEF) & which.min(BEF)<26))){ 
    
    print(data.frame(abMu = abundance_mean
                , abDisp = abundance_dispersion
                , abEF_slope = function_slope
                , abEF_err = function_error))
    return(NULL)}
  else{ 
   return(
    befDat %>% ggplot( aes(ell, BEF, color = rb) )+
      geom_point()+
      labs(y = "BEF", title = paste("abMu =", abundance_mean
                                    , "abDisp =", abundance_dispersion
                                    , "\nabEF_slope =", function_slope
                                    , "abEF_err =", function_error
                                    , "\nell-BEF raw cor (log-log)")) +
      scale_color_manual(values = rainbow(ncol(Dq))) +
      geom_vline(xintercept = 1, linetype = 2 ) +
      theme_classic()+ 
      theme(legend.position="none")
    )}

    
  
  


      })
    })
  })
})
dev.off()
print(toc())

# pdf("figures/test_looping.pdf")
# future_map(1:100, function(neplo){
#   if(neplo>0){
#     par(mar = c(5,4,4,5)+.1)
#   print(BEF_plot(q, BEF))}
# })
# dev.off()

