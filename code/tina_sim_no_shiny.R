# Tina sim without shiny

# params
sites <- 20
species_number <- 250
abundance_mean <- 12
abundance_dispersion <- .1
function_slope <- 5
function_error <- 0.3



# helpers
Ln = function(x) {ifelse(x != 0, log(x), 0)}
Exp = function(x, y) {ifelse(x != 0, x^y, 0)}
get.Dq = function(vec, q) {
  if(sum(vec) == 0) {return(0)}
  ifelse(q != 1,
         sum(Exp(vec/sum(vec), q))^(1/(1 - q)),
         exp(-sum((vec/sum(vec))*Ln(vec/sum(vec))))) }


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

ab_df<-data.frame(t(a)) %>% add_rownames("site") %>% pivot_longer(cols = 2:species_number+1, names_to= "sp", values_to = "Abundance")
ab_df<-ab_df %>% mutate(Biomass = rgamma(1, rate = function_error/exp(5+function_slope*Ln(Abundance)), shape = function_error))
# 
# amu<- apply(a, 2, mean)
# avar<-apply(a, 2, sd)
# rich<- apply(a, 2, function(x)sum(x >0))
# 
# hist(rich)
# plot(amu, avar)
# z <- sapply(1:as.numeric(sites), function(i) {
#   y_true = exp(5 + function_slope*Ln(a[,i])) # exp link glm
#   shape = function_error # larger shape values = lower error
#   rgamma(species_number, rate = shape/y_true, shape = shape)
# }) # close sapply
# 
# plot(a, z)
# 
# q <- seq(-2, 4, by = 0.2)
# # derived function and diversity
# f <- colSums(a*z)
# Dq <- sapply(q, function(q) apply(a, 2, get.Dq, q)) 
# # sites in rows, q' in columns
# # x.ax.val = seq(min(Dq), max(Dq), by = 1)
# 
# a
# 
# 
# 
# 
#   plot(a, z, # col=alpha("gray", 0.2),
#        main = "site selection function")
#   for(i in 1:as.numeric(sites)) {
#     X = a[,i]
#     mod = glm(z[,i] ~ X, family=Gamma(link="log"))
#     x = seq(1, max(a, by=1))
#     y = predict(mod, newdata=data.frame(X=x), type="response")
#     points(x, y, type="l")
#   } 
# 
# 
# output$plotgraph2 = renderPlot({
#   plot(range(q), range(Dq), type="n",
#        xlab = "q", ylab = "Dq",
#        main = "site diversity profiles")
#   apply(Dq, 1, function(x) points(q, x, 
#                                     type="l", lwd=2))
#   abline(v=0, lty=2)
#   evenness = log(Dq[,q==2])/log(Dq[,q==0])
#   mtext(paste("mean even =", round(mean(evenness), 2)), side=3, line=-2)
#   mtext(paste("sd even =", round(sd(evenness), 2)), side=3, line=-3)
#   
# })
# 
# output$plotgraph3 = renderPlot({
#   plot(f ~ Dq[,1], ylim = range(f), xlim=range(Dq), 
#        type = "n", xlab = "Dq",
#        main = "q-BEFs")
#   Col=rainbow(ncol(Dq))
#   sapply(1:length(q), function(x) points(Dq[, x], f, col=Col[x]))
#   ### linear mod
#   # plot(f ~ Dq[,1], ylim = range(f), xlim=range(Dq), type = "n")
#   sapply(1:length(q), function(x) abline(lm(f ~ Dq[, x]), col=Col[x], lwd=2))
# }) 
# 
# output$plotgraph4 = renderPlot({
#   plot(q, apply(Dq, 2, function(x) {
#     cor(f, x)
#   }), ylab = "f", col=rainbow(ncol(Dq())), pch=19, cex=2,
#   main = "q-BEF raw cor")
#   abline(v=0, lty=2)
#   
#   
# })
# } # close server
