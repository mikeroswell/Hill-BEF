# Tina's partition code
library(MeanRarity)
library(tidyverse)
# format the bee data
a.list = readRDS("data/a.list.rds")
z.list = readRDS("data/z.list.rds")
a.data <- a.list %>% lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% as.matrix}) 
z.t.data <- z.list %>% lapply(function(x) {rownames(x) = x$X; x[,-1] %>% 
  as.data.frame.matrix %>% as.matrix}) 
z.data = sapply(1:length(z.t.data),
                function(se) {z = z.t.data[[se]]/a.data[[se]]
                z[is.nan(z)] = 0; z})


# analysis workflow

n_data <-3
map(1:n_data, function(se){
  ell_vals = c(-10, -2, 1, 2, 10) 
# matrices
  a = a.data[[se]] # select one of three focal plants
  z = z.data[[se]]
# vector summaries
  p = apply(a, 2, function(a) a/sum(a)) 
  t = colSums(a*z)
  A = colSums(a)
  zbar = colSums(p*z)
  # diversities
  D = sapply(ell_vals, function(ell) {
    apply(a, 2, rarity, l = ell) })
  # fit models, get coefficients
  # model predictions (mp)
  mp <- lapply(1:5, function(ell) { 
    x = D[, ell] # ell=3
    Amod = lm(log(A) ~ log(x)) 
    zmod = lm(log(zbar) ~ log(x))
    # each model has an intercept and slope (coef 1 and coef 2)
    cbind( A.predict = exp(coef(Amod)[1])*x^(coef(Amod)[2]) %>% as.numeric
         , zbar.predict = exp(coef(zmod)[1])*x^(coef(zmod)[2]) %>% as.numeric
         # derived
         , t.predict = (exp(coef(Amod)[1] + 
                              coef(zmod)[1]) * x^(coef(Amod)[2] +
                              coef(zmod)[2])) %>% 
           as.numeric
        
         # nulls
         , t.A.predict = (exp(coef(Amod)[1] + 
                                coef(zmod)[1]) * x^(coef(Amod)[2] + 0)) %>% 
           as.numeric
         , t.zmod.predict = (exp(coef(Amod)[1] + coef(zmod)[1]) * 
                               x^(0 + coef(zmod)[2])) %>% 
           as.numeric
         ) %>% 
      .[order(x), ] # CAREFUL reorder by x
         }
  ) 
  # get ylim range for total function models
  # want to keep yaxis the same across ell values
  Ylim <- sapply(mp, function(p)
    range(p[,c("t.predict", "t.A.predict", "t.zmod.predict")])) %>%
    cbind(c(0,0)) %>% range 
  # plot
  par(mfrow=c(5,3), mar=c(2,2,1,1)) 
  for(ell in 1:length(ell_vals)){ # in this ell loop the ells are just indices
    x = D[, ell]
    y = mp[[ell]]
    # plot abundance, per-capita function 
    plot(x, A, col="red") 
      points(x[order(x)], y[, "A.predict"]
               , col="red", type="l", lwd=2) 
    plot(x, zbar, col="blue")
      points(x[order(x)], y[, "zbar.predict"], col="blue", type="l", lwd=2)
      # plot total function and factors
      plot(x, t, ylim=Ylim)
      points(x[order(x)], y[, "t.predict"], type="l", lwd=2) 
      points(x[order(x)], y[, "t.A.predict"], type="l", col="red", lwd=2) 
      points(x[order(x)], y[, "t.zmod.predict"], type="l", col="blue", lwd=2)
  }
  
  
 # basically the same but now for the correlations
  ell_vals = seq(-10, 10, by=1)
  D2 = sapply(ell_vals, function(ell) {
    apply(a, 2, rarity, l = ell) })
  mb<- sapply(1:length(ell_vals), function(ell) { 
    x = D2[, ell] # ell=3
    Amod2 = lm(log(A) ~ log(x))
    zmod2 = lm(log(zbar) ~ log(x))
    beta.t = (coef(Amod2)[2] + coef(zmod2)[2]) %>% as.numeric
    beta.A = coef(Amod2)[2] %>% as.numeric
    beta.z = coef(zmod2)[2] %>% as.numeric
    t.int = as.numeric(confint(lm(log(t) ~ log(x)), 2)) 
    A.int = confint(Amod2, 2)
    z.int = confint(zmod2, 2)
    # tmod = lm(log(t) ~ log(x)) # double-checking
    # coef(Amod)[2] + coef(zmod)[2]; coef(tmod)[2]
    # cor(log(t), log(x)); (beta.A+beta.z)*sd(log(x))/sd(log(t)) 
    c(
      beta.t
      , beta.A
      , beta.z
      , t.lower = t.int[1]
      , t.upper = t.int[2]
      , A.lower = A.int[1]
      , A.upper = A.int[2]
      , z.lower = z.int[1]
      , z.upper = z.int[2]
      , t.cor = cor(log(t)
                    , log(x))
      , A.cor = beta.A * sd(log(x))/sd(log(t))
      , z.cor = beta.z * sd(log(x))/sd(log(t))
    )
    }) # close sapply; mp for "model predictions"
    # plot
    Col=c("black", "red", "blue")
    par(mfrow=c(1,2)) # slope plot
    plot(range(ell_vals), range(mb[1:9,]), type="n",
         xlab = "ell value", ylab = "BEF slope")
      abline(v=1, lty=3)
      abline(h=0, lty=3)
      polygon(c(rev(ell_vals), ell_vals)
              , c(rev(mb["t.upper",]), (mb["t.lower", ]))
              , col=t_col("gray", 30), border=NA)
      polygon(c(rev(ell_vals), ell_vals)
              , c(rev(mb["A.upper", ]), mb["A.lower",])
              , col=t_col("pink", 30), border=NA)
      polygon(c(rev(ell_vals), ell_vals)
              , c(rev(mb["z.upper", ]), mb["z.lower",])
              , col=t_col("lightblue", 30), border=NA)
    for(i in 1:3) {points(ell_vals, mb[i,], type="l", col=Col[i], lwd=2)}
    # correlation plot
    plot(range(ell_vals), range(mb[10:12,])
         , type="n", xlab="ell value", ylab = "BEF correlation")
      abline(v=1, lty=3); abline(h=0, lty=3)
      points(ell_vals, mb["t.cor", ], col="black", type="l", lwd=2) 
      points(ell_vals, mb["A.cor", ], col="red", type="l", lwd=2) 
      points(ell_vals, mb["z.cor", ], col="blue", type="l", lwd=2)

})

# helper function: 
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk
t_col <- function(color, percent = 50, name = NULL) {
  # # #
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255,
               alpha = (100 - percent) * 255 / 100, names = name)
  ## Save the color
  invisible(t.col) }
## END
