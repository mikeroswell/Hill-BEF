library(MeanRarity)

com1 = c(80, 50, 20, 5, 3, NA, NA)
com2 = c(20, 1, 1, 1, 1, 1, 1)

ell_vals = seq(-2, 5, by=0.05)


D1 = divpro(com1, -2, 5, 0.05)
D2 = divpro(com2, -2, 5, 0.05)

par(mfrow=c(2,1), mar=c(3,4,.5,.5))
barplot(rbind(com1, com2), beside=T, col=c("purple", "orange"), log="y",
        ylab="Abundance")

plot(D1$ell, D1$d, type="l", col="purple", lwd=3, ylim=c(1,25), log="y",
     xlab="ell value", ylab="Diversity")
points(D2$ell, D2$d, type="l", col="orange", lwd=3)
abline(v=c(-1,0,1), lty=3:1)
