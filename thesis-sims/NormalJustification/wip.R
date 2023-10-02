.r <- function(){rm(list=ls());source(".Rprofile")}

# Gaussian ----------------------------------------------------------------
.r()
sim <- create.simfn(arguments = list(n = 250, theta = c(-3, 0.6), gamma = 0.25))
(X_ <- dataGen(sim))
(X <- Sample(X_, 3, T))
(X2 <- Sample(X_, 5.25, T, include.survival = F))
.plot.both(X,X2)

psis1 <- prop.by.mi(X, 0.2, 0.3)
psis2 <- prop.by.mi(X2, 0.2, 0.3)
quick.psi.plot(psis1,psis2)

ga <- compare.theor.distn(X, X2)

plot(ga)
