.r <- function(){rm(list=ls());gc();source(".Rprofile")}

# Updated 29/9/23 ----
# Should by `final` set constructing thesis appendices.
# These are all under "stock" conditions


# Updated 02/01/23 ----
# No longer doing separate no-/survival plots since
# another investigation covers this in a more meaningful way 
# than purely eyeballing

# 1. Gaussian -------------------------------------------------------------
.r()
sim <- create.simfn()
X_ <- dataGen(sim)
(X <- Sample(X_, 5.5, T))
# (X2 <- Sample(X_, 5.5, T, include.survival = FALSE))
.plot.both(X, X2)

plot.staircase(X,  40, file.name = "gaussian.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-gaussian-nosurv.png", N = 1e5)

# Poisson -----------------------------------------------------------------
.r()
# Under "stock" conditions

sim <- create.simfn(list("poisson"))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.75, T))
# (X2 <- Sample(X_, 5.75, T, include.survival = FALSE))
# .plot.both(X, X2)

plot.staircase(X,  40, file.name = "poisson.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-poisson-nosurv.png", N = 1e5)

# Negative Binomial -------------------------------------------------------
.r()
# "Stock" conditions
sim <- create.simfn(list("negbin"))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.75, T))
# (X2 <- Sample(X_, 5.75, T, include.survival = FALSE))
# .plot.both(X, X2)

plot.staircase(X,  40, file.name = "negbin.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-negbin-nosurv.png", N = 1e5)

# Gamma -------------------------------------------------------------------
.r()
sim <- create.simfn(list("Gamma"))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.75, T))
# (X2 <- Sample(X_, 5.75, T, include.survival = FALSE))
# .plot.both(X, X2)

plot.staircase(X,  40, file.name = "Gamma.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-Gamma-nosurv.png", N = 1e5)

# Binomial ----------------------------------------------------------------
.r()
sim <- create.simfn(list("binomial"), arguments = list(D = diag(c(2,.5))))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.5, T, T))
# (X2 <- Sample(X_, 5.75, T, T, include.survival = FALSE))
# .plot.both(X, X2)

plot.staircase(X,  40, file.name = "binomial.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-binomial-intslope-nosurv.png", N = 1e5)

# Generalised Poisson -----------------------------------------------------
.r()         # Might have to run this `block` a few times.
sim <- create.simfn(list("genpois"), arguments = list(D = diag(c(0.30,0.05))))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.33, T, T)) # Might have to retry a few times 
# (X2 <- Sample(X_, 5.33, T, T, include.survival = FALSE)) # Might have to retry a few times 
# .plot.both(X, X2)

plot.staircase(X,  40, file.name = "genpois.png", N = 10000)
# plot.staircase(X2, 40, file.name = "default-genpois-intslope-nosurv.png", N = 1e5)
