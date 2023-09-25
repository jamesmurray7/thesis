.r <- function(){rm(list=ls());source(".Rprofile")}

# 1. Gaussian -------------------------------------------------------------
.r()
# Under "stock" conditions
sim <- create.simfn()
X_ <- dataGen(sim)
X <- Sample(X_, 5.5, T)

plot.staircase(X, 40, file.name = "default.png")

# Inflated Covariance matrix D
sim2 <- create.simfn(arguments = list(D = diag(c(2.5, 0.9))))
X_2 <- dataGen(sim2)
(X2 <- Sample(X_2, 5.5, T))

plot.staircase(X, 40, file.name = "largeD.png")

# Poisson -----------------------------------------------------------------
.r()
# Under "stock" conditions

sim <- create.simfn(list("poisson"))
X_ <- dataGen(sim)
X <- Sample(X_, 5.75, T)

plot.staircase(X, 40, file.name = "default.png")


# Negative Binomial -------------------------------------------------------
.r()
# "Stock" conditions
sim <- create.simfn(list("negbin"))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.75, T))

plot.staircase(X, 40, file.name = "default.png")


# Gamma -------------------------------------------------------------------
.r()
sim <- create.simfn(list("Gamma"))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.75, T))

plot.staircase(X, 40, file.name = "default.png")

# Binomial ----------------------------------------------------------------
.r()
sim <- create.simfn(list("binomial"), arguments = list(D = diag(c(2,.5))))
X_ <- dataGen(sim)
(X <- Sample(X_, 5.5, T, T))

plot.staircase(X, 40, file.name = "default-intslope.png")

# Generalised Poisson -----------------------------------------------------
.r()         # Might have to run this `block` a few times.
sim <- create.simfn(list("genpois"), arguments = list(D = diag(c(0.30,0.05))))
X_ <- dataGen(sim)
X_
(X <- Sample(X_, 5.33, T, T))

plot.staircase(X, 40, file.name = "default-intslope.png")
