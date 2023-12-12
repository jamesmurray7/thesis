#' #################################
#' EllipseCompare.
#' ----
#' Comparing the shape of ellipses produced
#' by normal approximations with/out the
#' survival density included

.r <- function(){rm(list=ls());gc();source(".Rprofile"); source("compareNapproxes.R")}
.xyz <- function(f) file.path(save.dir.file.path(family.dir.name(f)), paste0(f,"-ellipse-new.png"))

# Gaussian ----------------------------------------------------------------
.r()
sim <- create.simfn(list("gaussian"), arguments = list(n = 1000L, theta = c(-3,.1)))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("gaussian"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()

# Poisson -----------------------------------------------------------------
.r()
sim <- create.simfn(list("poisson"), arguments = list(n = 1000L, theta = c(-3,.1)))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("poisson"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()


# Binomial ----------------------------------------------------------------
.r()
sim <- create.simfn(list("binomial"), arguments = list(n = 500L, theta = c(-3,.1), beta = t(c(2.0,-0.1,0.1,0.2)),
                                                       D = matrix(c(.5, .125, .125, .09), 2, 2)))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("binomial"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()

# Gamma -------------------------------------------------------------------
.r()
sim <- create.simfn(list("Gamma"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("Gamma"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()


# Negbin ------------------------------------------------------------------
.r()
sim <- create.simfn(list("negbin"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("negbin"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()

# Genpois -----------------------------------------------------------------
.r()
sim <- create.simfn(list("genpois"), arguments = list(n = 500L, theta = c(-3,.1),
                                                      D = matrix(c(0.5, .125, .125, .05), 2, 2),
                                                      sigma = list(0.2), beta = t(c(.5,-.2,.05,.4))))
X_ <- dataGen(sim)
surv <- get.b(X_, T, T)
no.surv <- get.b(X_, F, T)
compare <- compare.N.approx(surv, no.surv)
png(filename = .xyz("genpois"), width = 140, height = 80, units = "mm", res = 500)
plot(compare)
dev.off()

