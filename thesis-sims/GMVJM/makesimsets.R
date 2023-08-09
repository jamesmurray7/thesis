rm(list=ls())
source(".Rprofile")
N <- 300

# "basic" trivariate ------------------------------------------------------
out <- save.dir.file.path("data")
out.file <- save.dir.file.path("Triv.RData", out)
fn <- create.simfn(arguments = list(n=500, random.formulas = list(~time, ~time, ~1)))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)


# Univariate Gamma --------------------------------------------------------
out.file <- save.dir.file.path("UnivGa.RData", out)
fn <- create.simfn(family = list("Gamma"),
                   arguments = list(n=500))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)

# Univariate Negbin -------------------------------------------------------
out.file <- save.dir.file.path("UnivNB.RData", out)
fn <- create.simfn(family = list("negbin"),
                   arguments = list(n=500))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)

# Univariate GP1 ----------------------------------------------------------
out.file <- save.dir.file.path("UnivGa.RData", out)
fn <- create.simfn(family = list("genpois"), # This works nice enough...
                   arguments = list(n=500, sigma = list(-0.2), beta = c(0, .5, -.1, .2), D = diag(c(.2,.05))))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)
