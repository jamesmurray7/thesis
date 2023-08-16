rm(list=ls())
source(".Rprofile")
N <- 300

# "basic" trivariate ------------------------------------------------------
out <- save.dir.file.path("data")
out.file <- save.dir.file.path("Triv.RData", out)
fn <- create.simfn(arguments = list(n=500, random.formulas = list(~time, ~time, ~1),
                                    theta = c(-3, .1)))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)


# Univariate Gamma --------------------------------------------------------
out.file <- save.dir.file.path("UnivGa.RData", out)
fn <- create.simfn(family = list("Gamma"),
                   arguments = list(n=500, theta = c(-2.8, .1)))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)

# Univariate Negbin -------------------------------------------------------
out.file <- save.dir.file.path("UnivNB.RData", out)
fn <- create.simfn(family = list("negbin"),
                   arguments = list(n=500, theta = c(-2.8, .1)))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)

# Univariate GP1 ----------------------------------------------------------
out.file <- save.dir.file.path("UnivGP.RData", out)
fn <- create.simfn(family = list("genpois"),
                   arguments = list(n=500, theta = c(-2.8, .1), 
                                    beta = c(1,.05,-.05,.1),
                                    sigma = list(-0.3),
                                    random.formulas = list(~1), D = matrix(.30,1,1)))
sim.sets <- createNsims(fn, N)
DataSummary(sim.sets)
save(sim.sets, file = out.file)
