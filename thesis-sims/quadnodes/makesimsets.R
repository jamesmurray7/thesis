rm(list=ls())
source(".Rprofile")
N <- 200 # 200 might be overkill but can just shed down to first 100!

# Trivariate sim ----------------------------------------------------------
# (using randomfail.sim() from datagen1.R)
# which creates a set of trivariate joint data with failure rate between
# ~20ish-50%ish...

sim.sets <- replicate(N, randomfail.sim(), simplify = F)
save(sim.sets, file = save.dir.file.path("data/data.RData", save.dir))