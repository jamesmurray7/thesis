rm(list=ls()); gc()
source(".Rprofile")
sourceCpp("src/testMC.cpp")
load("/data/c0061461/THESIS/gmvjoint-MC/data_n500.RData")

N <- NROW(sim.sets)
MCfits <- vector("list", N)
MCfits <- lapply(MCfits, function(x) setNames(vector("list", 3), c("Ordinary", "Antithetic", "Sobol")))

pf <- parent.frame()
long.formulas <- list(
  as.formula("Y.1 ~ time + cont + bin + (1 + time|id)", pf),
  as.formula("Y.2 ~ time + cont + bin + (1 + time|id)", pf),
  as.formula("Y.3 ~ time + cont + bin + (1|id)", pf)
)

surv.formula <- Surv(survtime, status) ~ bin
environment(surv.formula) <- pf

disp.formulas <- replicate(3, as.formula("~1", pf), F)

fitfn <- function(data, MCtype){
  fit <- tryCatch(suppressWarnings(
    joint(long.formulas = long.formulas, surv.formula = surv.formula,
          data = data, family = list("gaussian", "poisson", "binomial"),
          disp.formulas = disp.formulas, MCtype = MCtype, N = 1e2,
          sf = matrix(1,5,5),
          control = list(return.dmats = FALSE))
  ), error = function(e) NULL)
  fit
}

pb <- utils::txtProgressBar(max = N, style = 3)
for(j in 1:N){
  # MCfits[[j]]$Ordinary <- fitfn(sim.sets[[j]], "ordinary")
  MCfits[[j]]$Antithetic <- fitfn(sim.sets[[j]], "antithetic")
  MCfits[[j]]$Sobol <- fitfn(sim.sets[[j]], "sobol")
  utils::setTxtProgressBar(pb, j)
}
save(MCfits, file = "/data/c0061461/THESIS/gmvjoint-MC/fitsMC4-nosf-n500.RData")
