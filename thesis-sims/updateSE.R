load("~/Downloads/fits5variate.RData")
load("/data/c0061461/GLMM_Paper_Sims/Revision1/simsets_5variate_Thu_Mar__9_17_06_09_2023.RData")

long.formulas <- fits[[1]]$ModelInfo$long.formulas
formulas <- lapply(long.formulas, parseFormula)
surv.formula <- fits[[1]]$ModelInfo$surv.formulas
family <- fits[[1]]$ModelInfo$family
beta.inds <- fits[[1]]$ModelInfo$inds$beta
b.inds <- fits[[1]]$ModelInfo$inds$b
K <- 5

GH <- statmod::gauss.quad.prob(3, 'normal')
w <- GH$w; v <- GH$n

redoSE <- function(ind){
  fit <- fits[[ind]]
  data <- sim.sets[[ind]]
  Omega <- fit$coeffs
  # Re-get data matrices
  dmats <- createDataMatrices(data, formulas)
  l0 <- fit$haz[,2]
  surv <- parseCoxph(surv.formula, data)
  sv <- surv.mod(surv, formulas, l0)
  # REs
  n <- fit$ModelInfo$n
  b <- ranef(fit)
  b <- lapply(1:n, function(i) b[i,])
  # Get SE
  II <- obs.emp.I2(Omega, dmats, surv, sv, b, sv$l0i, sv$l0u, 
                   w, v, n, family, 
                   K, sv$q, beta.inds, b.inds)
  setNames(sqrt(diag(solve(II$Hessian))),
           names(fit$SE))
}

args(obs.emp.I2)

a$SE
redoSE(1)
