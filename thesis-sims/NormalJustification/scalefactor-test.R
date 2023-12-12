#' ##############################################
#' Dummy program to try and find a scale-factor #
#' to place on Sigma.hat to reduce psi_i        #  
#' ##############################################
wrap <- function(f = "gaussian"){
  obj <- function(a, sim){
    df <- do.call(c, lapply(1:NROW(sim$Walks), function(i){
      this.df <- sim$df[sim$df$id==i,]
      this.b.hat <- c(sim$b.hats[[i]])
      psi <- check.walks(sim$Walks[[i]], this.b.hat, sim$Sigmas[[i]] * a)
      psi
    }))
    abs(mean(df - 0.95))
  }
  
  sim <- create.simfn(list(f), arguments = list(n = 100L, theta = c(-3,.1)))
  X_ <- dataGen(sim)
  S <- Sample(X_, return.walks = T, TUNE = 3, force.intslope = T, b.dist = "t", df = 4,burnin = 1000, NMC = 5000L)
  
  optim(1, obj, NULL, sim = S, method = 'Brent', lower = 1e-2, upper = 2,
        control = list(abstol = 1e-3, reltol = 1e-2))$par
}
r.gaussian <- replicate(100, wrap("gaussian"))
wrap("poisson")
wrap("Gamma")