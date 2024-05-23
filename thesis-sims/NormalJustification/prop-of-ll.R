# Corrections -------------------------------------------------------------
# Working out proportion of log-likelihood explained by Y / b / {T_\D}.
source('Sample.R')

# Let's rip from Sample.R and return ll related things...
cSample <- function(X_, TUNE = 1., b.dist = 'normal', df = 4, force.intslope = FALSE, burnin = 1000, NMC = 10000){
  # Check
  stopifnot(inherits(X_, "dataGen"))
  
  # Get the data, survival data (maybe not used?) and true random effects ----
  data <- X_$data; surv.data <- X_$surv.data; btrue <- X_$btrue
  
  # Unpack true values used to generate {data, btrue} ----
  ids <- X_$ids; args <- X_$args
  family <- args$family[[1]]
  beta <- c(args$beta[1,])
  sigma <- args$sigma[[1]]
  D <- args$D
  gamma <- args$gamma; zeta <- args$zeta
  theta <- args$theta
  gamma.rep <- rep(gamma, ncol(btrue))
  
  # Create req. data objects ----
  W <- X <- Y <- Z <- setNames(vector('list', length(ids)), paste0("id: ", ids))
  for(i in seq_along(ids)){
    X[[i]] <- W[[i]] <- Y[[i]] <- Z[[i]] <- list() # This for ease of use with `gmvjoint`.
    X[[i]][[1]] <- model.matrix(~time+cont+bin, data[data$id==ids[i],,drop=F])
    if(force.intslope){
      Z[[i]][[1]] <- model.matrix(~time, data[data$id==ids[i],,drop=F])
    }else{
      if(!family%in%c("binomial", "genpois"))
        Z[[i]][[1]] <- model.matrix(~time, data[data$id==ids[i],,drop=F])
      else
        Z[[i]][[1]] <- model.matrix(~1, data[data$id==ids[i],,drop=F])
    }
    Y[[i]][[1]] <- data[data$id==ids[i], 'Y.1']
    W[[i]][[1]] <- matrix(1, nrow = nrow(X[[i]][[1]]), ncol = 1)
  }
  # Random effects list
  b <- lapply(seq_along(ids), function(i) btrue[ids[i],,drop=F])
  b.inds.cpp <- list(0:(ncol(btrue) - 1)) # And indices
  # m
  m <- sapply(unique(X_$data$id), function(i) nrow(X_$data[X_$data$id == i,,drop=F]))
  
  # Survival part ----
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data)
  fts <- sort(surv$ft)
  l0 <- exp(theta[1] + theta[2] * fts) # exp(log(nu)) * exp(alpha * t)
  # Hardcoding for genpois and binomial case, where only random intercept is
  # fitted/modelled.
  if(family %in% c("genpois", "binomial") && !force.intslope){
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1 ~ time + cont + bin + (1|id)), gmvjoint:::parseFormula), l0)
  }else{
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1 ~ time + cont + bin + (1 + time|id)), gmvjoint:::parseFormula), l0)  
    b.inds <- list(0:1)
    gamma.rep <- rep(gamma, 2)
  }
  
  # Omega^{(TRUE)} ---->
  Omega <- list(D = D, beta = beta, sigma = list(sigma), gamma = gamma, gamma.rep = gamma.rep, zeta = zeta[zeta!=0L])
  
  # Tuning parameters; aiming for about 22.5% acceptance across all subjects.
  tune <- TUNE
  
  # Obtain _all_ generated survival data objects ----
  Delta <- lapply(seq_along(ids), function(x) surv$Delta[[ids[x]]])
  S <- lapply(seq_along(ids), function(x) sv$S[[ids[x]]])
  Fi <- lapply(seq_along(ids), function(x) sv$Fi[[ids[x]]])
  l0i <- lapply(seq_along(ids), function(x) sv$l0i[[ids[x]]])
  SS <- lapply(seq_along(ids), function(x) sv$SS[[ids[x]]])
  Fu <- lapply(seq_along(ids), function(x) sv$Fu[[ids[x]]])
  l0u <- lapply(seq_along(ids), function(x) sv$l0u[[ids[x]]])
  
  # Get Sigma
  Sigma <- Map(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    Sigma <- getSigma(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds.cpp, Omega, TRUE)
  }, b = b, Y = Y, X = X, Z = Z, W = W, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, l0u = l0u)
  
  b.hat <- lapply(Sigma, el, 1)
  Sigma <- lapply(Sigma, el, 2)
  
  cli::cli_progress_bar(name = "Sampling...", total = length(ids))
  norm.dens <- cond.dens <- MVN.dens <- vector("list", length(ids))
  Walks <- norm.dens
  q <- ncol(D)
  Acc <- numeric(length(ids))
  for(a in 1:length(ids)){
    cond.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    norm.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    # I think it should be MH scheme, not just Met; give _very_ similar results though.
    cond.sim <- Metropolis_Hastings(b[[a]], Omega, Y[[a]], X[[a]], Z[[a]], W[[a]],
                                    list(family), Delta[[a]], S[[a]], Fi[[a]], l0i[[a]], 
                                    SS[[a]], Fu[[a]], l0u[[a]], Omega$gamma.rep,
                                    list(0:3), b.inds.cpp, 1L, q, burnin, NMC, Sigma[[a]], 
                                    b.dist, df, tune)
    Acc[a] <- cond.sim$AcceptanceRate
    for(j in 1:q){
      cond.dens[[a]][[j]] <- density(t(cond.sim$walks)[,j])
      norm.dens[[a]][[j]] <- dnorm(cond.dens[[a]][[j]]$x, mean = b.hat[[a]][j], 
                                   sd = sqrt(Sigma[[a]][j,j]))
    }
    MVN.dens[[a]] <- mvtnorm::dmvnorm(t(cond.sim$walks), mean = b.hat[[a]], sigma = Sigma[[a]])
    Walks[[a]] <- t(cond.sim$walks)
    rm(cond.sim)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  b.post <- t(sapply(Walks, function(x) apply(x,2,median)))
  b.post <- lapply(1:nrow(b.post), function(i) b.post[i,,drop=F])
  # Work out logSurv -> value of log-likelihood of survival part
  logSurv <- mapply(function(b,S, SS, Fi, Fu, l0i, l0u, Delta){
      gmvjoint:::logfti(b, S, SS, Fi, Fu, l0i, l0u, Delta, Omega$gamma.rep, Omega$zeta)
  }, b = b.post, S = S, SS = SS, Fi = Fi, Fu = Fu, l0i = l0i, l0u = l0u, Delta = Delta)
  logOverall <- mapply(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    -1 * gmvjoint:::joint_density(b, Y, X, Z, W, Omega$beta, Omega$D, Omega$sigma,
                                  as.list(family), Delta, S, Fi, l0i, SS, Fu, l0u, 
                                  Omega$gamma.rep, Omega$zeta, list(0:3), b.inds.cpp, 1L)
  }, b = b.post, Y = Y, X = X, Z = Z, W = W, Delta = Delta, 
     S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, l0u = l0u)
  # What proportion of the joint density is the survival part accounting for?
  prop.surv <- logSurv/logOverall * 100
  prop.long <- (logOverall-logSurv)/logOverall * 100
  cat(sprintf("Family: %s\n", family)) # just to check!
  cat("Median acceptance rate: ", round(median(Acc), 3), '\n')
  cbind(logSurv, logOverall, prop.surv, prop.long, m)
}


# Investigation per family ------------------------------------------------
# Gaussian
sim <- create.simfn()
X_ <- dataGen(sim)
X <- cSample(X_, 5.5)
apply(X,2,quantile)
#         logSurv  logOverall prop.surv   prop.long  m
# 0%   -4.1503850 -13.9112533  0.000000   0.7385872  1
# 25%  -0.9627334  -7.1962178  6.234129  77.7458925  8
# 50%  -0.5857602  -5.0846069 12.061549  87.9384511 10
# 75%  -0.3177149  -3.8242518 22.254108  93.7658712 10
# 100%  0.0000000  -0.1572439 99.261413 100.0000000 10

# Poisson
sim <- create.simfn(list('poisson'))
X_ <- dataGen(sim)
X <- cSample(X_, 5.75)
apply(X,2,quantile)
#         logSurv logOverall prop.surv prop.long  m
# 0%   -3.9666111 -36.856400  0.000000  35.10684  1
# 25%  -1.0976296 -27.969180  2.152718  96.12354 10
# 50%  -0.6624869 -23.924107  2.745927  97.25407 10
# 75%  -0.4594787 -17.815643  3.876456  97.84728 10
# 100%  0.0000000  -2.483314 64.893160 100.00000 10

# Binomial
sim <- create.simfn(list("binomial"), arguments = list(D = diag(c(2,.5))))
X_ <- dataGen(sim)
X <- cSample(X_, 5.5, force.intslope = T)
apply(X,2,quantile)
#         logSurv logOverall prop.surv prop.long  m
# 0%   -4.9664099 -14.047860  0.000000  37.20965  1
# 25%  -2.9585731  -8.032236  6.439259  67.81841  7
# 50%  -0.5258749  -5.855739 11.212765  88.78724 10
# 75%  -0.2391967  -3.843612 32.181593  93.56074 10
# 100%  0.0000000  -2.183144 62.790347 100.00000 10
