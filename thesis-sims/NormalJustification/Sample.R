# Obtain b.hat and Sigma.hat given observed data at TRUE parameter estimates.
getSigma <- function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, Omega){
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta[Omega$zeta!=0L]
  gamma.rep <- Omega$gamma.rep
  ln <- nrow(X[[1]])
  uu <- optim(b, gmvjoint:::joint_density, gmvjoint:::joint_density_ddb,
              Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma,
              family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
              gamma_rep = gamma.rep, zeta = zeta, beta_inds = list(0:3), b_inds = b.inds, K = 1L, 
              method = 'BFGS', hessian = T)
  list(bhat = uu$par, Sigma = solve(uu$hessian))
}

# Main "Sample" function --------------------------------------------------
# Function to sample given data, true random effects, known family and target ids.
# `X` is an object of class `DataGen`
Sample <- function(X_, TUNE = 1., return.walks = FALSE){
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
    if(!family%in%c("binomial", "genpois"))
      Z[[i]][[1]] <- model.matrix(~time, data[data$id==ids[i],,drop=F])
    else
      Z[[i]][[1]] <- model.matrix(~1, data[data$id==ids[i],,drop=F])
    Y[[i]][[1]] <- data[data$id==ids[i], 'Y.1']
    W[[i]][[1]] <- matrix(1, nrow = nrow(X[[i]][[1]]), ncol = 1)
  }
  # Random effects list
  b <- lapply(seq_along(ids), function(i) btrue[ids[i],,drop=F])
  b.inds.cpp <- list(0:(ncol(btrue) - 1)) # And indices
  
  # Survival part ----
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data)
  fts <- sort(surv$ft)
  l0 <- exp(theta[1] + theta[2] * fts) # exp(log(nu)) * exp(alpha * t)
  # Hardcoding for genpois and binomial case, where only random intercept is
  # fitted/modelled.
  if(family %in% c("genpois", "binomial")){
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1 ~ time + cont + bin + (1|id)), gmvjoint:::parseFormula), l0)
  }else{
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1 ~ time + cont + bin + (1 + time|id)), gmvjoint:::parseFormula), l0)  
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
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
    Sigma <- getSigma(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds.cpp, Omega)
  }, b = b, Y = Y, X = X, Z = Z, W = W, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, 
     Fu = Fu, l0u = l0u)
  
  b.hat <- lapply(Sigma, el, 1)
  Sigma <- lapply(Sigma, el, 2)
  
  cli::cli_progress_bar(name = "Sampling...", total = length(ids))
  norm.dens <- cond.dens <- MVN.dens <- vector("list", length(ids))
  if(return.walks) Walks <- norm.dens
  q <- ncol(D)
  Acc <- numeric(length(ids))
  for(a in 1:length(ids)){
    cond.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    norm.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    # Walks
    cond.sim <- gmvjoint:::metropolis(b[[a]], Omega, Y[[a]], X[[a]], Z[[a]], W[[a]],
                                  list(family), Delta[[a]], S[[a]], Fi[[a]], l0i[[a]], 
                                  SS[[a]], Fu[[a]], l0u[[a]], Omega$gamma.rep,
                                  list(0:3), b.inds.cpp, 1L, q, 1000, 10000, Sigma[[a]], tune)
    Acc[a] <- cond.sim$AcceptanceRate
    for(j in 1:q){
      cond.dens[[a]][[j]] <- density(t(cond.sim$walks)[,j])
      norm.dens[[a]][[j]] <- dnorm(cond.dens[[a]][[j]]$x, mean = b.hat[[a]][j], 
                                   sd = sqrt(Sigma[[a]][j,j]))
    }
    MVN.dens[[a]] <- mvtnorm::dmvnorm(t(cond.sim$walks), mean = b.hat[[a]], sigma = Sigma[[a]])
    if(return.walks) Walks[[a]] <- t(cond.sim$walks)
    rm(cond.sim)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  # Make a data.frame to export.
  mi <- sapply(Z, function(x) nrow(x[[1]]))
  
  dfs <- lapply(1:length(ids), function(i){
    id <- i
    m <- mi[i]
    b0 <- data.frame(condx = cond.dens[[i]]$b0$x, condy = cond.dens[[i]]$b0$y,
                     normy = norm.dens[[i]]$b0, hatb = b.hat[[i]][1],
                     family = family, var = "b[0]")
    if(q > 1){
      b1 <- data.frame(condx = cond.dens[[i]]$b1$x, condy = cond.dens[[i]]$b1$y,
                       normy = norm.dens[[i]]$b1, hatb = b.hat[[i]][2],
                       family = family, var = "b[1]")
      i.dens <- rbind(b0,b1)
    }else{
      i.dens <- b0
    }
    i.dens$id <- i; i.dens$m <- mi[i]
    i.dens
  })
  
  out <- do.call(rbind, dfs) # About 5MB at 10,000 iterations on 100 subjects.
  out$ApproxBias <- out$normy - out$condy
  out <- list(df = out, MVN.dens = MVN.dens, mi = unname(mi),
              family = family, Sigmas = Sigma, true.b = btrue,
              Acc = Acc)
  if(return.walks) out$Walks <- Walks
  class(out) <- "Sample"
  out
}

print.Sample <- function(x){
  stopifnot(inherits(x, "Sample"))
  qAcc <- quantile(x$Acc, probs = c(.5, .25, .75))
  cat(sprintf("\nAcceptance: %.4f [%.4f, %.4f]\n", qAcc[1], qAcc[2], qAcc[3]))
  print(head(x$df))
  cat("\n")
}

trim.Walks <- function(x){
  if(inherits(x, "Sample")){ # Either supply a `Sample` object...
    stopifnot(!is.null(x$Walks))
    W <- x$Walks
    if(class)
      out <- lapply(W, function(w){
        w[!duplicated.matrix(w),]
      })
  }else{ # OR a standalone matrix of random effect walks...
    out <- w[!duplicated.matrix(w),]
  }
  out
}
