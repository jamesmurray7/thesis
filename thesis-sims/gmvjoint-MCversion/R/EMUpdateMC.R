EMupdate <- function(Omega, family, dmats, b,                # Params; families; dmats; REs;
                     sv, surv, MCtype, N,                    # Survival; MC stuff;
                     con, inds, XTX){                        # control arguments; indices
  
  # Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  
  # Find b.hat and Sigma ==================================================
  b.update <- Map(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma, family = family, 
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
          beta_inds = inds$Cpp$beta, b_inds = inds$Cpp$b, K = dmats$K,
          method = 'BFGS', hessian = TRUE)
  }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, W = dmats$W, Delta = surv$Delta, 
     S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, Fu = sv$Fu, l0u = sv$l0u)
  
  b.hat <- lapply(b.update, el, 1)
  Sigma <- lapply(b.update, function(x) solve(x$hessian))
  
  # Split into K constituent sub-matrices.
  SigmaSplit <- lapply(Sigma, function(x) lapply(inds$R$b, function(y) as.matrix(x[y,y])))
  
  # Draws for E-step ======================================================
  if(MCtype == "ordinary"){
    draws <- Map(function(b.hat, Sigma.hat){
      MASS::mvrnorm(N, b.hat, Sigma.hat)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else if(MCtype == "antithetic"){
    draws <- joineRML:::bSim(floor(N/2), b.hat, Sigma)
  }else if(MCtype == "sobol"){
    zzz <- suppressWarnings(randtoolbox::sobol(N, sum(dmats$q), normal = T, scrambling = 1))
    draws <- Map(function(b.hat, Sigma.hat){
      C <- chol(Sigma.hat)
      matrix(rep(b.hat, N), N, byrow = T) + (zzz %*% C)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else{
    stop("999")
  }
  
  # E-step ================================================================
  
  # E[b_i] ====
  Eb <- lapply(draws, colMeans)
  
  # E[\eta] ====
  Eeta <- Map(function(X, Z, b) make_eta(X, Z, Omega$beta, b, inds$Cpp$beta, inds$Cpp$b), 
              X = dmats$X, Z = dmats$Z, b = Eb)
  
  # E[exp{\eta}] ====
  Eexpeta <- Map(function(X, Z, draws) make_Expeta(X, Z, Omega$beta, draws, inds$Cpp$beta, inds$Cpp$b),
                 X = dmats$X, Z = dmats$Z, draws = draws)
  Eexpeta.overetas <- Map(function(Y, X, Z, draws) make_binExp(Y[[3]], X[[3]], Z[[3]], beta[inds$R$beta[[3]]], draws[,inds$R$b[[3]], drop = F]),
                          Y = dmats$Y, X = dmats$X, Z = dmats$Z, draws = draws)
  
  # E[(Y-\eta)^\top(Y-\eta)] ====
  # Hardcoded to be the first response...
  a <- mapply(function(Y, X, Z, draws) Eymeta(Y[[1]], X[[1]], Z[[1]], Omega$beta[1:4], draws[,1:2]),
              Y = dmats$Y, X = dmats$X, Z = dmats$Z, draws = draws)
  
  # E[<survival part>] ====
  Esurvpart <- Map(function(draws, SS, Fu) make_Esurvpart(draws, SS, Fu, gamma.rep, zeta), draws = draws, SS = sv$SS, Fu = sv$Fu)
  
  # M - step ==============================================================
  
  # D ======================================
  D.update <- Map(function(Sigma, Eb) Sigma + tcrossprod(Eb), Sigma = Sigma, Eb = Eb)
  D.new <- Reduce('+', D.update)/dmats$n
  
  # \beta ==================================
  # (k=1, gaussian)
  Gauss.beta.update.rhs <- rowSums(mapply(function(X,Y,Z,b){
    crossprod(X[[1]], Y[[1]] - Z[[1]] %*% b[1:2])
  },X=dmats$X, Y = dmats$Y, Z = dmats$Z, b = Eb))
  Gauss.beta.update <- solve(XTX, Gauss.beta.update.rhs)
  
  # (k=2, poisson)
  Poisson.score <- rowSums(mapply(function(X, Y, Eexpeta) crossprod(X[[2]], Y[[2]] - Eexpeta[[2]]), X = dmats$X, Y = dmats$Y, Eexpeta = Eexpeta))
  Poisson.diaghess <- lapply(Eexpeta, function(e) c(e[[2]] * -1))
  Poisson.Hessian <- Reduce("+", Map(function(X, d){
    crossprod(X[[2]], diag(d, nrow = length(d), ncol = length(d)) %*% X[[2]])
  }, X = dmats$X, d = Poisson.diaghess))
  Poisson.beta.update <- beta[inds$R$beta[[2]]] - solve(Poisson.Hessian, Poisson.score)
  
  # (k=3, binomial)
  Bin.score <- rowSums(mapply(function(Y, X, aa) crossprod(X[[3]], Y[[3]] - aa[[1]]), Y = dmats$Y, X = dmats$X, aa = Eexpeta.overetas))
  Bin.diaghess <- lapply(Eexpeta.overetas, function(e) c(e[[2]] * -1))
  Bin.Hessian <- Reduce("+", Map(function(X, d){
    crossprod(X[[2]], diag(d, nrow = length(d), ncol = length(d)) %*% X[[2]])
  }, X = dmats$X, d = Bin.diaghess))
  Bin.beta.update <- beta[inds$R$beta[[3]]] - solve(Bin.Hessian, Bin.score)
  
  beta.new <- setNames(c(Gauss.beta.update, Poisson.beta.update, Bin.beta.update), names(beta))
  
  # sigma^2_\varepsilon ====================
  sigma2.update <- sum(a)/dmats$m[1]
  
  # Survival parameters (\gamma, \zeta) ====
  # (for profile update to \Phi).
  uu <- length(sv$ft)
  Esurvpart2 <- do.call(cbind, lapply(Esurvpart, function(e){
    e <- c(e)
    ui <- length(e)
    if(ui == uu) 
      return(e)
    else
      return(c(e, rep(0, uu-ui)))
  }))
  
  l0.hat <- sv$nev/rowSums(Esurvpart2)
  l0u.hat <- lapply(sv$l0u, function(ll){
    l0.hat[1:length(ll)]
  })
  
  Psurv <- length(c(gamma, zeta))
  
  Sgz <- mapply(function(draws, Eb, S, SS, Fu, Fi, l0u, Delta, ST){
    if(length(ST))
      return(Sgammazeta(c(gamma, zeta), draws, Eb, S, SS, Fu, Fi, l0u, Delta,
                        inds$Cpp$b, dmats$K, .Machine$double.eps^(1/3)))
    else
      return(rep(0, Psurv))
  }, draws = draws, Eb = Eb, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
     Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)
  
  .SS <- rowSums(Sgz)
  .SiSiT <- Reduce('+', lapply(1:dmats$n, function(i) tcrossprod(Sgz[, i])))
  .H <- .SiSiT - tcrossprod(.SS)/dmats$n
  
  # Hgz is really slow, so use a GN update is okay here.
  # Hgz <- Map(function(draws, Eb, S, SS, Fu, Fi, l0u, Delta, ST){
  #   if(length(ST))
  #     return(Hgammazeta(c(gamma, zeta), draws, Eb, S, SS, Fu, Fi, l0u, Delta, 
  #                         inds$Cpp$b, dmats$K, .Machine$double.eps^(1/4)))
  #   else
  #     return(matrix(0, Psurv, Psurv))
  # }, draws = draws, Eb = Eb, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  # Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)

  # Survival parameters (gamma, zeta)
  # gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz))
  gammazeta.new <- c(gamma, zeta) + 0.5 * solve(.H, .SS)
  gamma.new <- gammazeta.new[1:dmats$K]; zeta.new <- gammazeta.new[(dmats$K+1):length(gammazeta.new)]
  
  # The baseline hazard and related objects
  lambda.update <- c(lambda_update(draws, sv$Fu, sv$SS, sv$surv.times, 
                                   rep(gamma.new, sapply(inds$R$b, length)), zeta.new, sv$nev))
  
  l0u.new <- lapply(sv$l0u, function(ll){
    lambda.update[1:length(ll)]
  })
  l0i.new <- lambda.update[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Return
  list(
    D = D.new, beta = beta.new, sigma = list(sigma2.update,0,0),#   K responses;
    gamma = gamma.new, zeta = zeta.new,                         #   survival;
    l0 = lambda.update, l0u = l0u.new, l0i = as.list(l0i.new),  #   hazard;
    l0u.hat = l0u.hat,                                          #   (+ profile);
    b = b.hat,                                                  #   REs;
    Sigma = Sigma                                               #   + their variance.
  )
  
}
