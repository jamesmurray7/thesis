rm(list=ls())
# Trivariate fit as per Sec 4.5.3
# beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
beta <- rbind(
  c(-2, 0.1, -0.1, 0.2),
  c(2, -0.1, 0.1, 0.2),
  c(1, -1, 1, -1)
)
gamma <- c(0.5, -0.5, .5)
D <- diag(c(0.25, 0.09, 0.50, 0.09, 2))
D[1,3] <- D[3,1] <- D[1,5] <- D[5,1] <- D[3,5] <- D[5,3] <- .125
family <- list('gaussian', 'poisson', "binomial")
data <- simData(ntms = 10, beta = beta, D = D, n = 250,
                family = family, zeta = c(0, -0.2),
                random.formulas = list(~time,~time,~1),
                sigma = list(0.16, 0, 0), gamma = gamma, theta = c(-3,.1))$data

# Specify formulae and target families
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id),  # Poisson
  Y.3 ~ time + cont + bin + (1|id)          # Binomial
)
surv.formula <- Surv(survtime, status) ~ bin

control <- list()
disp.formulas = NULL

# Run joint.R up to EMUpdate call
update.GH <- update # update in ls() is from joint.R/EMUpdate.R obtained by GHQ.

# MC version below --------------------------------------------------------
N <- 5e3
# sourceCpp("testMC.cpp")
# Unpack Omega, the parameter vector
# D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
XTX <- lapply(dmats$X, function(x) crossprod(x[[1]]))
XTX <- Reduce("+", XTX)

get.update <- function(MCtype = "montecarlo", Omega, l0u, l0i){
  # Unpack Omega ====
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  # Find b.hat and Sigma ====================================================
  b.update <- Map(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, W = W, beta = Omega$beta, D = D, sigma = Omega$sigma, family = family, 
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, 
          gamma_rep = gamma.rep, zeta = Omega$zeta,
          beta_inds = inds$Cpp$beta, b_inds = inds$Cpp$b, K = dmats$K,
          method = 'BFGS', hessian = TRUE)
  }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, W = dmats$W, Delta = surv$Delta, 
  S = sv$S, Fi = sv$Fi, l0i = as.list(l0i), SS = sv$SS, Fu = sv$Fu, l0u = l0u)
  
  b.hat <- lapply(b.update, el, 1)
  Sigma <- lapply(b.update, function(x) solve(x$hessian))
  
  # E step ==================================================================
  # Draws for E-step
  if(MCtype == "montecarlo"){
    draws <- Map(function(b.hat, Sigma.hat){
      MASS::mvrnorm(N, b.hat, Sigma.hat)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else if(MCtype == "antithetic"){
    draws <- joineRML:::bSim(floor(N/2), b.hat, Sigma)
  }else if(MCtype == "sobol"){
    zzz <- randtoolbox::sobol(N, sum(dmats$q), normal = T, scrambling = 1)
    draws <- Map(function(b.hat, Sigma.hat){
      C <- chol(Sigma.hat)
      matrix(rep(b.hat, N), N, byrow = T) + (zzz %*% C)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else{
    stop("999")
  }
  
  # E[b_i b_i^\top] ====
  EbbT <- lapply(draws, sum_tcrossprods) 
  EbbT <- lapply(EbbT, '/', N)
  # Don't really need this as I wrote down M steps wrong, oops
  
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
  
  # M - step ================================================================
  
  # D ====
  D.update <- Map(function(Sigma, Eb) Sigma + tcrossprod(Eb), Sigma = Sigma, Eb = Eb)
  D.update <- Reduce('+', D.update)/n
  
  # \beta ====
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
  
  beta.update <- setNames(c(Gauss.beta.update, Poisson.beta.update, Bin.beta.update), names(beta))
  
  # sigma^2_\varepsilon ====
  sigma2.update <- sum(a)/dmats$m[1]
  
  # \hat{\lambda} ====
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
  
  # \Phi ====
  Esurvpart <- function(gammazeta){
    gamma <- gammazeta[1:dmats$K]
    gamma.repp <- rep(gamma, sapply(inds$R$b, length))
    zeta <- gammazeta[(dmats$K+1):length(gammazeta)]
    Esp <- Map(function(draws, SS, Fu) make_Esurvpart(draws, SS, Fu, gamma.repp, zeta), draws = draws, SS = sv$SS, Fu = sv$Fu)
    sum(mapply(function(Del, S, Fi, Eb, haz, Esp){
      Del * (S %*% zeta + Fi %*% (gamma.repp * Eb)) - crossprod(haz, Esp)
    }, Del = surv$Delta, S = sv$S, Fi = sv$Fi, Eb = Eb, haz = l0u.hat, Esp = Esp))
  }
  
  Sgz <- pracma::grad(Esurvpart, c(gamma,zeta))
  Hgz <- pracma::hessian(Esurvpart, c(gamma,zeta))
  gammazeta.update <- c(gamma, zeta) - solve(Hgz, Sgz)
  gamma.update <- gammazeta.update[1:dmats$K]
  zeta.update <-  gammazeta.update[(dmats$K+1):length(gammazeta.update)]
  
  # The baseline hazard and related objects
  Esurvpart <- Map(function(draws, SS, Fu) make_Esurvpart(draws, SS, Fu, rep(gamma.update, sapply(inds$R$b, length)), zeta.update), draws = draws, SS = sv$SS, Fu = sv$Fu)
  Esurvpart2 <- do.call(cbind, lapply(Esurvpart, function(e){
    e <- c(e)
    ui <- length(e)
    if(ui == uu) 
      return(e)
    else
      return(c(e, rep(0, uu-ui)))
  }))
  
  lambda.update <- sv$nev/rowSums(Esurvpart2)
  l0u.new <- lapply(sv$l0u, function(ll){
    lambda.update[1:length(ll)]
  })
  l0i.new <- lambda.update[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Some comparisons ========================================================
  cat("\nvech(D)----\n")
  print(rbind("GH" = setNames(vech(update.GH$D),
                              paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')), 
              "MC" = setNames(vech(D.update),
                              paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']'))))
  cat("\nbeta----\n")
  print(rbind("GH" = c(update.GH$beta), "MC" = beta.update))
  cat("\nsigma^2----\n")
  print(c("GH" = update.GH$sigma[[1]], "MC" = sigma2.update))
  cat("\nzeta----\n")
  print(rbind("GH" = c(update.GH$gamma, update.GH$zeta),
        "MC" = gammazeta.update))
  cat("\n----\n")
  
  return(list(
    Omega = list(D = D.update,
                 beta = beta.update,
                 sigma = list(sigma2.update,0,0),
                 gamma = gamma.update, zeta = zeta.update),
    l0u = l0u.new,
    l0i = l0i.new
  ))
}

Omega0 <- Omega
OMC <- get.update(Omega = Omega0, l0u = sv$l0u, l0i = sv$l0i)
OMC2 <- get.update(Omega = OMC$Omega, l0u = OMC$l0u, l0i = OMC$l0i)
AMC <- get.update("antithetic", Omega = Omega0, l0u = sv$l0u, l0i = sv$l0i)
AMC2 <- get.update("antithetic", Omega = AMC$Omega, l0u = AMC$l0u, l0i = AMC$l0i)
QMC <- get.update("sobol", Omega = Omega0, l0u = sv$l0u, l0i = sv$l0i)
QMC2 <- get.update("sobol", Omega = QMC$Omega, l0u = QMC$l0u, l0i = QMC$l0i)

all.updates <- list(OMC = OMC$Omega, OMC2 = OMC2$Omega,
                    AMC = AMC$Omega, AMC2 = AMC2$Omega,
                    QMC = QMC$Omega, QMC2 = QMC2$Omega,
                    GH = list(D = update.GH$D,
                              beta = update.GH$beta,
                              sigma = update.GH$sigma,
                              gamma = update.GH$gamma,
                              zeta = update.GH$zeta),
                    GH2 = list(D = update.GH2$D,
                               beta = update.GH2$beta,
                               sigma = update.GH2$sigma,
                               gamma = update.GH2$gamma,
                               zeta = update.GH2$zeta))

tt <- function(x){
  vd <- vech(x$D); beta <- x$beta; sigma <- x$sigma[[1]]; gamma <- x$gamma; zeta <- x$zeta
  format(round(c(vd, beta, sigma, gamma, zeta), 4), nsmall = 4, justify = 'right')
}
tt2 <- function(x){
  vd <- vech(x$D); beta <- x$beta; sigma <- x$sigma[[1]]; gamma <- x$gamma; zeta <- x$zeta
  c(vd, beta, sigma, gamma, zeta)
}

all.updates2 <- lapply(all.updates, tt)
ttt <- cbind(p = names(params), GHQ = all.updates2$GH, OMC = all.updates2$OMC, AMC = all.updates2$AMC, QMC = all.updates2$QMC,
             p = names(params), GHQ = all.updates2$GH2, OMC = all.updates2$OMC2, AMC = all.updates2$AMC2, QMC = all.updates2$QMC2)

library(xtable)
xt <- xtable(ttt)
print(xt, include.rownames = F, sanitize.text.function = identity)

# Biggest difference at each iteration
iter0 <- cbind(GH = tt2(all.updates$GH), OMC = tt2(all.updates$OMC), AMC = tt2(all.updates$AMC), QMC = tt2(all.updates$QMC))
iter1 <- cbind(GH = tt2(all.updates$GH2), OMC = tt2(all.updates$OMC2), AMC = tt2(all.updates$AMC2), QMC = tt2(all.updates$QMC2))

# Quick absolute difference
absdiff0 <- apply(iter0, 1, function(x){
  omc <- abs(x[1]-x[2]); amc <- abs(x[1]-x[3]); qmc <- abs(x[1]-x[4])
  max(c(omc,amc,qmc))
})
ind <- which.max(absdiff0)
absdiff0[ind]
params[ind]
iter0[ind,]

absdiff0 <- apply(iter1, 1, function(x){
  omc <- abs(x[1]-x[2]); amc <- abs(x[1]-x[3]); qmc <- abs(x[1]-x[4])
  max(c(omc,amc,qmc))
})
ind <- which.max(absdiff0)
absdiff0[ind]
params[ind]
iter1[ind,]


