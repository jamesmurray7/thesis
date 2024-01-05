# 1) Fit on simulated bivariate data, (1x gaussian, 1x poisson) --------
beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
gamma <- c(0.3, -0.3)
D <- diag(c(0.25, 0.09, 0.25, 0.05))
family <- list('gaussian', 'poisson')
data <- simData(ntms = 10, beta = beta, D = D, n = 250,
                family = family, zeta = c(0, -0.2),
                sigma = list(0.16, 0), gamma = gamma, theta = c(-3,.1))$data

# Specify formulae and target families
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id)   # Poisson
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
D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
gamma.rep <- rep(gamma, sapply(inds$R$b, length))
XTX <- lapply(dmats$X, function(x) crossprod(x[[1]]))
XTX <- Reduce("+", XTX)

# Find b.hat and Sigma ====================================================
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

# Split into K constituent sub-matrices (pointless for this e.g.)
# SigmaSplit <- lapply(Sigma, function(x) lapply(inds$R$b, function(y) as.matrix(x[y,y])))

# E step ==================================================================
# Draws for E-step
draws <- Map(function(b.hat, Sigma.hat){
  MASS::mvrnorm(N, b.hat, Sigma.hat) # Sobol?
}, b.hat = b.hat, Sigma.hat = Sigma)

# E[b_i b_i^\top] ====
EbbT <- lapply(draws, sum_tcrossprods) 
EbbT <- lapply(EbbT, '/', N)
# Don't really need this as I wrote down M steps wrong, oops

# E[b_i] ====
Eb <- lapply(draws, colMeans)

# E[\eta] ====
Eeta <- Map(function(X, Z, b) make_eta(X, Z, beta, b, inds$Cpp$beta, inds$Cpp$b), 
            X = dmats$X, Z = dmats$Z, b = Eb)

# E[exp{\eta}] ====
Eexpeta <- Map(function(X, Z, draws) make_Expeta(X, Z, beta, draws, inds$Cpp$beta, inds$Cpp$b),
               X = dmats$X, Z = dmats$Z, draws = draws)

# E[(Y-\eta)^\top(Y-\eta)] ====
# Hardcoded to be the first response...
a <- mapply(function(Y, X, Z, draws) Eymeta(Y[[1]], X[[1]], Z[[1]], beta[1:4], draws[,1:2]),
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

beta.update <- setNames(c(Gauss.beta.update, Poisson.beta.update), names(beta))

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

# Some comparisons ========================================================
rbind("GH" = vech(update.GH$D), "MC" = vech(D.update))
rbind("GH" = c(update.GH$beta), "MC" = beta.update)
c("GH" = update.GH$sigma[[1]], "MC" = sigma2.update)
rbind("GH" = c(update.GH$gamma, update.GH$zeta),
      "MC" = gammazeta.update)
