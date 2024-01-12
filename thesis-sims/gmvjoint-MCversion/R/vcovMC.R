#' @keywords internal
obs.emp.I <- function(Omega, dmats, surv, sv, family,
                      b, l0i, l0u, MCtype, N, sf, inds, con){
  # Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma); gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  zeta <- c(Omega$zeta)
  
  # Calculate b.hat, Sigma.hat at MLEs
  b.update <- Map(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma, family = family, 
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
          beta_inds = inds$Cpp$beta, b_inds = inds$Cpp$b, K = dmats$K,
          method = 'BFGS', hessian = TRUE)
  }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, W = dmats$W, Delta = surv$Delta, 
  S = sv$S, Fi = sv$Fi, l0i = l0i, SS = sv$SS, Fu = sv$Fu, l0u = l0u)
  
  Sigma <- lapply(b.update, function(x) solve(x$hessian))
  b.hat <- lapply(b.update, el, 1)
  
  if(MCtype == "ordinary"){
    draws <- Map(function(b.hat, Sigma.hat){
      MASS::mvrnorm(N, b.hat, Sigma.hat * sf)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else if(MCtype == "antithetic"){
    draws <- joineRML:::bSim(floor(N/2), b.hat, lapply(Sigma, '*', sf))
  }else if(MCtype == "sobol"){
    zzz <- randtoolbox::sobol(N, sum(dmats$q), normal = T, scrambling = 1)
    draws <- Map(function(b.hat, Sigma.hat){
      C <- chol(Sigma.hat * sf)
      matrix(rep(b.hat, N), N, byrow = T) + (zzz %*% C)
    }, b.hat = b.hat, Sigma.hat = Sigma)
  }else{
    stop("999")
  }
  
  # Some expectations
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
  
  # For sigma
  a <- mapply(function(Y, X, Z, draws) Eymeta(Y[[1]], X[[1]], Z[[1]], Omega$beta[inds$R$beta[[1]]], draws[,inds$R$b[[1]]]),
              Y = dmats$Y, X = dmats$X, Z = dmats$Z, draws = draws)
  
  # E[<survival part>] ====
  Esurvpart <- Map(function(draws, SS, Fu) make_Esurvpart(draws, SS, Fu, gamma.rep, zeta), draws = draws, SS = sv$SS, Fu = sv$Fu)
  
  
  # Profile estimate for (gamma, zeta) at MLEs and b.hat
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
  
  # Scores ------------------------------------------------------------------
  # D =========================================
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i i=1,...,length(vech(D))
    out
  })
  
  sDi <- function(i) {
    mapply(function(b, S) {
      EbbT <- S + tcrossprod(b)
      invDdD <- Dinv %*% delta.D[[i]]
      invDdDinvD <- invDdD %*% Dinv
      0.5 * sum(diag(invDdDinvD %*% EbbT)) - 0.5 * sum(diag(invDdD))
    },
    b = Eb, S = Sigma,
    SIMPLIFY = TRUE)
  }

  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
  # Below gets the same result
  # postmult <- diag(1, nrow = sv$q, ncol = sv$q)
  # postmult[postmult == 0] <- 2
  # EbbT <- Map(function(b, S) S + tcrossprod(b), b = b.hat, S = Sigma)
  # term <- lapply(EbbT, function(x) 0.5 * Dinv %*% x %*% Dinv - 0.5 * Dinv)
  # sD <- lapply(term, function(x) vech(0.5 * (t(x) + x)) * vech(postmult))
  
  # \beta =====================================
  Sb <- Map(function(X, Y, Z, W, Eb, Eeta, Eexpeta, Eexpeta.overetas){
    # Gaussian
    V <- diag(nrow = length(Y[[1]]), ncol = length(Y[[1]])) * sigma[[1]]
    Sg <- crossprod(X[[1]], solve(V, Y[[1]] - Eeta[[1]]))
    # Poisson
    Sp <- crossprod(X[[2]], Y[[2]] - Eexpeta[[2]])  
    # Binomial
    Sb <- crossprod(X[[3]], Y[[3]] - Eexpeta.overetas[[1]])
    c(Sg, Sp, Sb)
  }, X = dmats$X, Y = dmats$Y, Z = dmats$Z, W = dmats$W, Eb = Eb,
  Eeta = Eeta, Eexpeta = Eexpeta, Eexpeta.overetas = Eexpeta.overetas)
  
  # Dispersion ('\sigma') =====================
  sig <- sigma[[1]]
  Ssig <- Map(function(a, mi){
    -mi[[1]]/(2 * sig) + a/(2 * sig^2)
  }, a = a, mi = dmats$mi)
  
  # Survival parameters (\gamma, \zeta) =======
  Psurv <- length(c(gamma, zeta))
  
  Sgz <- Map(function(draws, Eb, S, SS, Fu, Fi, l0u, Delta, ST){
    if(length(ST))
      return(Sgammazeta(c(gamma, zeta), draws, Eb, S, SS, Fu, Fi, l0u, Delta,
                        inds$Cpp$b, dmats$K, .Machine$double.eps^(1/3)))
    else
      return(rep(0, Psurv))
  }, draws = draws, Eb = Eb, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, Sb, Ss, Sgz)
  }, sD = sD, Sb = Sb, Ss = Ssig, Sgz = Sgz)
  
  SS <- rowSums(S) # sum 
  #  observed empirical information matrix (Mclachlan and Krishnan, 2008).
  SiSiT <- Reduce('+', lapply(1:dmats$n, function(i) tcrossprod(S[, i])))
  H <- SiSiT - tcrossprod(SS)/dmats$n
  
  return(list(Score = S,
              Hessian = H,
              b.hat = b.hat,
              Sigma = Sigma))
}

#' Extract the variance-covariance matrix from a \code{joint} fit.
#' 
#' @details Uses the observed-empirical \strong{approximation} of information matrix 
#' (Mclachlan & Krishnan, 2008). The standard errors for the baseline hazard are not estimated. 
#' 
#' @param object a joint model fit by the \code{joint} function.
#' @param corr should the correlation matrix be returned instead of the variance-covariance?
#' @param ... extra arguments (none used).
#' 
#' @return A variance-covariance matrix for the joint model object.
#'
#' @author James Murray \email{j.murray7@@ncl.ac.uk}
#' 
#' @section Methodology: 
#' 
#' Many competing ways exist for obtaining the observed information matrix in an EM algorithm. 
#' In the context of joint modelling, the observed empirical approximation of the information 
#' matrix has been used previously (\code{joineRML}, Hickey et al. 2018). Elsewhere,
#' estimation of the observed information in a semi-parametric setting is outlined neatly in
#' Xu et al. (2014). Here, they advocate for approximation of this information matrix by 
#' numerical differentiation of the profile Fisher Score vector. We do not consider this 
#' methodology owing to its computational expense. That is, for each element of \eqn{\Omega} 
#' which is perturbed by some small amount \eqn{\tilde{\Omega}^{p}}, we must re-calculate
#' \eqn{\hat{b}_i} and \eqn{\hat{\Sigma}_i}.
#' 
#' @references 
#' 
#' Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. \code{joineRML}: a joint model and
#' software package for time-to-event and multivariate longitudinal outcomes.
#' \emph{BMC Med. Res. Methodol.} 2018; \strong{50}
#' 
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second Edition. 
#' Wiley-Interscience; 2008.
#' 
#' Xu C, Baines PD, Wang J. Standard error estimation using the EM algorithm for the joint 
#' modeling of survival and longitudinal data. \emph{Biostatistics} 2014; \strong{15}(4).
#' 
#' @method vcov joint
#' @export
#' 
#' @examples
#' # Univariate fit on PBC data -------------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify univariate fit
#' long.formulas <- list(
#'   albumin ~ time + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian'))
#' 
#' vcov(fit)
vcov.joint <- function(object, corr = FALSE, ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  v <- object$vcov
  if(corr) 
    return(cov2cor(v))
  else
    return(v)
}

