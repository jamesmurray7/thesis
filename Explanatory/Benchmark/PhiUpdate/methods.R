# Way 1: Pure C++ implementation
# (as appears in gmvjoint)
gmvj <- function(){ 
  
  S.gammazeta <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    gmvjoint:::Sgammazeta(c(Omega$gamma, Omega$zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 
               Info$inds$Cpp$b, Info$K, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  H.gammazeta <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    gmvjoint:::Hgammazeta(c(Omega$gamma, Omega$zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 
                          Info$inds$Cpp$b, Info$K, .Machine$double.eps^(1/4))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  c(Omega$gamma, Omega$zeta) - solve(Reduce('+', H.gammazeta), rowSums(S.gammazeta))
  
}

# Way 2: R call onto C++ implementation
R.on.Cpp <- function(){
  
  S.gammazeta <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::grad(gmvjoint:::Egammazeta, c(Omega$gamma, Omega$zeta), heps = .Machine$double.eps^(1/3),
                 b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                 w = w, v = v, inds = Info$inds$Cpp$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
     Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)

  H.gammazeta <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::hessian(gmvjoint:::Egammazeta, c(Omega$gamma, Omega$zeta), h = .Machine$double.eps^(1/4),
                 b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                 w = w, v = v, inds = Info$inds$Cpp$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
     Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  c(Omega$gamma, Omega$zeta) - solve(Reduce('+', H.gammazeta), rowSums(S.gammazeta))
  
}

# Way 3: Defining conditional expectation in R (should be a lot slower)
# Exclusively R --> Mimicking C++ implementation _exactly_
Egammazeta.R <- function(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, 
                         w, v, inds, K){
  q <- ncol(Fu); gh <- length(w);
  gam <- gammazeta[1:K]; zet <- gammazeta[(K+1):length(gammazeta)]
  store <- numeric(gh); gammas <- numeric(q)
  for(k in 1:K){
    gk <- gam[k]
    gammas[inds[[k]]] <- gk
  }
  gmat <- diag(gammas)
  Q <- Fu %*% gmat
  A <- tcrossprod(Q %*% Sigma, Q)
  mu <- SS %*% zet + Q %*% b
  tau <- sqrt(diag(A))
  for(l in 1:gh){
    store[l] <- w[l] * crossprod(haz, exp(mu + tau * v[l]))
  }
  Delta * (S %*% zet + Fi %*% (b * gammas)) - sum(store)
}

# Exclusively R --> Avoiding loops (should be a little quicker?)
Egammazeta.R.noloop <- function(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, 
                                w, v, inds, K){
  q <- ncol(Fu); gh <- length(w);
  gam <- gammazeta[1:K]; zet <- gammazeta[(K+1):length(gammazeta)]
  gammas <- do.call(c, lapply(1:K, function(k){
    rep(gam[k], length(inds[[k]]))
  }))
  gmat <- diag(gammas)
  Q <- Fu %*% gmat
  A <- tcrossprod(Q %*% Sigma, Q)
  mu <- SS %*% zet + Q %*% b
  tau <- sqrt(diag(A))
  store <- sapply(1:gh, function(l){
    w[l] * crossprod(haz, exp(mu + tau * v[l]))
  })
  Delta * (S %*% zet + Fi %*% (b * gammas)) - sum(store)
}

# And define the two wrapper functions
R.way1 <- function(){
  S.gammazeta <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::grad(Egammazeta.R, c(Omega$gamma, Omega$zeta), heps = .Machine$double.eps^(1/3),
                 b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                 w = w, v = v, inds = Info$inds$R$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  H.gammazeta <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::hessian(Egammazeta.R, c(Omega$gamma, Omega$zeta), h = .Machine$double.eps^(1/4),
                    b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                    w = w, v = v, inds = Info$inds$R$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  c(Omega$gamma, Omega$zeta) - solve(Reduce('+', H.gammazeta), rowSums(S.gammazeta))
}

R.way2 <- function(){
  S.gammazeta <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::grad(Egammazeta.R.noloop, c(Omega$gamma, Omega$zeta), heps = .Machine$double.eps^(1/3),
                 b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                 w = w, v = v, inds = Info$inds$R$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  H.gammazeta <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    pracma::hessian(Egammazeta.R.noloop, c(Omega$gamma, Omega$zeta), h = .Machine$double.eps^(1/4),
                    b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, haz = l0u, Delta = Delta,
                    w = w, v = v, inds = Info$inds$R$b, K = Info$K)
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = sv$l0u, Delta = surv$Delta)
  
  c(Omega$gamma, Omega$zeta) - solve(Reduce('+', H.gammazeta), rowSums(S.gammazeta))
}



