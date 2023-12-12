# Obtain b.hat and Sigma.hat given observed data at TRUE parameter estimates.
getSigma <- function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, Omega, include.survival){
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  gamma.rep <- Omega$gamma.rep
  ln <- nrow(X[[1]])
  
  if(!include.survival){ # `Remove' RHS of complete data likelihood...
    Delta <- Delta * 0   # i.e. find tilde{b}_i, tilde{Sigma}_i
    l0u <- l0u * 0
    l0i <- l0i * 0
    zeta <- zeta * 0
    gamma <- gamma * 0
    gamma.rep <- gamma.rep * 0
  }
  
  uu <- optim(b, gmvjoint:::joint_density, gmvjoint:::joint_density_ddb,
              Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma,
              family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
              gamma_rep = gamma.rep, zeta = zeta, beta_inds = list(0:3), b_inds = b.inds, K = 1L, 
              method = 'BFGS', hessian = T)
  list(bhat = uu$par, Sigma = solve(uu$hessian))
}

# This is a rip from `Sample.R` but only keeps b.hat and Sigma.hat (i.e. no sampling)
get.b <- function(X_, include.survival = TRUE, force.intslope = TRUE){
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
    Sigma <- getSigma(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds.cpp, Omega, include.survival)
  }, b = b, Y = Y, X = X, Z = Z, W = W, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, l0u = l0u)
  
  b.hat <- lapply(Sigma, el, 1)
  Sigma <- lapply(Sigma, el, 2)
  mi <- sapply(ids, function(x) nrow(Z[[x]][[1]]))
  structure(list(b.hat = b.hat, Sigma.hat = Sigma, mi = mi, surv = include.survival),
            class = 'Napprox')
}

compare.N.approx <- function(X_surv, X_nosurv){
  stopifnot(inherits(X_surv, "Napprox") && inherits(X_nosurv, "Napprox"))
  if(isFALSE(X_surv$surv)) stop("First argument should be survival one")
  if(isTRUE(X_nosurv$surv)) stop("First argument should be survival one")
  
  b.hats1 <- X_surv$b.hat; b.hats2 <- X_nosurv$b.hat
  Sigma.hats1 <- X_surv$Sigma.hat
  # Ellipse for survival approximations
  eig.S1 <- lapply(Sigma.hats1, eigen)
  theta.S1 <- sapply(eig.S1, function(x){
    xx <- atan2(x$vec[2,1], x$vectors[1,1])
    if(xx < 0) xx <- xx + 2 * pi
    xx
  })
  semi.maj1 <- sapply(eig.S1, function(x) sqrt(x$values[1] * qchisq(0.95, 2)))
  semi.min1 <- sapply(eig.S1, function(x) sqrt(x$values[2] * qchisq(0.95, 2)))
  
  # Ellipse for non-survival approximations
  Sigma.hats2 <- X_nosurv$Sigma.hat
  eig.S2 <- lapply(Sigma.hats2, eigen)
  theta.S2 <- sapply(eig.S2, function(x){
    xx <- atan2(x$vec[2,1], x$vectors[1,1])
    if(xx < 0) xx <- xx + 2 * pi
    xx
  })
  semi.maj2 <- sapply(eig.S2, function(x) sqrt(x$values[1] * qchisq(0.95, 2)))
  semi.min2 <- sapply(eig.S2, function(x) sqrt(x$values[2] * qchisq(0.95, 2)))
  
  # Differences between the two approximations {1}-{2}
  loc.diff <- do.call(rbind, Map(function(b1,b2) b1 - b2, b1 = b.hats1, b2 = b.hats2))
  colnames(loc.diff) <- c("b[0]", "b[1]")
  ang.diff <- theta.S1 - theta.S2
  maj.diff <- semi.maj1 - semi.maj2
  min.diff <- semi.min1 - semi.min2
  
  if(!all(X_surv$mi == X_nosurv$mi)) stop("Two different samples?")
  
  structure(list(mi = X_surv$mi,
                 loc.diff = loc.diff,
                 ang.diff = ang.diff,
                 maj.diff = maj.diff,
                 min.diff = min.diff,
                 which.surv = "X_surv"),
            class = "dist.compare.new")
}

# S3 ----------------------------------------------------------------------
print.dist.compare.new <- function(x){
  stopifnot(inherits(x, "dist.compare.new"))
  cat(sprintf("\n----\n%s is the survival density\n----\n", x$which.surv))
  cat("\nSummary of b.hat differences:\n ")
  print(summary(x$loc.diff))
  cat("\nSummary of angle differences:\n")
  print(summary(x$ang.diff))
  cat("\nSummary of semi-major axes differences:\n")
  print(summary(x$maj.diff))
  cat("\nSummary of semi-minor axes differences:\n")
  print(summary(x$min.diff))
}

plot.dist.compare.new <- function(x, show.angles = FALSE){
  stopifnot(inherits(x, "dist.compare.new"))
  loc.df <- as.data.frame(x$loc.diff)
  loc.df$mi <- x$mi
  P1 <- ggplot(loc.df, aes(x = `b[0]`, y = `b[1]`, colour = mi)) + 
    geom_point(size=.10) + 
    labs(x = expression(hat(b)[0]-tilde(b)[0]),
         y = expression(hat(b)[1]-tilde(b)[1]),
         colour = expression(m[i])) + 
    scale_colour_gradient(low = .nice.orange, high = "red3",
                          breaks = c(2,4,6,8,10)) + 
    theme_csda()+
    theme(
      legend.position = "none",
      axis.text = element_text(size = 5)
    )
  
  r.df <- data.frame(mi = x$mi, ry = x$min.diff, rx = x$maj.diff)
  P2 <- ggplot(r.df, aes(x = rx, y = ry, colour = mi)) + 
    geom_point(size=.10) + 
    labs(x = expression(r[x]^{"(S)"}-r[x]^{"(NS)"}),
         y = expression(r[y]^{"(S)"}-r[y]^{"(NS)"}),
         colour = expression(m[i])) + 
    scale_color_gradient(low = .nice.orange, high = "red3",
                         breaks = c(2,4,6,8,10)) + 
    theme_csda() + 
    theme(
      legend.key.width = unit(3, "mm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      axis.text = element_text(size = 5)
    )
  
  ang.df <- data.frame(mi = x$mi, theta = x$ang.diff)
  P3 <- ggplot(ang.df, aes(x = mi, y = theta, fill = mi, group = mi)) + 
    geom_hline(yintercept = 0, lty = 5, alpha = 0.33) +
    geom_boxplot(outlier.alpha = 0.33, lwd = 0.25, 
                 fatten = 2, outlier.size = 0.50) + 
    labs(x = expression(m[i]),
         y = expression(vartheta^{"(S)"}-vartheta^{"NS"})) + 
    scale_fill_gradient(low = .nice.orange, high = "red3") + 
    scale_x_continuous(breaks = 1:10) + 
    theme_csda() + 
    theme(
      legend.position = 'none',
    )
  
  if(show.angles)
    gridExtra::grid.arrange(P1, P2, P3,
                            layout_matrix = rbind(
                              c(1,1,2,2),
                              c(1,1,2,2),
                              c(3,3,3,3),
                              c(3,3,3,3)
                            ))
  else
    gridExtra::grid.arrange(P1, P2,
                            layout_matrix = rbind(
                              c(1,1,2,2),
                              c(1,1,2,2),
                              c(1,1,2,2),
                              c(1,1,2,2)
                            ))
}

