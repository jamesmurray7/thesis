# For density of exp(eta)|b;bO (i.e. showing it's log normal).
# Below is copied from sim.R in same directory.

source('../theme_csda.R')
library(gmvjoint)
theta <- c(-1, 0.0)
.sim <- function(family, n = 100, D = NULL){
  if(family != "binomial"){
    D <- if(is.null(D)) diag(c(0.25, 0.09)) else D
    random.formula <- NULL
  }else{
    D <- if(is.null(D)) matrix(2, 1, 1) else D
    random.formula <- list(~1)
  }
  
  a <- gmvjoint::simData(n = n, ntms = 15, theta = theta,
                         beta = t(c(2, -0.1, 0.1, -0.2)),
                         sigma = c(0.16),
                         D =  D,
                         zeta = c(0, -0.2),
                         family = as.list(family),
                         random.formula = random.formula,
                         gamma = 0.5,
                         return.ranefs = TRUE)
  list(a$data, a$ranefs, D)
}

# Obtain b.hat and Sigma.hat given observed data at TRUE parameter estimates.
getSigma <- function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D){
  ln <- nrow(X[[1]])
  uu <- optim(b, gmvjoint:::joint_density, gmvjoint:::joint_density_ddb,
              Y = Y, X = X, Z = Z, W = matrix(1, nr = ln, nc = 1), 
              beta = c(2, -0.1, 0.1, -0.2), D = D, sigma = list(0.16),
              family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
              gamma_rep = rep(0.5, ncol(Z[[1]])), zeta = -0.2, beta_inds = list(0:3), b_inds = b.inds,
              K = 1L, method = 'BFGS', hessian = T)
  list(bhat = uu$par, Sigma = solve(uu$hessian))
}

# Function to sample given data, true random effects, known family and target ids.
Sample <- function(data, btrue, family, ids, D){
  # Longit.
  X <- Y <- Z <- setNames(vector('list', length(ids)), paste0("id: ", ids))
  for(i in seq_along(ids)){
    X[[i]] <- Y[[i]] <- Z[[i]] <- list()
    X[[i]][[1]] <- model.matrix(~time+cont+bin, data[data$id==ids[i],,drop=F])
    if(family!='binomial') 
      Z[[i]][[1]] <- model.matrix(~time, data[data$id==ids[i],,drop=F])
    else
      Z[[i]][[1]] <- model.matrix(~1, data[data$id==ids[i],,drop=F])
    Y[[i]][[1]] <- data[data$id==ids[i],'Y.1']
  }
  b <- lapply(seq_along(ids), function(x) btrue[ids[x],,drop=F])
  # Survival part
  fts <- sort(unique(data[data$status==1,'survtime']))
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data)
  l0 <- exp(theta[1] + theta[2] * fts)
  if(family == "binomial"){
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1|id)), gmvjoint:::parseFormula), l0)
    b.inds <- list(0)
    gamma.rep <- 0.5
  }else if(family == "gaussian"){
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), gmvjoint:::parseFormula), l0)  
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }else{
    sv <- gmvjoint:::surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), gmvjoint:::parseFormula), l0)  
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }
  
  Omega <- list(D = D, beta = c(2, -0.1, 0.1, -0.2), sigma = list(0.16),
                gamma = 0.5, zeta = -0.2)
  
  # Tuning parameters, given D, these return about 22.5% acceptane.
  tune <- if(family == 'gaussian') 6 else if(family == 'poisson') 6 else 23
  
  Delta <- lapply(seq_along(ids), function(x) surv$Delta[[ids[x]]])
  S <- lapply(seq_along(ids), function(x) sv$S[[ids[x]]])
  Fi <- lapply(seq_along(ids), function(x) sv$Fi[[ids[x]]])
  l0i <- lapply(seq_along(ids), function(x) sv$l0i[[ids[x]]])
  SS <- lapply(seq_along(ids), function(x) sv$SS[[ids[x]]])
  Fu <- lapply(seq_along(ids), function(x) sv$Fu[[ids[x]]])
  l0u <- lapply(seq_along(ids), function(x) sv$l0u[[ids[x]]])
  
  Sigma <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    Sigma <- getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, 
  Fu = Fu, l0u = l0u)
  
  b.hat <- lapply(Sigma, el, 1)
  Sigma <- lapply(Sigma, el, 2)
  
  # 5 from each group of follow-up length
  mi <- sapply(Z, function(x) nrow(x[[1]]))
  groups <- cut(mi, c(1, 5, 10, 15), include.lowest = T)
  ids.to.sample <- c(sapply(levels(groups), function(i) sample(which(groups==i), 5)))
  groups <- as.list(as.character(groups[ids.to.sample]))
  
  cli::cli_progress_bar(name = "Sampling...", total = length(ids))
  b.out <- vector("list", length(ids.to.sample))
  Acc <- numeric(length(ids.to.sample))
  for(i in seq_along(ids.to.sample)){
    a <- ids.to.sample[i]
    W <- gmvjoint:::metropolis(b[[a]], Omega, Y[[a]], X[[a]], Z[[a]], list(matrix(1, nr = nrow(X[[a]][[1]]), nc = 1)),
                               list(family), Delta[[a]], S[[a]], Fi[[a]], l0i[[a]], 
                               SS[[a]], Fu[[a]], l0u[[a]], gamma.rep,
                               list(0:3), b.inds, 1L, length(b.inds[[1]]), 500, 20000, Sigma[[a]], tune)
    Acc[i] <- W$AcceptanceRate
    b.out[[i]] <- t(W$walks)
    rm(W)
    cli::cli_progress_update(inc = 1L)
  }
  cli::cli_progress_done()
  
  linps <- lapply(seq_along(ids.to.sample), function(i){
    id <- ids.to.sample[i]
    mi <- groups[[i]]
    b.out.i <- b.out[[i]]
    hat.b.i <- b.hat[[id]]; Sigma.i <- Sigma[[id]]; 
    X.i.beta <- X[[id]][[1]] %*% c(2, -0.1, 0.1, -0.2)
    Z.i <- Z[[id]][[1]]
    # The linear predictors
    linps <- apply(b.out.i, 1, function(x) c(X.i.beta + Z.i %*% c(x)))
    if(!"matrix"%in%class(linps)) linps <- t(linps)
    # Expected quantities
    Med <- exp(X.i.beta + Z.i %*% c(hat.b.i))
    A <- tcrossprod(Z.i %*% Sigma.i, Z.i)
    tau2 <- diag(A)
    Mean <- exp(X.i.beta + Z.i %*% c(hat.b.i) + tau2/2)
    list(id = id, migroup = mi, 
         linps = linps, Median = Med, Mean = Mean)
  })
  
  return(linps)
}

# Show t = {0,5,10,15} terms only?
getLinps <- function(family, n = 100, D = NULL){ # Wrapper for simulation + Sampling
  d <- .sim(family, n, D)
  btrue <- d[[2]]; D <- d[[3]]; data <- d[[1]]
  Sample(data, btrue, family, 1:n, D)
}

# get closest (for arrow drawing)
fc <- function(point, x) which.min(abs(point-x$x))
# function to parse + plot all
parseLinps <- function(x, times){
  grps <- sapply(x, '[[', "migroup") # Randomly sample ids from these later.
  ids <- sapply(x, '[[', "id")       # and their "actual" ids
  t.surv <- sapply(x, function(x) length(x$Med)) # t's "survived".
  ref <- data.frame(id = ids, grp = grps, l = t.surv)
  
  for(t in times){
    cand <- ref[ref$l >= t,]  
    # Now want to select _one_ profile to show for each 
    unq.grp <- unique(cand$grp)
    cand.t <- do.call(rbind, lapply(unq.grp, function(y){
      cand.grp <- cand[cand$grp==y,]
      samp <- sample(1:nrow(cand.grp), 1, FALSE)
      as.data.frame(cand.grp[samp,])
    }))
    # Wittled down, now which id's are we choosing?
    inds.t <- which(ids%in%cand.t$id)
    
    step <- lapply(seq_along(inds.t), function(i){
      ii <- inds.t[i]
      linp <- exp(x[[ii]]$linps[(t+1),])
      Med <- x[[ii]]$Med[(t+1)]
      Mean <- x[[ii]]$Mean[(t+1)]
      this <- cand.t[cand.t$id==cand.t$id[i],]
      list(linp = linp, Med = Med, Mean = Mean, this = this)
    })
    
    # This is all for same plot, so work out limits.
    Dens <- sapply(step, function(x) density(x$linp))
    xlims <- sapply(Dens["x",,drop=F], range); ylims <- sapply(Dens["y",,drop=F], range)
    xlim <- c(min(xlims[1,]), max(xlims[2,]))
    ylim <- c(min(ylims[1,]), max(ylims[2,]))
    for(p in 1:ncol(Dens)){
      grp.p <- step[[p]]$this$grp
      p.col <- ifelse(grp.p=="[1,5]", "black",
                      ifelse(grp.p=="(5,10]", "blue2", "red2"))
      # plot/draw line
      if(p == 1){
        plot(x = Dens["x", p]$x, y = Dens["y", p]$y, type = "l",
             col = p.col, xlab = expression(exp(eta)), ylab = "Density",
             xlim = xlim, ylim = ylim, main = bquote(t[.(t)]))
        # Draw median line
        arrows(x0 = step[[p]]$Med, x1 = step[[p]]$Med,
               y0 = 0, y1 = Dens["y", p]$y[fc(step[[p]]$Med, Dens["x", p])],
               angle = 0, col = p.col, lty = 3)
        # Draw mean line
        arrows(x0 = step[[p]]$Mean, x1 = step[[p]]$Mean,
               y0 = 0, y1 = Dens["y", p]$y[fc(step[[p]]$Mean, Dens["x", p])],
               angle = 0, col = p.col, lty = 5)
      }else{
        lines(x = Dens["x", p]$x, y = Dens["y", p]$y, col = p.col)
        # Draw median line
        arrows(x0 = step[[p]]$Med, x1 = step[[p]]$Med,
               y0 = 0, y1 = Dens["y", p]$y[fc(step[[p]]$Med, Dens["x", p])],
               angle = 0, col = p.col, lty = 3)
        # Draw mean line
        arrows(x0 = step[[p]]$Mean, x1 = step[[p]]$Mean,
               y0 = 0, y1 = Dens["y", p]$y[fc(step[[p]]$Mean, Dens["x", p])],
               angle = 0, col = p.col, lty = 5)
      }
    }
  }
}

parseLinps(P, 3)
