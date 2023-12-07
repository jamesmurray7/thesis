#' #################################
#' Shape of integrands             #
#' -------------------             #
#' 2. Exponent in e.g. poisson     #
#'    linear predictor             #
#' #################################

rm(list=ls())
source(".Rprofile")
# Generate `f` data
create.f.plot <- function(f, show.mode = F){
  file.name <- save.dir.file.path(paste0(f, "_ExpLinp.png"))
  if(show.mode) file.name <- gsub("\\.png", "_with_mode.png", file.name)
  if(f != "genpois"){
    dd <- simData(n = 500, zeta = c(0,-0.2), sigma = list(1.0),
                  beta = c(2.0,-0.1,0.1,0.2),
                  D = matrix(c(.5, .125, .125, .09), 2, 2),
                  gamma = .5, family = list(f))$data
  }else{
    dd <- simData(n = 500, zeta = c(0,-0.2), sigma = list(0.2),
                  beta = c(0.5, -0.2, 0.05, 0.4),
                  gamma = 0.5, family = list(f),
                  D = matrix(c(0.5, .125, .125, .05), 2, 2),theta = c(-2.95, .1))$data
  }
  
  cli::cli_alert("Creating {f} E[exp(eta)|...] plots\n")
  cli::cli_alert("which will be saved in {file.name}.\n")
  
  jj <- joint(
    list(
      Y.1 ~ time + cont + bin + (1 + time|id)
    ),
    Surv(survtime, status) ~ bin,
    dd, list(f)
  )
  
  # Take first `id` with a full profile
  id <- min(which(with(dd, tapply(time, id, length))==10))
  
  L <- jj$dmats$long
  S <- jj$dmats$surv
  M <- jj$ModelInfo
  Samp <- Metropolis_Hastings(
    b = jj$REs[id,], Omega = jj$coeffs, Y = L$Y[[id]], X = L$X[[id]], Z = L$Z[[id]],
    W = L$W[[id]], family = jj$ModelInfo$family, Delta = jj$dmats$ph$Delta[[id]],
    S = S$S[[id]], Fi = S$Fi[[id]], l0i = S$l0i[[id]], SS = S$SS[[id]], Fu = S$Fu[[id]],
    haz = S$l0u[[id]], gamma_rep = rep(jj$coeffs$gamma, sapply(M$inds$R$b, length)), 
    beta_inds = M$inds$Cpp$beta, b_inds = M$inds$Cpp$b, K = M$K, q = M$Pcounts$q,
    burnin = 1000L, N = 50000L, Sigma = gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[id,], M$Pcounts$q),
    b_dist = "t", df = 4L, tune = 3.5
  )
  
  cli::cli_alert_warning("Sampling complete! Acceptance rate: {Samp$Acc}\n")
  
  W <- t(Samp$walks)
  
  # Shape of E[\exp{SS^T\zeta + \sum_k\gamma_k\b{_k)}]
  b <- jj$REs[id,]
  Sig <- gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[id,], M$Pcounts$q)
  tau <- sqrt(diag(tcrossprod(L$Z[[id]][[1]] %*% Sig, L$Z[[id]][[1]])))
  I <- apply(W,1,function(x){
    exp(L$X[[id]][[1]]%*%jj$coeffs$beta + L$Z[[id]][[1]] %*% x)
    # exp(L$Z[[id]][[1]] %*% x)
  })
  
  dens <- apply(I, 1, density)
  med <- exp(L$X[[id]][[1]]%*%jj$coeffs$beta + L$Z[[id]][[1]] %*% b)
  mode <- med * exp(-tau^2)
  # med <- exp(L$Z[[id]][[1]] %*% b)
  mean <- med * exp(tau^2/2)
  
  plot.dens.with.point <- function(i){
    di <- dens[[i]]
    x <- di$x; y <- di$y
    plot(x, y ,'l', 
         xlab = expression(exp(X[i]*beta*"+"*Z[i]*b[i])),
         ylab = 'Density')
    abline(v = med[i], col = 'blue', lty = 'dotted')
    abline(v = mean[i], col = 'red', lty = 'dotted')
    if(show.mode) abline(v = mode[i], col = 'darkgreen', lty = 'dotted')
  }
  
  png(file.name,
      width = 140, height = 160, units = 'mm', res = 5e2)
  par(mfrow = c(3,1))
  for(pp in c(1, 5, 10)) # start, midway, end
    plot.dens.with.point(pp)
  dev.off()
  
  cli::cli_alert_success("Done for {f}!\n\n\n")
  
}

create.f.plot("poisson")
create.f.plot("negbin")
create.f.plot("Gamma")
create.f.plot("binomial")
create.f.plot("genpois")
