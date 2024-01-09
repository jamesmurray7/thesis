# Fit a joint model to time-to-event and multivariate longitudinal data
# `joint` is untouched, but EMUpdate EXPECTS gaussian/poisson/binomial response in that order.

joint <- function(long.formulas, surv.formula, 
                  data, family,
                  disp.formulas = NULL, 
                  MCtype = 'ordinary', N = 5e3,
                  control = list()){
  
  con <- list(correlated = T, gh.nodes = 3, gh.sigma = 1, center.ph = T,
              tol.abs = 1e-3, tol.rel = 1e-2, tol.thr = 1e-1, tol.den = 1e-3,
              grad.eps = .Machine$double.eps^(1/3), hess.eps = .Machine$double.eps^(1/4), inits = NULL,
              maxit = 200, conv = 'sas', verbose = F, return.inits = F, return.dmats = T, post.process = T)
  conname <- names(con)
  if(any(!names(control)%in%conname)){
    warning("Supplied control arguments do not match with possible names:\n", 
            paste0(sapply(conname, sQuote), collapse=', '), '.')
  }
  con[(conname <- names(control))] <- control
  
  if(!MCtype%in%c("ordinary", "antithetic", "sobol")) stop("Invalid MCtype")
  
  # Ensure supplied families are character, not functions
  if(!is.list(family) | !all(sapply(family, function(x) is.character(x) & !is.function(x))))
    stop(sQuote("family"), " must be supplied as a list of character strings")
  
  start.time <- proc.time()[3]
  
  # Initial parsing ----
  pf <- parent.frame()
  if("factor"%in%class(data$id)) data$id <- as.numeric(as.character(data$id))
  formulas <- lapply(long.formulas, parseFormula)
  surv <- parseCoxph(surv.formula, data, con$center.ph)
  n <- surv$n; K <- length(family)
  if(K!=length(long.formulas)) stop('Mismatched lengths of ', sQuote("family"), " and ", sQuote("long.formulas"),'.')
  # Parse dispersion formulas
  if(is.null(disp.formulas)){
    disp.formulas <- replicate(K, ~1, simplify = FALSE)
  }else{
    if(sum(!sapply(disp.formulas, is.null)) != K)
      stop("Need to supply dispersion formulas for all responses even if not required for all K responses\n",
           "(You can just fill-in ", sQuote("~1"), " for those with no wanted dispersion model.)")
  }
  # Reduce memory overheads, particularly for returned object(s) ?
  disp.formulas <- lapply(disp.formulas, function(x){
    environment(x) <- pf
    x
  })
  
  # Initial conditons ----
  inits.long <- Longit.inits(long.formulas, disp.formulas, data, family)
  id.assign <- AssignIds(data)
  dmats <- getdmats(inits.long$fits, id.assign)
  # Suss out indices of b_k and beta_k.
  b.inds <- lapply(seq_along(dmats$q), function(x){
    xx <- seq(dmats$q[x])
    if(x > 1) 
      return(xx + sum(dmats$q[1:(x-1)]))
    else
      return(xx)
  })
  beta.inds <- lapply(seq_along(dmats$P), function(x){
    xx <- seq(dmats$P[x])
    if(x > 1) 
      return(xx + sum(dmats$P[1:(x-1)]))
    else
      return(xx)
  })
  inds <- list(
    R = list(b = b.inds, beta = beta.inds),
    Cpp = list(b = lapply(b.inds, function(x) x - 1), beta = lapply(beta.inds, function(x) x - 1))
  )
  
  inits.surv <- TimeVarCox(data, inits.long$b, surv, formulas, b.inds, inits.long)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  names.beta <- names(beta)
  D <- inits.long$D.init
  sigma <- inits.long$sigma.init # dispersions
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  # zeta <- inits.surv$inits[match(names(surv$ph$assign), names(inits.surv$inits))]
  zeta <- inits.surv$inits[1:surv$pS]
  names(zeta) <- paste0('zeta_', names(zeta))
  # gamma <- inits.surv$inits[grepl('gamma\\_', names(inits.surv$inits))]
  gamma <- inits.surv$inits[(surv$pS+1):length(inits.surv$inits)]
  
  # Survival data objects 
  sv <- surv.mod(surv, formulas, inits.surv$l0.init, inits.long)
  
  # Parameter vector and list ----
  Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, unlist(sigma)[inits.long$sigma.include], gamma, zeta)
  sigma.include <- inits.long$sigma.include
  if(!is.null(con$inits)){
    parsed.inits <- parseInits(con$inits, params, inds, inits.long, Omega)
    params <- parsed.inits$params
    Omega <- parsed.inits$Omega
  }
  if(!con$return.inits) rm(inits.surv)
  
  # step-sizes for grad/hess in dispersion updates
  if(length(con$grad.eps) == 1L & !is.list(con$grad.eps)) con$grad.eps <- replicate(K, con$grad.eps, simplify = FALSE)
  if(length(con$hess.eps) == 1L & !is.list(con$hess.eps)) con$hess.eps <- replicate(K, con$hess.eps, simplify = FALSE)
  if(length(con$grad.eps) != K) stop("Control argument ", sQuote("grad.eps"), " is not of appropriate length.\n")
  if(length(con$hess.eps) != K) stop("Control argument ", sQuote("hess.eps"), " is not of appropriate length.\n")
  
  # Gauss-Hermite Quadrature ----
  GH <- statmod::gauss.quad.prob(con$gh.nodes, 'normal', sigma = con$gh.sigma)
  w <- GH$w; v <- GH$n

  # Begin EM ----
  diff <- 100; iter <- 0;
  # Convergence criteria setup
  if(!con$conv%in%c('abs', 'rel', 'either', 'sas')){
    warning("Convergence criteria must be one of ", sQuote('abs'), ', ', sQuote('rel'), ', ',
            sQuote('either'), ', or ', sQuote('sas'), '. Using ', sQuote('sas'), '.')
    con$conv <- "sas"
  }
    
  convergence.criteria <- list(type = con$conv, tol.abs = con$tol.abs, tol.rel = con$tol.rel, tol.den = con$tol.den,
                               threshold = con$tol.thr)
  
  if(con$verbose){
    cat("Initial conditions: \n")
    converge.check(0, 0, convergence.criteria, 0, Omega, TRUE)
    cat("\nStarting EM algorithm...\n")
  }
  converged <- FALSE
  XTX <- lapply(dmats$X, function(x) crossprod(x[[1]]))
  XTX <- Reduce("+", XTX)
  EMstart <- proc.time()[3]
  while((!converged) && (iter < con$maxit)){
    update <- EMupdate(Omega, family, dmats, b, sv, 
                       surv, MCtype, N, con, inds, XTX)
    
    if(!con$correlated) update$D[inits.long$off.inds] <- 0
    params.new <- c(vech(update$D), update$beta, unlist(update$sigma)[sigma.include], 
                    update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    # Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; sigma <- update$sigma
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    sv$l0 <- l0; sv$l0i <- l0i; sv$l0u <- l0u
    iter <- iter + 1
    Omega <- list(D = D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
    
    convcheck <- converge.check(params, params.new, convergence.criteria, iter, Omega, con$verbose)
    if(iter >= 4) converged <- convcheck$converged # Allow to converge after 3 iterations
    params <- params.new
  }
  
  EMend <- proc.time()[3]
  coeffs <- Omega
  coeffs$beta <- setNames(c(Omega$beta), names.beta)
  if(is.not.SPD(coeffs$D)) warning("Covariance matrix D is not positive semi-definite at convergence,",
                                   " potential model misspecification or lower tolerance options required.")
  
  out <- list(coeffs = coeffs,
              hazard = cbind(ft = sv$ft, haz = l0, nev = sv$nev))
  
  ModelInfo <- list()
  ModelInfo$ResponseInfo <- sapply(1:K, function(k){
    paste0(inits.long$responses[k], ' (', family[k], ')')
  })
  ModelInfo$Resps <- inits.long$responses
  ModelInfo$family <- family
  ModelInfo$K <- dmats$K
  ModelInfo$Pcounts <- list(P = dmats$P, Pd = dmats$Pd, 
                            q = sv$q, vD = length(vech(D)))
  ModelInfo$long.formulas <- long.formulas
  ModelInfo$disp.formulas <- disp.formulas
  ModelInfo$surv.formula <- surv.formula
  ModelInfo$survtime <- surv$survtime
  ModelInfo$status <- surv$status
  ModelInfo$control <- con
  ModelInfo$convergence.criteria <- convergence.criteria
  ModelInfo$inds <- inds
  ModelInfo$n <- dmats$n
  ModelInfo$nobs <- setNames(dmats$m, inits.long$responses)
  ModelInfo$mi <- sapply(dmats$mi, unlist)
  ModelInfo$nev <- sum(sv$nev)
  ModelInfo$MCtype <- MCtype
  ModelInfo$N <- N
  ModelInfo$id.assign <- list(original.id = id.assign$id,
                              assigned.id = id.assign$assign)
  out$ModelInfo <- ModelInfo
  
  # Post processing ----
  if(con$post.process){
    if(con$verbose) cat('Post-processing...\n')
    gamma.rep <- rep(gamma, sapply(b.inds, length))
    pp.start.time <- proc.time()[3]
    
    II <- obs.emp.I(coeffs, dmats, surv, sv, family, b, 
                    l0i, l0u, MCtype, N, inds, con)
    H <- structure(II$Hessian,
                   dimnames = list(names(params), names(params)))
    
    I <- tryCatch(solve(H), error = function(e) e)
    if(inherits(I, 'error')) I <- structure(MASS::ginv(H), dimnames = dimnames(H))
    out$Hessian <- H
    out$vcov <- I
    out$SE <- sqrt(diag(I))
   
    postprocess.time <- round(proc.time()[3] - pp.start.time, 2)
    
    # Calculate log-likelihood. Done separately as EMtime + postprocess.time is for EM + SEs.
    out$logLik <- joint.log.lik(coeffs, dmats, surv, sv, family, II$b.hat, l0i, l0u, inds, II$Sigma)
    # Collate RE and their variance
    REs <- do.call(rbind, II$b.hat)
    attr(REs, 'Var') <- do.call(rbind, lapply(II$Sigma, diag))
    attr(REs, 'vcov') <- do.call(rbind, lapply(II$Sigma, vech))
    out$REs <- REs
  }else{
    REs <- do.call(rbind, b)
    attr(REs, 'Var') <- do.call(rbind, lapply(update$Sigma, diag))
    attr(REs, 'vcov') <- do.call(rbind, lapply(update$Sigma, vech))
    out$REs <- REs
  }
  comp.time <- round(proc.time()[3] - start.time, 3)
  out$elapsed.time <- c(`EM time` = unname(round(EMend - EMstart, 3)),
                        `Post processing` = if(con$post.process) unname(postprocess.time) else NULL,
                        `Total Computation time` = unname(comp.time),
                        `iterations` = iter)
  if(con$post.process) sv <- surv.mod(surv, formulas, l0, inits.long)
  dmats <- list(long = dmats, surv = sv, ph = surv)
  if(con$return.dmats) out$dmats <- dmats
  
  if(con$return.inits) out$inits = list(inits.long = inits.long,
                                        inits.surv = inits.surv)
  class(out) <- 'joint'
  return(out)
}

#' @method print joint
#' @keywords internal
#' @export
print.joint <- function(x, ...){
  if(!inherits(x, 'joint')) stop('x must be a "joint" object!')
  
  M <- x$ModelInfo
  K <- M$K # Number of responses
  fams <- M$family
  dpsL <- sapply(M$long.formulas, long.formula.to.print, 1)
  dpsD <- sapply(M$disp.formulas, deparse)
  dpsS <- deparse(M$surv.formula)
  
  # Data information
  cat(sprintf("Number of subjects: %d\n", M$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", M$nev, 100 * M$nev/M$n))
  
  # Longitudinal information: 
  cat("\n===================\nModel specification\n===================\n")
  if(K == 1)
    cat("Univariate longitudinal process specification:\n")
  else
    cat("Multivariate longitudinal process specifications: \n")
  
  for(k in 1:K){
    cat(sprintf("%s: %s\n", M$ResponseInfo[k], dpsL[k]))
    if(fams[[k]] %in% c("Gamma", "negbin", "genpois"))
      cat(sprintf("Dispersion model: %s\n",
                  if(dpsD[k]=="~1") "(Intercept) model" else dpsD[k]
      ))
  }
  
  cat("\nSurvival sub-model specification: \n")
  cat(deparse(M$surv.formula))
  
  cat("\n\nAssociation parameter estimates: \n")
  print(setNames(x$coeffs$gamma,
        M$Resps))
  cat("\n")
  invisible(x)
}
