pf <- parent.frame()
surv.formula <- Surv(survtime, status) ~ bin

# Creating necessary formula objects --------------------------------------
create.K.long.formulas <- function(method = c("gmvjoint", "joineRML", "JMbayes2"), K = 3){
  method <- match.arg(method)
  if(method == "JMbayes2"){
    # *: for some reason JMbayes2 runs my .Rprofile on initialisation, so temporarily move it?
    #fs::file_move("./.Rprofile", "../.Rprofile")
    fixed <- lapply(1:K, function(k){
      f <- paste0("Y.", k, ' ~ time + cont + bin')
      f <- as.formula(f)
      environment(f) <- pf
      return(f)
    })
    random <- lapply(1:K, function(k){
      f <- ~ time|id
      environment(f) <- pf
      return(f)
    })
    return(list(method = method, fixed = fixed, random = random))
    #fs::file_move("../.Rprofile", "./.Rprofile") # Move it back to this wd
  }else if(method == "joineRML"){
   formLongFixed <- lapply(1:K, function(k){
     f <- paste0("Y.", k, ' ~ time + cont + bin')
     f <- as.formula(f)
     environment(f) <- pf
     return(f)
   })
   formLongRandom <- lapply(1:K, function(k){
     f <- ~ 1 + time|id
     environment(f) <- pf
     return(f)
   })
   return(list(method = method, 
               formLongRandom = formLongRandom,
               formLongFixed = formLongFixed))
  }else{
    long.formulas <- lapply(1:K, function(k){
      f <- paste0("Y.", k, ' ~ time + cont + bin + (1 + time|id)')
      f <- as.formula(f)
      environment(f) <- pf
      return(f)
    })
    disp.formulas <- lapply(1:K, function(k){
      f <- ~1
      environment(f) <- pf
      return(f)
    })
    return(list(method = method, 
                long.formulas = long.formulas,
                disp.formulas = disp.formulas))
  }
}

# Function to fit ONE Kvariate joint model --------------------------------
fit.one <- function(data, method = c("gmvjoint", "joineRML", "JMbayes2"),
                    formula.list = list(), ...){
  method <- match.arg(method)
  if(method == "JMbayes2"){
    # *: for some reason JMbayes2 runs my .Rprofile on initialisation, so temporarily move it?
    fs::file_move("./.Rprofile", "../.Rprofile")
    fixed <- formula.list$fixed; random <- formula.list$random
    # ...
    M <- lapply(seq_along(fixed), function(k){
      f1 <- nlme::lme(fixed = fixed[[k]], random = random[[k]], data = data)
    })
    S <- coxph(Surv(survtime, status) ~ bin, data = data[!duplicated(data$id), ])
    f <- tryCatch(
      JMbayes2::jm(S, M, "time", control = list(
        n_chains = 1L, cores = 1L,
        ...
      ), parallel = "multicore")
      ,error= function(e) NULL
    )
    fs::file_move("../.Rprofile", "./.Rprofile")
    return(f)
  }else if(method == "gmvjoint"){
    longs <- formula.list$long.formulas; disps <- formula.list$disp.formulas
    family <- as.list(rep("gaussian", length(longs)))
    a <- tryCatch(joint(
      long.formulas = longs, surv.formula = surv.formula, data = data,
      family = family, disp.formulas = disps,
      control = list(
        return.dmats = FALSE,
        tol.rel = 5e-3, ...
      )
    ),
    error = function(e) NULL)
  }else{
    fLF <- formula.list$formLongFixed; fLR <- formula.list$formLongRandom
    f <- tryCatch(suppressMessages(mjoint(
      formLongFixed = fLF, formLongRandom = fLR, formSurv = surv.formula,
      data = data, timeVar = "time", 
      control = list(
        type = "sobol", tol.em = 1e-2, tol2 = 5e-3, convCrit = "sas"
      )
    )), error = function(e) NULL)
    
    if(!is.null(f)){
      a <- list(
        coeffs = f$coefficients,
        elapsed = f$comp.time,
        SE = sqrt(diag(solve(f$Hessian))),
        REs = f$Eb,
        REs.var = f$Vb
      )
      return(a)
    }else{
      return(NULL)
    }
  }
}

# Fit ALL N data sets -----------------------------------------------------
fit.all <- function(X, # X list of data
                    method = c("gmvjoint", "joineRML", "JMbayes2"), K = 3,
                    ...){
  stopifnot("list"%in%class(X))
  method <- match.arg(method)
  pb <- utils::txtProgressBar(max = N, style = 3)
  fits <- vector("list", N)
  formula.list <- create.K.long.formulas(method = method, K = K)
  for(j in 1:N){
    fits[[j]] <- fit.one(X[[j]], method, formula.list, ...)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  return(fits)
}


