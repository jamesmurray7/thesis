# Create a function to simulate (G)MVJM data ----
# (this same as GMVJM, but with greater amount of objects returned)
create.simfn <- function(family = list("gaussian"),
                    #    family = list("gaussian", "poisson", "binomial"),
                         arguments = list()){
  # Check family
  stopifnot(any(sapply(family, .check.family)))
  # The arguments as default
  base.args <- as.list(args(simData))   # These are calls
  base.args <- lapply(base.args, eval)  # These are the actual simData() objects
  
  # By default simData is bivariate, so need to overwrite this with more apt defaults
  # for the trivariate scenario. Defining replacements below
  K <- length(family)
  if(K > 1) warning("K > 1; originally considering univariate case remember!?")
  # By default we output the trivariate scenario
  beta <- .makebeta(family)
  D <- .makeD(family)
  sigma <- .makesigma(family)
  gamma <- .makegamma(K)
  zeta <- c(0, -0.2)
  # and assigning...
  base.args$n <- 100 # Make a bit smaller by default!
  base.args$return.ranefs <- TRUE # Always return...
  base.args$beta <- beta
  base.args$family <- family
  base.args$sigma <- sigma
  base.args$gamma <- gamma
  base.args$zeta <- zeta
  base.args$D <- D
  # Replace base.args with arguments
  nm <- names(base.args)
  supplied <- names(arguments)
  
  if(any(!supplied%in%nm)){
    stop("At least one supplied ", sQuote("arguments"), " does not match with possible names: ",
         paste0(sapply(nm[-length(nm)], sQuote), collapse=', '), '.')
  }
  
  base.args[(nm <- supplied)] <- arguments # Overwrite with supplied argument(s).
  
  # Create a function which calls simData(...) using the supplied `arguments`
  this.sim.fn <- function(){
    simData(n = base.args$n, ntms = base.args$ntms, fup = base.args$fup, family = base.args$family,
            sigma = base.args$sigma, beta = base.args$beta, D = base.args$D, gamma = base.args$gamma, 
            zeta = base.args$zeta, theta = base.args$theta, cens.rate = base.args$cens.rate,
            regular.times = base.args$regular.times, dof = base.args$dof, random.formulas = base.args$random.formulas, 
            disp.formulas = base.args$disp.formulas, return.ranefs = base.args$return.ranefs)
  }
  
  # And return it
  return(list(fun = this.sim.fn,
         args = base.args, ids = 1:base.args$n))
}

dataGen <- function(X){
  data <- .quietly(X$fun)
  btrue <- data$ranefs
  surv.data <- data$surv.data
  data <- data$data
  out <- list(data = data, 
              surv.data = surv.data,
              btrue = btrue, 
              args = X$args, ids = X$ids)
  class(out) <- "dataGenOut"
  out
}