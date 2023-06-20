# Create a function to simulate MVJM data.
create.simfn <- function(arguments = list(), K = 3){
  # The arguments as default
  base.args <- as.list(args(simData))   # These are calls
  base.args <- lapply(base.args, eval)  # These are the actual simData() objects
  
  # By default simData is bivariate, so need to overwrite this with more apt defaults
  # for the trivariate scenario. Defining replacements below
  beta <- .makebeta(K)
  family <- as.list(rep("gaussian", K))
  sigma <- as.list(rep(0.16, K))
  gamma <- .makegamma(K)
  zeta <- c(0, -0.2)
  D <- .makeD(K)
  # and assigning...
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
            disp.formulas = base.args$disp.formulas, return.ranefs = base.args$return.ranefs)$data
  }
  
  # And return it
  return(this.sim.fn)
}

# 'fun' from above, N the number of simulations.
createNsims <- function(fun, N){
  cli::cli_progress_bar(name = "Generating simulated joint data...", total = N)
  out <- vector("list", N)
  for(j in 1:N){
    cli::cli_progress_update(inc = 1)
    out[[j]] <- .quietly(fun)
  }
  cli::cli_progress_done()
  return(out)
}

# Report failure rate(s) for a generated set of data
DataSummary <- function(data){
  N <- NROW(data)
  #cat("Number of sets:", N, "\n")
  out <- sapply(data, function(x){
    unq <- x[!duplicated(x$id), c("id", "status", "survtime")]
    n <- nrow(unq)
    d <- sum(unq$status)/n * 100
    qn <- median(unq[unq$status == 1, "survtime"])
    c(n, d, qn)
  })
  cat("Data summary ---------\n")
  cat(sprintf("n: %d\n", out[1,1]))
  cat(sprintf("Average failure rate: %.2f%%\n", mean(out[2,])))
  cat(sprintf("Median failure time: %.2f\n", median(out[3,])))
  cat("\n")
}
