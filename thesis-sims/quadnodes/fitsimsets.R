rm(list=ls())
source(".Rprofile")
N <- 200L
out.dir <- save.dir.file.path("fits")
data.dir <- save.dir.file.path("data")
files <- dir(data.dir, pattern = "\\.RData")

# Setting out candidate quadrature nodes ----------------------------------
rhos <- c(2, 3, 5, 7, 9, 15, 
          97, 98, 99) # soft-key 97 for 10%; 98 for 25%; and 99 for 33% 

# Setting out (sub-)model formulae ----------------------------------------
pf <- parent.frame()
surv.formula <- as.formula(Surv(survtime, status) ~ bin,pf)


make.long.formula <- function(){
 # Creating formula objects.
  long.formulas <- list(
    Y.1 ~ time + cont + bin + (1 + time|id),
    Y.2 ~ time + cont + bin + (1 + time|id),
    Y.3 ~ time + cont + bin + (1|id)
  )
  lapply(long.formulas, `environment<-`, pf) # Neat little workaround!
  
  disp.formulas <- lapply(1:3, function(f){
    f <- ~1
    environment(f) <- pf
    f
  })
  
  return(list(long.formulas = long.formulas, 
              disp.formulas = disp.formulas))
}

# Fit one/ALL datasets ----------------------------------------------------
fit.one <- function(data, lf, df, gh.nodes){
  a <- tryCatch(joint(
    long.formulas = lf, surv.formula = surv.formula, data = data,
    family = list("gaussian", "poisson", "binomial"), disp.formulas = df,
    control = list(return.dmats = FALSE, tol.rel = 5e-3, gh.nodes = gh.nodes)
  ), error = function(e) NULL)
  return(a)
}

fit.all <- function(X, quad.vec){
  stopifnot("list"%in%class(X))
  
  forms <- make.long.formula()
  longs <- forms$long.formulas; disps <- forms$disp.formulas
  
  # The N fits
  pb <- utils::txtProgressBar(max = N, style = 3)
  fits <- vector("list", N)
  for(j in 1:N){
    fits[[j]] <- fit.one(X[[j]], longs, disps, quad.vec[j])
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  return(fits)
}


# Load data, work out nodes and fit ---------------------------------------
load(save.dir.file.path(files, data.dir)) # `sim.sets`

# Work out quad.node.grid
quad.node.grid <- sapply(rhos, function(r){
  if(r < 90){
    return(rep(r, N))
  }else{
    if(r == 97) # probably a more elegant way to do this!
      return(sapply(sim.sets, stringer.gh, 10))
    else if(r == 98)
      return(sapply(sim.sets, stringer.gh, 25))
    else if(r == 99)
      return(sapply(sim.sets, stringer.gh, 33))
  }
})
colnames(quad.node.grid) <- rhos
nrhos <- ncol(rhos)

file.name <- function(r) save.dir.file.path(paste0("gh", r, ".RData"), save.dir.file.path("fits"))

for(r in rhos){
  gh.r <- quad.node.grid[,as.character(r)]
  out.file <- file.name(r)
  
  cli::cli_alert_info("Starting fits on rho = {r}. Quadrature node vector: ")
  cat(gh.r, sep = ", ")
  cat("\n\n")
  fits <- fit.all(sim.sets, gh.r)
  cli::cli_alert_success("Done, saving in {out.file}")
  save(fits, file = out.file)
  rm(gh.r)
}







