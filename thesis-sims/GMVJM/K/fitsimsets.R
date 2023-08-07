N <- 100
# n -----------------------------------------------------------------------
load.dir <- save.dir.file.path('K/data') # location of .RData files to load
out.dir <- save.dir.file.path('K/fits')  # location of resulting .RData joint models
fs::dir_tree(load.dir) # Check this looks right

surv.formula <- Surv(survtime, status) ~ bin

pf <- parent.frame()
fn <- function(K){
  long.formulas <- lapply(1:K, function(k){
    f <- as.formula(paste0('Y.', k, ' ~ time + cont + bin + (1 + time|id)'))
    environment(f) <- pf
    f
  }) 
  disp.formulas <- lapply(1:K, function(k){
    f <- ~1
    environment(f) <- pf
    f
  })
  return(list(long.formulas = long.formulas, disp.formulas = disp.formulas))
}

# Function to fit ONE model
fit.one <- function(data, lf, df, fam){
  a <- tryCatch(joint(
    long.formulas = lf, surv.formula = surv.formula, data = data, 
    family = fam, disp.formulas = df,
    control = list(return.dmats = FALSE)
  ), error = function(e) NULL)
  return(a)
}

# Function to fit ALL N data sets
fit.all <- function(X, K){
  stopifnot("list"%in%class(X))
  pb <- utils::txtProgressBar(max = N, style = 3)
  fits <- vector("list", N)
  # long/disp formulas
  forms <- fn(K)
  longs <- forms$long.formulas; disps <- forms$disp.formulas
  fam <- as.list(rep("gaussian", K))
  
  # The N fits
  for(j in 1:N){
    fits[[j]] <- fit.one(X[[j]], longs, disps, fam)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  return(fits)
}

files <- dir(load.dir)
for(f in files){
  assign("sim.sets", get(load(file.path(load.dir, f))))
  K <- as.numeric(gsub("^K", "", gsub("\\.RData",'',f)))
  out.file <- save.dir.file.path(f, out.dir)
  cli::cli_alert_info("Starting fits on {f}, this will be saved in {out.file}.")
  fits <- fit.all(sim.sets, K)
  cli::cli_alert_success("Done, saving in {out.file}")
  save(fits, file = out.file)
  rm(sim.sets)
}


