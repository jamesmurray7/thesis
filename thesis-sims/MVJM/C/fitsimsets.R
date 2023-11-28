N <- 100
# n -----------------------------------------------------------------------
load.dir <- save.dir.file.path('C/data') # location of .RData files to load
out.dir <- save.dir.file.path('C/fits')  # location of resulting .RData joint models
fs::dir_tree(load.dir) # Check this looks right

# Fit ONE trivariate joint model
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ bin
disp.formulas <- list(~1, ~1, ~1) 

# Function to fit ONE model
fit.one <- function(data){
  a <- tryCatch(joint(
    long.formulas = long.formulas, surv.formula = surv.formula, data = data, 
    family = list("gaussian","gaussian","gaussian"), disp.formulas = disp.formulas,
    control = list(return.dmats = FALSE, tol.rel = 5e-3)
  ), error = function(e) NULL)
  return(a)
}

# Function to fit ALL N data sets
fit.all <- function(X){
  stopifnot("list"%in%class(X))
  pb <- utils::txtProgressBar(max = N, style = 3)
  fits <- vector("list", N)
  for(j in 1:N){
    fits[[j]] <- fit.one(X[[j]])
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  return(fits)
}

files <- dir(load.dir)
for(f in files){
  assign("sim.sets", get(load(file.path(load.dir, f))))
  out.file <- save.dir.file.path(f, out.dir)
  cli::cli_alert_info("Starting fits on {f}, this will be saved in {out.file}.")
  fits <- fit.all(sim.sets)
  cli::cli_alert_success("Done, saving in {out.file}")
  save(fits, file = out.file)
}


