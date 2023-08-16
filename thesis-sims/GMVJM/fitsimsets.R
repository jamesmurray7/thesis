rm(list=ls())
source(".Rprofile")
N <- 300L
out.dir <- save.dir.file.path("fits")
data.dir <- save.dir.file.path("data")
files <- dir(data.dir, pattern = "\\.RData")


# Setting out (sub-)model formulae ----------------------------------------
surv.formula <- Surv(survtime, status) ~ bin

pf <- parent.frame()
make.long.formula <- function(key){
  # Trivariate case ---
  if(key == "Triv"){
    long.formulas <- list(
      Y.1 ~ time + cont + bin + (1 + time|id),
      Y.2 ~ time + cont + bin + (1 + time|id),
      Y.3 ~ time + cont + bin + (1|id)
    )
    disp.formulas <- lapply(1:3, function(f){
      f <- ~1
      environment(f) <- pf
      f
    })
  }else{ # All others are univariate...
    if(key != "GP"){
      long.formulas <- list(
        Y.1 ~ time + cont + bin + (1 + time|id)
      )
    }else{
      long.formulas <- list(
        Y.1 ~ time + cont + bin + (1|id)
      )
    }
    disp.formulas <- lapply(1:1, function(f){
      f <- ~1
      environment(f) <- pf
      f
    })
  }
  
  long.formulas <- lapply(long.formulas, function(f){
    environment(f) <- pf
    f
  })
  
  return(list(long.formulas = long.formulas, disp.formulas = disp.formulas))
}

# Fit one/ALL datasets ----------------------------------------------------
fit.one <- function(data, lf, df, fam){
  a <- tryCatch(joint(
    long.formulas = lf, surv.formula = surv.formula, data = data,
    family = fam, disp.formulas = df,
    control = list(return.dmats = FALSE, tol.rel = 5e-3)
  ), error = function(e) NULL)
  return(a)
}

fit.all <- function(X, key){
  stopifnot("list"%in%class(X))
  
  # Get family from key
  fam <- switch(key,
                Triv = list("gaussian", "poisson", "binomial"),
                Ga = list("Gamma"),
                GP = list("genpois"),
                NB = list("negbin"))
  forms <- make.long.formula(key)
  longs <- forms$long.formulas; disps <- forms$disp.formulas
  
  # The N fits
  pb <- utils::txtProgressBar(max = N, style = 3)
  fits <- vector("list", N)
  for(j in 1:N){
    fits[[j]] <- fit.one(X[[j]], longs, disps, fam)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  return(fits)
}

# >>hotfix 16/08/23 refitting safe GP1.
for(f in files[3]){#files[-c(1:2)]){
  assign("sim.sets", get(load(save.dir.file.path(f, data.dir))))
  out.file <- save.dir.file.path(f, out.dir)
  key <- gsub("\\.RData|^Univ", "", f)
  cli::cli_alert_info("Starting fits on {f}. Key: {key}.")
  fits <- fit.all(sim.sets, key)
  cli::cli_alert_success("Done, saving in {out.file}")
  save(fits, file = out.file)
}







