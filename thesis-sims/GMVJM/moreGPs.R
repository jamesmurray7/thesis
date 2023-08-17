# Need 75 more GP fits
fits <- vector("list", 75L)
fn <- create.simfn(family = list("genpois"),
                   arguments = list(n=500, theta = c(-2.8, .1), 
                                    beta = c(1,.05,-.05,.1),
                                    sigma = list(-0.3),
                                    random.formulas = list(~1), D = matrix(.30,1,1)))
successful.fits <- 0; p <- 1
repeat{
  this.data <- try(fn(), silent = T)
  if(inherits(this.data, "try-error")) next # Data generation can fail sometimes, GP is finicky!
  
  this.fit <- try(joint(
    long.formulas = list(Y.1 ~ time + cont + bin + (1|id)),
    surv.formula = Surv(survtime, status) ~ bin,
    data = this.data, family = list("genpois"),
    disp.formulas = list(~1),
    control = list(tol.rel = 5e-3, return.dmats = FALSE)
  ), silent = T)
  
  if(inherits(this.fit, "try-error")) next # Skip to next attempt if fit fails.
  
  # Successful if got this far!
  fits[[p]] <- this.fit
  successful.fits <- successful.fits + 1
  p <- p + 1
  
  # Print out short message with current status
  cat(sprintf("\n===\nCurrent number of successful fits: %d\n===\n\n", successful.fits))
  if(successful.fits == 75) break # Stop at 75 successful fits (# needed to make 300!)
}
save(fits, file = "/data/c0061461/THESIS/GMVJM/fits/TEMP_UnivGP.RData")

# surgery on existing GP fits
assign("newfits", get(load("/data/c0061461/THESIS/GMVJM/fits/TEMP_UnivGP.RData")))
load("/data/c0061461/THESIS/GMVJM/fits/UnivGP.RData")

nn <- NROW(fits)
p <- 1
for(i in 1:nn){
  if(is.null(fits[[i]])){
    fits[[i]] <- newfits[[p]]
    cat(sprintf("Replaced element %d of `fits` with element %d of `newfits`\n", i, p))
    p <- p + 1
  }
}

save(fits, file = "/data/c0061461/THESIS/GMVJM/fits/UnivGP_Surgery.RData")