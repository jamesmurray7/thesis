# Dispersion models for count data ----------------------------------------
rm(list=ls())
source('.Rprofile')
log.dir <- save.dir.file.path("log")
library(glmmTMB)

# ##################
# In LongModels.R we determined:
# MMSE: 
#   Poisson:
#     intercept only.
# ADAS11:
#   Poisson:
#     Intercept only;
#     Linear. 
#   Genpois:
#     Intercept only;
#     Cubic splines.
# ADAS13:
#   Poisson: 
#     Intercept only;
#     Linear. 
#   Genpois:
#     Intercept only;
#     Cubic splines.
# FAQ: 
#   Poisson:
#     ALL.
#   Genpois:
#     ALL.
# ##################

# "Base" model fit statistics ---------------------------------------------
# For means of comparison, what were model fit statistics of these "null" dispmodels?
# (i.e. dispersion model of ~1).
gp.fits <- dir(log.dir, pattern = 'genpois')
test <- readLines(file.path(log.dir,gp.fits[1]))
aictab <- el(strsplit(trimws(test[which(grepl("BIC",test))+1]),"\\s+"))
ExtractFrom.log <- function(fn){ # one of gp.fits
  # Read it in line-by-line
  test <- readLines(file.path(log.dir, fn)) 
  # extract line containing AIC/BIC/logLik/dev/df.resid
  aictab <- el(strsplit(trimws(test[which(grepl("BIC",test))+1]), "\\s+"))
  aictab <- setNames(as.numeric(aictab), c("AIC", "BIC", "logLik", "deviance", "dfResid"))
  
  return(aictab)
}

ints <- sapply(gp.fits, ExtractFrom.log)


# Copying across functions to fit glmmTMB (lazy) --------------------------
## Function to fit time specs //
pf <- parent.frame()
get.fit.formula <- function(response = 'ADAS13', time.spec = 'linear', extras = c("age_scaled", "gender")){
  extras <- paste0(extras, collapse = ' + ')
  # RHS
  time <- switch(time.spec,
                 linear = " APOE4 * time + (1 + time|id)",
                 int = " APOE4 * time + (1|id)",
                 quad = " APOE4 * (time + I(time^2)) + (1 + time + I(time^2)|id)",
                 ns = " APOE4 * splines::ns(time, knots = c(0.5, 1.5)) + (1 + splines::ns(time, knots = c(0.5, 1.5))|id)")
  
  form <- if(nchar(extras > 1)) paste0(response, " ~ ", extras, ' +', time) else paste0(response, " ~ ", time)
  form <- as.formula(form)
  environment(form) <- pf
  return(form)
}

## Function to fit the GLMM //
fit.glmmTMB <- function(response, time.spec, family, extras = c("age_scaled", "gender"),
                        disp.formula = NULL, summ = T){
  family <- switch(family, 
                   poisson = poisson(),
                   negbin = glmmTMB:::nbinom2(),
                   genpois = glmmTMB::genpois()
  )
  
  fm <- get.fit.formula(response, time.spec, extras)
  
  if(is.null(disp.formula)) stop("dispformula unspecified!\n\n")
  
  disp.formula <- as.formula(disp.formula)
  environment(disp.formula) <- pf
  
  fit <- tryCatch(glmmTMB(formula = fm, data = adni, family = family, dispformula = disp.formula), 
                  warning = function(e) NULL)#
  if(summ)
    if(is.null(fit)) return(NULL) else return(summary(fit))
  else
    if(is.null(fit)) return(NULL) else return(fit)
}

# Three candidate disp.formulas -------------------------------------------
disp.formulas <- list("~time", "~APOE4", "~APOE4 * time") # Presence of gene -> more/less variable?
nd <- length(disp.formulas)
# Work out responses/time.spec we need to re-fit with different disp.formulae
to.refit <- do.call(rbind, strsplit(gsub("\\_gen.*$","", gsub("^ADNI\\_", "", gp.fits)),"\\_"))

refits <- setNames(apply(to.refit, 1, function(x){
  resp <- x[1]; time.spec <- x[2]
  family <- "genpois"
  cat(cli::rule(center = cli::col_blue(paste0(resp, ", ", time.spec)), line_col = "orange"))
  cat("\n")
  out <- setNames(vector("list", nd),  gsub("\\~","", disp.formulas))
  for(i in 1:nd){
    this <- suppressMessages(fit.glmmTMB(resp, time.spec, family, disp.formula = disp.formulas[[i]]))
    if(is.null(this)){
      cli::cli_alert_warning(sprintf("Failed for %s.", disp.formulas[[i]]))
    }else{
      cli::cli_alert_success(sprintf("Successful fit for %s.", disp.formulas[[i]]))
      out[[i]] <- this
    }
  }
  out
}), apply(to.refit, 1, paste0, collapse = '_'))


# Big print-out -----------------------------------------------------------
nms <- names(refits)
invisible(lapply(seq_along(refits), function(i){
  r <- refits[[i]]; nm <- nms[i]
  ss <- el(strsplit(nm, '_'))
  resp <- ss[1]; ts <- ss[2]
  # Get intercept model AIC tab to print --
  .lookup <- paste0("ADNI_", nm, "_genpois.log")
  .lookup <- ints[,colnames(ints)==.lookup]
  cat(sprintf("Response: %s; time spec: %s\n", resp, ts))
  cat("AICtab on ~1 dispersion model:\n")
  print(.lookup); cat("\n")
  # Check _any_ completed fits:
  if(sum(sapply(r, is.null)) == nd){
    cat("No other dispersion models successfully fit for this time spec!\n\n")
  }else{
    tabs <- lapply(r, '[[', "AICtab")
    succ.fit <- paste0(names(tabs)[!sapply(tabs, is.null)], collapse = ', ')
    cat(sprintf("Successful fits for dispersion models: %s\n", succ.fit))
    print(do.call(rbind, tabs))
    cat("\n\n")
  }
  cat(cli::rule(line_col = "orange", line = 2))
  cat("\n")
  NULL
}))

# #####################
# Summary for GENPOIS fits (checking Poisson after!)
# ADAS11: 
#   ts: ns --
#    -> ~1 best.
# ADAS13: 
#   ts: int --
#     -> ~time best (12134.56/12184.93)
#   ts: linear --
#     -> ~time best (11928.72/11990.29)
#   ts: ns --
#     -> ~1 best (11840.2/11957.7)
# FAQ:
#   ts: int --
#     -> ~time best (9307.677/9358.143)
#   ts: linear --
#     -> ~1 best (9184.6/9240.7)
#   ts: ns --
#     -> ~1 best (9168.9/9286.6)
#   ts: quad --
#     -> ~1 best (9166.8/9250.9)
# #####################

aa <- (model.response(ADAS13.best.GP2$frame)-fitted(ADAS13.best.GP2))/sqrt(V(fitted(ADAS13.best.GP2), sigma(ADAS13.best.GP2)))
bb <- (model.response(ADAS13.best.GP2$frame)-fitted(ADAS13.best.GP2))/sqrt(V(fitted(ADAS13.best.GP2), predict(ADAS13.best.GP2, type = 'disp')))

# Function (i think) which obtains Pearson residuals
get.Pearson <- function(fit){
  V <- family(fit)$variance
  # mu <- predict(fit, type = 'conditional') # E[Y] (same as fitted(fit)!)
  return(
    (model.response(fit$frame) - fitted(fit))/sqrt(V(fitted(fit), predict(fit, type = 'disp')))
  )
}

plot.Resids <- function(x, main = ""){
  if(main == "Poisson"){
    plot(resid(x, "pearson") ~ fitted(x), main = main, pch = 20, cex = 0.25)
  }else{
    plot(get.Pearson(x) ~ fitted(x), main = main, pch = 20, cex = 0.25)
  }
}


# Comparing with Poisson fits ---------------------------------------------
# ADAS11 //
# Poisson linear vs genpois: ns (~1)
par(mfrow=c(1,2))
ADAS11.best.Po <- fit.glmmTMB("ADAS11", "linear", "poisson", disp.formula = "~1", summ = F)
plot.Resids(ADAS11.best.Po, "Poisson")
ADAS11.best.GP <- fit.glmmTMB("ADAS11", "ns", "genpois", disp.formula = "~1", summ = F)
plot.Resids(ADAS11.best.GP, "GP")
par(mfrow=c(1,1)) # Look extremely similar; Poisson slightly better?

# ADAS13 //
ADAS13.best.Po <- fit.glmmTMB("ADAS13", "linear", "poisson", disp.formula = "~1", summ = F)
ADAS13.best.GP1 <- fit.glmmTMB("ADAS13", "linear", "genpois", disp.formula = "~time", summ = F)
ADAS13.best.GP2 <- fit.glmmTMB("ADAS13", "ns", "genpois", disp.formula = "~1", summ = F)
par(mfrow=c(1,3))
plot.Resids(ADAS13.best.Po, "Poisson")
plot.Resids(ADAS13.best.GP1, "GP - linear")
plot.Resids(ADAS13.best.GP2, "GP - ns")
par(mfrow=c(1,1)) # Similar again

# FAQ
FAQ.best.Po <- fit.glmmTMB("FAQ", "ns", "poisson", disp.formula = "~1", summ = F)
FAQ.best.GP <- fit.glmmTMB("FAQ", "ns", "genpois", disp.formula = "~1", summ = F)
par(mfrow=c(1,2))
plot.Resids(FAQ.best.Po, "Poisson")
plot.Resids(FAQ.best.GP, "GP")
par(mfrow=c(1,1)) # Similar again

