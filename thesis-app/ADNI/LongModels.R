# Longitudinal fits for count data ----------------------------------------
rm(list=ls())
source('.Rprofile')
log.dir <- save.dir.file.path("log")
library(glmmTMB)

# Establish responses
nm.adni <- c("ADAS11", "ADAS13", "MMSE", "FAQ")

# Function to quickly look @ dispersion ------------------------------------
plot.disp <- function(resp){
  sset <- adni[,c("id", resp)]
  sset <- na.omit(sset)
  ids.with.prof <- which(table(sset$id) > 1)
  sset <- sset[sset$id %in% ids.with.prof, ]
  v <- tapply(sset[,resp], sset$id, var)
  m <- tapply(sset[,resp], sset$id, mean)
  plot(m,v,xlab = "E[Y]", ylab = 'Var[Y]', main = resp, pch = 20, cex=.4)
  abline(0,1)
}

par(mfrow=c(2,2))
sapply(nm.adni, plot.disp)
par(mfrow=c(1,1))

# MMSE -> underdispersed
# FAQ: More over- than under- (but both present)
# ADASxx: Roughly equal over-vs-under.

# Function to fit time specs ----------------------------------------------
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

baseline.covars <- c("PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY", "age_scaled")
# combn(baseline.covars, 3) ? Would be way too many!
# combn(baseline.covars, 6)

# Function to fit glmmTMB -------------------------------------------------
log.name <- function(response, time.spec, family) save.dir.file.path(paste0("ADNI_", response, "_", time.spec, "_", family, ".log"), log.dir)

fit.glmmTMB <- function(response, time.spec, family, extras = c("age_scaled", "gender"),
                        disp.formula = ~1){
  family <- switch(family, 
    poisson = poisson(),
    negbin = glmmTMB:::nbinom2(),
    genpois = glmmTMB::genpois()
  )
  
  fm <- get.fit.formula(response, time.spec, extras)
  
  fit <- tryCatch(glmmTMB(formula = fm, data = adni, family = family, dispformula = disp.formula), 
                  warning = function(e) NULL)
  if(is.null(fit)) return(NULL) else return(summary(fit))
}

# Function to output fits to .log file ------------------------------------
fit.to.log <- function(response, time.spec, family){
  f <- fit.glmmTMB(response, time.spec, family)
  if(!is.null(f)){
    sink(log.name(response, time.spec, family))  
    print(f)
    sink()
  }
}

# All possible longit. specs and families (Think negbin will be unwise!)
all.poss <- expand.grid(time = c("int", "linear", "quad", "ns"),
                        family = c("poisson", "genpois"))

for(b in nm.adni){
  cli::cli_alert_info("Starting {b}")
  apply(all.poss, 1, function(x){
    time.spec <- x[1]; family <- x[2]
    fit.to.log(b, time.spec, family)
    cli::cli_alert_success("time spec: {time.spec} for {family} done!")
  })
  cat("\n")
  cat(cli::rule(col = 'green', background_col = 'black'))
  cat("\n")
}

# Summarise which parameterisations fit to each response: -----------------
successful.fits <- gsub("\\.log$", "", gsub("^ADNI\\_", "", dir(log.dir, pattern = "log")))[-1]
sink(save.dir.file.path("_SuccessfulFits.log", log.dir))
cli::rule(center = cli::col_red(" * Successful model fits * "), line_col = "red", width = 81)
print(do.call(rbind, strsplit(successful.fits, '\\_')))
sink()

# ##################
# SUMMARY: 
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