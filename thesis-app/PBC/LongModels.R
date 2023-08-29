# Longitudinal fits for count data ----------------------------------------
rm(list=ls())
source('.Rprofile')
library(glmmTMB)
log.dir <- save.dir.file.path("Longs")
pf <- parent.frame()

# Make a sub-directory for each response ----------------------------------
# (or else significant flooding will occur!)

(resps <- names(PBCredx)[8:15])
for(r in resps){
  .this <- save.dir.file.path(r, log.dir)
  if(!fs::dir_exists(.this)){
    cat(sprintf("\n==\nCreating log directory %s\n==\n\n", .this))
    fs::dir_create(.this)
  }else{
    cat(sprintf("The directory %s already exists, skipping...", .this))
  }
  cat("\n")
}

# Function which makes all we need ----------------------------------------
# Outputs: 
#  1) Formula for glmmTMB.
#  2) file name for output to be stored in.
#  3) The directory for 2).
#  4) Number of coefficients (i.e. beta terms) minus intercept.
#  5) R family object requested.
#  6) family name as supplied to function.

toSentence <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

file.name <- function(drugint, time.spec, extras, family, cc){
  ll <- if(drugint) "_DI" else NULL
  ex <- paste0(toSentence(extras), collapse = "_")
  paste0(cc, toSentence(family), "_", time.spec, ll, "_", ex, ".log")
}

make.fit.objects <- function(response = 'serBilir', 
                             time.spec = 'linear', 
                             extras = c("sex", "age", "histologic"),
                             drugint = TRUE,  # This will be a LINEAR combination only.
                             family = "gaussian"){ 
  # Checking for typos/improper specs.
  stopifnot(response%in%resps) 
  extras <- unname(sapply(extras, match.arg, extras))
  # Creating long.formula -------------------------------------------------
  # "Extra" baseline covariates
  rhs.extras <- paste0(extras, collapse = " + ")
  cc <- length(extras)
  # Time specification (linear/quadratic/cubic splines).
  # Also counting up number of coefficients (sans intercepts).
  time <- switch(time.spec,
                 linear = {
                   if(drugint)
                     cc <- cc + 3
                   else
                     cc <- cc + 1
                   "time + (1 + time|id)"
                  },
                 int = {
                   if(drugint)
                     cc <- cc + 3
                   else
                     cc <- cc + 1
                   "time + (1|id)"
                 },
                 quad = {
                   if(drugint)
                     cc <- cc + 5
                   else
                     cc <- cc + 2
                   "(time + I(time^2)) + (1 + time + I(time^2)|id)"
                 },
                 ns = {
                   if(drugint)
                     cc <- cc + 7
                   else
                     cc <- cc + 3
                   "splines::ns(time, knots = c(1, 4)) + (1 + splines::ns(time, knots = c(1, 4))|id)"
                 })
  # Model with drug interaction on above time specification?
  if(drugint) time <- paste0("drug * ", time)
  
  # Put everything together
  form <- as.formula(paste0(response, " ~ ", rhs.extras, " + ", time), pf)
  
  # Creating file.name ----------------------------------------------------
  out.dir <- save.dir.file.path(response, log.dir)
  fn <- file.name(drugint, time.spec, extras, family, cc)
  
  # Parsing family input for direct use later -----------------------------
  family.txt <- family
  family <- switch(family, 
                   poisson = poisson(),
                   gaussian = gaussian(),
                   binomial = binomial(),
                   Gamma = Gamma(link = "log"),
                   negbin = glmmTMB:::nbinom2(),
                   genpois = glmmTMB::genpois()
            )
  # Return everything
  list(
    longformula = form,
    output.directory = out.dir,
    filename = fn,
    Pcount = cc,
    family = family,
    familytext = family.txt
  )
  
}

# Function to fit glmmTMB -------------------------------------------------
# Input: longformula and family from `make.fit.objects` above.
# Output: glmmTMB object or NULL (if fit fails).
fit.glmmTMB <- function(longformula, family){
  fit <- tryCatch(suppressMessages(
                    suppressWarnings(
                      glmmTMB(formula = longformula, 
                              data = PBCredx, 
                              family = family)
                    )
                  ), 
                  warning = function(e) NULL,
                  error = function(e) NULL)
  if(is.null(fit)) return(NULL) else return(fit)
}

# Sorting out all models --------------------------------------------------
resps
extras <- c("sex", "age", "histologic")
var2 <- combn(extras, 2) # Bivariate combos;
var1 <- combn(extras, 1) # univariates
time.specs <- c("int", "linear", "quad", "ns")

for(r in resps){ # Probably a nicer way to do this but it doesn't matter.
  # Entering families
  xx <- readline(sprintf("Enter family for %s --> ", r))
  if(xx == "next" | xx == "skip") next
  if(xx == "stop") break
  
  # Creating an object to store results in (so can consider different families).
  eval(parse(text = paste0(r, xx, "<<-vector('list')")))
  for(j in 0:1){ # Outer loop: Drug interaction.
    druginteract <- as.logical(j); key <- if(druginteract) "DI" else ""
    for(s in time.specs){ # Inner loop: Time specification
      # Trivariate spec:
      mf3 <- make.fit.objects(r, s, extras, druginteract, xx)
      # Store
      eval(parse(text = paste0(r, xx,'[["',s,key,'triv"]]<-mf3'))) # Largely a junk lookup since it won't get used.
      # Bivariate extra items:
      mf2 <- apply(var2, 2, function(x) make.fit.objects(r,s,x,druginteract,xx))
      eval(parse(text = paste0(r, xx,'[["',s,key,'biv"]]<-mf2'))) 
      # Univariate:
      mf1 <- apply(var2, 2, function(x) make.fit.objects(r,s,x,druginteract,xx))
      eval(parse(text = paste0(r, xx,'[["',s,key,'univ"]]<-mf1'))) 
    }
  }
  cat("\n")
  cat(cli::rule(center = paste("List created for", toSentence(xx), r, "models!"),
                col = "white", background_col = "black"))
  cat("\n\n")
}


# Fitting the models ------------------------------------------------------

# This makes A LOT of log files!!!!
make.all.logs <- function(X){
  stopifnot("list"%in%class(X))  
  lapply(X, function(x){
    if(NROW(x)!=6){
      lapply(x, function(y){
        ff <- fit.glmmTMB(y$longformula, y$family)
        if(!is.null(ff)){
          sink(save.dir.file.path(y$filename, y$output.directory))
          print(summary(ff))
          sink()
        }
      })
    }else{
      ff <- fit.glmmTMB(x$longformula, x$family)
      if(!is.null(ff)){
        sink(save.dir.file.path(x$filename, x$output.directory))
        print(summary(ff))
        sink()
      }
    }
    capture.output(rnorm(1), file = nullfile())
  })
  objname <- deparse(substitute(X))
  cat(cli::rule(center = paste0("`", objname, "` done!"),
                col = "white", background_col = "grey2", line = 2))
}

# Albumin
make.all.logs(albumingaussian)
# Alkaline
make.all.logs(alkalinepoisson)
# Hepatomegaly
make.all.logs(hepatomegalybinomial)
# spiders
make.all.logs(spidersbinomial)
# SGOT/AST
make.all.logs(SGOTgaussian)
# platelets
make.all.logs(plateletspoisson)
# proth.
make.all.logs(prothrombinGamma)
# ser. bili
make.all.logs(serBilirgaussian)
# ser. alb.
make.all.logs(albumingaussian)


