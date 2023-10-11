rm(list=ls())
source(".Rprofile")
library(splines)
disp.dir <- save.dir.file.path("Disps")
# file.show("printout.txt", title = "models for reference") # Changed now!
resps <- c("prothrombin", "platelets", "alkaline")
pf <- parent.frame()
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)

# Intercept-only is already in the Longs dir!!
univs <- poss <- c("time", "drug", "age", "histologic2", "sex")
# biv  <- combn(poss, 2)
# triv <- combn(poss, 3)
# quad <- combn(poss, 4)
# five <- c(combn(poss, 5))

# Get long.formula --------------------------------------------------------
pf <- parent.frame()
get.long.stuff <- function(resp){
  form <- switch(resp,
                   prothrombin = as.formula(
                     "prothrombin ~ histologic2 + sex + (time + I(time^2)) + (1 + time + I(time^2)|id)",
                     pf
                   ),
                   platelets = as.formula(
                     "platelets ~ age + histologic2 + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)",
                     pf
                   ),
                   alkaline = as.formula(
                     "alkaline ~ age + histologic2 + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)",
                     pf
                   )
                 )
  fam <- switch(resp,
                prothrombin = Gamma(link="log"),
                platelets = glmmTMB::genpois(),
                alkaline = glmmTMB::nbinom2()
  )
  list(long.formula = form, family = fam)
}


# File-naming functions ---------------------------------------------------
toSentence <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

file.name <- function(resp, vars, cc){
  ex <- paste0(toSentence(vars), collapse = "_")
  save.dir.file.path(paste0(cc, "_", ex, ".log"), save.dir.file.path(resp, disp.dir))
}


# Creating dispersion model formula ---------------------------------------
make.disp.formula <- function(dispvars = c("time", "drug", "age", "histologic2", "sex")){
  dispvars <- unname(sapply(dispvars, match.arg, dispvars))
  # Construct this candidate dispersion model
  dispvars2 <- paste0(dispvars, collapse = " + ")
  as.formula(paste0("~", dispvars2), pf)
}

# Fitting via glmmTMB -----------------------------------------------------
fit.glmmTMB <- function(long.formula, family, disp.formula){
  fit <- tryCatch(suppressMessages(
                    suppressWarnings(
                      glmmTMB(formula = long.formula, data = PBCredx,
                              family = family, dispformula = disp.formula)
                    )
                  ),
                  warning = function(e) NULL,
                  error = function(e) NULL)
  if(is.null(fit)) return(NULL) else return(fit)
}

# All-in-one --------------------------------------------------------------
one.fn <- function(resp, long.formula, family, dispvars){
  cc <- length(dispvars)
  disp.formula <- make.disp.formula(dispvars)
  out <- file.name(resp, dispvars, cc)
  cat(sprintf("\nDoing %s...\n", deparse(disp.formula)))
  ff <- fit.glmmTMB(long.formula, family, disp.formula)
  if(!is.null(ff)){
    sink(file = out)
    print(summary(ff))
    sink()
    cat(sprintf("Successful, saved in %s\n", out))
  }else{
    cat(sprintf("Unsuccessful fit!\n"))
  }
}

for(r in resps){
  cat(sprintf("\n===\nStarting fits for %s\n===\n\n", r))
  aa <- get.long.stuff(r)
  long <- aa$long.formula; family <- aa$family
  cat(cli::rule(center = "Univariates", col = .nice.orange, background_col = "magenta4", line = 2))
  
  sapply(univs, function(x) one.fn(r, long, family, x))
    cat("\n")

    cat(cli::rule(center = paste0("`", r, "` done!"),
                  col = "white", background_col = "grey2", line = 2))

    xx <- readline("Press enter for next >>>>> ")
}

# OLD --> Only considering univariates.
# Loop over resps, fit and store summary if not NULL fit.
# for(r in resps){
#   cat(sprintf("\n===\nStarting fits for %s\n===\n\n", r))
#   aa <- get.long.stuff(r)
#   long <- aa$long.formula; family <- aa$family
#   
#   # Do five-variate first
#   cat(cli::rule(center = "Five-variate", col = .nice.orange, background_col = "magenta4", line = 2))
#   one.fn(r, long, family, five)
#   cat("\n")
#   
#   # Do all four-variates
#   cat(cli::rule(center = "Four-variates", col = .nice.orange, background_col = "magenta4", line = 2))
#   apply(quad, 2, function(x) one.fn(r, long, family, x))
#   cat("\n")
#   
#   # Do all three-variates
#   cat(cli::rule(center = "Trivariates", col = .nice.orange, background_col = "magenta4", line = 2))
#   apply(triv, 2, function(x) one.fn(r, long, family, x))
#   cat("\n")
#   
#   # Do all bivariates
#   cat(cli::rule(center = "Bivariates", col = .nice.orange, background_col = "magenta4", line = 2))
#   apply(biv, 2, function(x) one.fn(r, long, family, x))
#   cat("\n")
#   
#   # Finally, all univariates
#   cat(cli::rule(center = "Univariates", col = .nice.orange, background_col = "magenta4", line = 2))
#   sapply(univs, function(x) one.fn(r, long, family, x))
#   cat("\n")
#     
#   cat(cli::rule(center = paste0("`", r, "` done!"),
#                 col = "white", background_col = "grey2", line = 2))
#   
#   xx <- readline("Press enter for next >>>>> ")
# }
