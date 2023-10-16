# UPDATED 12/10/23
rm(list=ls())
source('.Rprofile')
library(splines)
log.dir <- save.dir.file.path("Joint/Univs")
head(PBCredx)
pf <- parent.frame()
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)

(resps <- names(PBCredx)[8:15])
resps <- resps[resps!="spiders"]

# Function to fit all based on best model
surv.formula <- Surv(survtime, status) ~ age + sex + histologic2

# Best-fitting longitudinal formula ---------------------------------------
get.formulae <- function(resp){
  stopifnot(resp%in%resps)
  out <- switch(resp,
                albumin = list(
                  long.formula = list(albumin ~ time + age + histologic2 + (1 + time|id)),
                  family = list("gaussian")
                ),
                alkaline = list(
                  long.formula = list(alkaline ~ ns(time, knots = c(1,4)) + age + histologic2 + (1 + ns(time, knots = c(1,4))|id)),
                  family = list("negbin"),
                  disp.formula = list(~histologic2)
                ),
                SGOT = list(
                  long.formula = list(SGOT ~ time + I(time^2) + age + histologic2 + (1 + time + I(time^2)|id)),
                  family = list("gaussian")
                ),
                hepatomegaly = list(
                  long.formula = list(hepatomegaly ~ time + histologic2 + sex + (1 + time|id)),
                  family = list("binomial")
                ),
                platelets = list(
                  long.formula = list(platelets ~ ns(time, knots = c(1,4)) + age + histologic2 + (1 + ns(time, knots = c(1,4))|id)),
                  family = list("genpois"),
                  disp.formula = list(~time)
                ),
                prothrombin = list(
                  long.formula = list(prothrombin ~ time + histologic2 + sex + (1 + time|id)),
                  family = list("Gamma"),
                  disp.formula = list(~time)
                ),
                serBilir = list(
                  long.formula = list(serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex + (1 + ns(time, knots = c(1,4))|id)),
                  family = list("gaussian")
                ))
  environment(out$long.formula[[1]]) <- pf
  if(!is.null(out$disp.formula)) environment(out$disp.formula[[1]]) <- pf
  out
}


# Fit the joint model -----------------------------------------------------
fit.joint <- function(resp){
  formulae <- get.formulae(resp)
  cat(sprintf("Fitting %s\n", resp))
  if(!is.null(formulae$disp.formula)){
    out <- joint(long.formulas = formulae$long.formula,
                 surv.formula = surv.formula,
                 data = PBCredx,
                 family = formulae$family,
                 disp.formulas = formulae$disp.formula,
                 control = list(tol.rel = 5e-3))
  }else{
    out <- joint(long.formulas = formulae$long.formula,
                 surv.formula = surv.formula,
                 data = PBCredx,
                 family = formulae$family,
                 control = list(tol.rel = 5e-3))
  }
  return(out)
}

resp.list <- setNames(as.list(resps), resps)
joints <- lapply(resp.list, fit.joint)

# save(joints, file = save.dir.file.path("AllUnivs.RData", log.dir))

# Get some stuffs
(et <- sapply(joints, '[[', "elapsed.time"))
(survs <- lapply(joints, function(x) summary(x)$Sur))


# 
# # Loop over univariate joint models and print the survival sub-model summaries.
# jms <- ls()[sapply(ls(), function(i) eval(parse(text=paste0("class(",i,")"))))=="joint"]
# 
# # So all significant by themselves.
# for(j in jms){
#   cat(j, "\n")
#   print(eval(parse(text = paste0(
#     "summary(", j, ")$Surv"
#   ))))
#   cat(cli::rule())
#   cat("\n")
# }

