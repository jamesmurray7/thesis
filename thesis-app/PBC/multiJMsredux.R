# Multivariate joint models
# Two disparate approaches
# This file _simply creates_ .RData objects for different model fits.

rm(list=ls())
source(".Rprofile")
source("helpers.R")
out.dir <- save.dir.file.path("Joint/Multivs")
library(splines)
pf <- parent.frame()
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)

(resps <- names(PBCredx)[8:15])
resps <- resps[resps!="spiders"]

# Function to fit all based on best model
surv.formula <- Surv(survtime, status) ~ age + sex + histologic2

# Longit. -----------------------------------------------------------------
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
  if(is.null(out$disp.formula)) out$disp.formula <- list(~1)
  environment(out$disp.formula[[1]]) <- pf
  out
}


# Fit all at once first ---------------------------------------------------
# .all <- lapply(resps, get.formulae)
# long.formulas <- sapply(.all, function(x){
#   x$long.formula
# })
# disp.formulas <- sapply(.all, function(x){
#   x$disp
# })
# family <- sapply(.all, '[[', 2)
# 
# OneBigModel <- joint(long.formulas = long.formulas, surv.formula = surv.formula,
#                      data = PBCredx, disp.formulas = disp.formulas,
#                      family = family, control = list(tol.rel = 5e-3))
# save(OneBigModel, file = save.dir.file.path("OneBigModel.RData", out.dir))

# This doesn't look good - not enough events for num. parameters
# But this is to be expected (see THESIS!)


# Type-specific ---------------------------------------------------------
# Continuous
conts <- c("SGOT", "serBilir", "albumin", "prothrombin")
.conts <- lapply(conts, get.formulae)
long.formulas.conts <- sapply(.conts, function(x){
  x$long.formula
})
disp.formulas.conts <- sapply(.conts, function(x){
  x$disp
})
family.conts <- sapply(.conts, '[[', 2)

joint.conts <- joint(
  long.formulas = long.formulas.conts, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.conts,
  family = family.conts, control = list(tol.rel = 5e-3)
)

save(joint.conts, file = save.dir.file.path("ContsOnly.RData", out.dir))

# Counts
counts <- c("platelets", "alkaline")
.counts <- lapply(counts, get.formulae)
long.formulas.counts <- sapply(.counts, function(x){
  x$long.formula
})
disp.formulas.counts <- sapply(.counts, function(x){
  x$disp
})
family.counts <- sapply(.counts, '[[', 2)

joint.counts <- joint(
  long.formulas = long.formulas.counts, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.counts,
  family = family.counts, control = list(tol.rel = 5e-3)
)

save(joint.counts, file = save.dir.file.path("CountsOnly.RData", out.dir))

# Bin
joint.bin <- joint(
  long.formulas = list(hepatomegaly ~ time + histologic2 + sex + (1 + time|id)),
  surv.formula = surv.formula,
  data = PBCredx, family = list("binomial"),
  control = list(tol.rel = 5e-3)
)

save(joint.bin, file = save.dir.file.path("BinOnly.RData", out.dir))

# Function-specific -------------------------------------------------------
# Liver Enzymes: ----
# Alkaline (Alkaline phosphatase)
# SGOT (Serum Glutamic Oxaloacetic Transaminase, or AST)
# Serum Bilirubin (SerBilir)

enzymes <- c("alkaline", "SGOT", "serBilir")

enzymes <- lapply(enzymes, get.formulae)
long.formulas.enzymes <- sapply(enzymes, function(x){
  x$long.formula
})
disp.formulas.enzymes <- sapply(enzymes, function(x){
  x$disp
})
family.enzymes <- sapply(enzymes, '[[', 2)

joint.enzymes <- joint(
  long.formulas = long.formulas.enzymes, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.enzymes,
  family = family.enzymes, control = list(tol.rel = 5e-3)
)

save(joint.enzymes, file = save.dir.file.path("enzymes.RData", out.dir))

# Blood Clotting and Coagulation: ----
# Prothrombin (Prothrombin time)
# Platelets

bloods <- c("prothrombin", "platelets")
bloods <- lapply(bloods, get.formulae)
long.formulas.bloods <- sapply(bloods, function(x){
  x$long.formula
})
disp.formulas.bloods <- sapply(bloods, function(x){
  x$disp
})
family.bloods <- sapply(bloods, '[[', 2)

joint.bloods <- joint(
  long.formulas = long.formulas.bloods, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.bloods,
  family = family.bloods, control = list(tol.rel = 5e-3)
)

save(joint.bloods, file = save.dir.file.path("bloods.RData", out.dir))

# Liver Health and Function: ----
# Hepatomegaly (Enlarged liver)
# Albumin (Liver-produced protein)

health <- c("hepatomegaly", "albumin")
health <- lapply(health, get.formulae)
long.formulas.health <- sapply(health, function(x){
  x$long.formula
})
disp.formulas.health <- sapply(health, function(x){
  x$disp
})
family.health <- sapply(health, '[[', 2)

joint.health <- joint(
  long.formulas = long.formulas.health, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.health,
  family = family.health, control = list(tol.rel = 5e-3)
)

save(joint.health, file = save.dir.file.path("health.RData", out.dir))


# Next model stage --------------------------------------------------------
triv <- c("prothrombin", "serBilir", "albumin")
triv <- lapply(triv, get.formulae)
long.formulas.triv <- sapply(triv, function(x){
  x$long.formula
})
disp.formulas.triv <- sapply(triv, function(x){
  x$disp
})
family.triv <- sapply(triv, '[[', 2)

joint.triv <- joint(
  long.formulas = long.formulas.triv, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.triv,
  family = family.triv, control = list(tol.rel = 5e-3)
)

save(joint.triv, file = save.dir.file.path("triv.RData", out.dir))


# Next (dropping proth) ---------------------------------------------------
biv <- c("serBilir", "albumin")
biv <- lapply(biv, get.formulae)
long.formulas.biv <- sapply(biv, function(x){
  x$long.formula
})
disp.formulas.biv <- sapply(biv, function(x){
  x$disp
})
family.biv <- sapply(biv, '[[', 2)

joint.biv <- joint(
  long.formulas = long.formulas.biv, surv.formula = surv.formula,
  data = PBCredx, disp.formulas = disp.formulas.biv,
  family = family.biv, control = list(tol.rel = 5e-3)
)

save(joint.biv, file = save.dir.file.path("biv.RData", out.dir))

# Dropping zeta_sex (this is the final model) -----------------------------
joint.biv2 <- joint(
  long.formulas = long.formulas.biv, 
  surv.formula = Surv(survtime, status) ~ age + histologic2,
  data = PBCredx, disp.formulas = disp.formulas.biv,
  family = family.biv, control = list(tol.rel = 5e-3)
)

save(joint.biv2, file = save.dir.file.path("biv2.RData", out.dir))



