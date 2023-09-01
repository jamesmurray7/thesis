rm(list=ls())
source(".Rprofile")
library(splines)


# Gaussians ---------------------------------------------------------------
Gs <- joint(
  long.formulas = list(
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    SGOT ~ age + histologic + drug * (time + I(time^2)) + (1 + time + I(time^2)|id)
    
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian", "gaussian")
)

cc <- cond.ranefs(Gs)

summary(Gs)
