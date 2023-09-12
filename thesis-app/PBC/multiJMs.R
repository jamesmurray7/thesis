rm(list=ls())
source(".Rprofile")
source("helpers.R")
out.dir <- save.dir.file.path("Joint")
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

cc <- cond.ranefs(Gs, tune = .5) # .5 ~27% (good enough)
cc
plot(cc)

summary(Gs)
Ps <- plotMultiJM(Gs)

# I think put the survival ones in main-text
Ps$Surv
ggsave(save.dir.file.path("GaussiansSurv.png", out.dir),
       width = 140, height = 90, units = "mm")
# and the longitudinals in appendix.
Ps$Long
ggsave(save.dir.file.path("GaussiansLong.png", out.dir),
       width = 140, height = 120, units = "mm")

# Counts ------------------------------------------------------------------
Cs <- joint(
  long.formulas = list(
    alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + 
      (1 + ns(time, knots = c(1, 4))|id),
    platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + 
      (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("negbin", "genpois")
)

Cs
Ccc <- cond.ranefs(Cs, tune = 1)
Ccc
plot(Ccc)

summary(Cs)
Pc <- plotMultiJM(Cs)

Pc$Surv
ggsave(save.dir.file.path("CountsSurv.png", out.dir),
       width = 140, height = 90, units = "mm")
Pc$Long
ggsave(save.dir.file.path("CountsLong.png", out.dir),
       width = 140, height = 90, units = "mm")

# + DispModels ------------------------------------------------------------
C.Disps <- joint(
  long.formulas = list(
    alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + 
      (1 + ns(time, knots = c(1, 4))|id),
    platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + 
      (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("negbin", "genpois"),
  disp.formulas = list(~time+age+drug+histologic,
                       ~time+age+histologic+sex) # This model is much better, carry forwards.
)

C.Disps
C.Disps.cc <- cond.ranefs(C.Disps, tune = 1)
C.Disps.cc
plot(C.Disps.cc)

summary(C.Disps)
PC.Disps <- plotMultiJM(C.Disps, 6)

PC.Disps$Surv
ggsave(save.dir.file.path("DispCountsSurv.png", out.dir),
       width = 140, height = 90, units = "mm")
PC.Disps$Long
ggsave(save.dir.file.path("DispCountsLong.png", out.dir),
       width = 140, height = 180, units = "mm")


# Testing - all `cont`inuous ----------------------------------------------
Conts <- joint(
  long.formulas = list(
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    SGOT ~ age + histologic + drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
    prothrombin ~ sex + histologic + time + (1 + time|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian", "gaussian", "Gamma"),
  disp.formulas = list(~1,~1,~1,~time+drug+age+histologic+sex)
)

Conts.cc <- cond.ranefs(Conts, N = 1e3, tune = .5)
Conts.cc
plot(Conts.cc)

Pconts <- plotMultiJM(Conts) # -> And then we drop {SGOT, prothrombin}
Pconts$Surv
ggsave(save.dir.file.path("ContsSurv.png", out.dir),
       width = 140, height = 90, units = "mm")
Pconts$Long
ggsave(save.dir.file.path("ContsLong.png", out.dir),
       width = 140, height = 180, units = "mm")


# `Cont`inuous -> reduced now ---------------------------------------------
Conts.reduced <- joint(
  long.formulas = list(
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian"),
  control = list(tol.rel = 5e-2)
)

Conts.reduced.cc <- cond.ranefs(Conts.reduced, tune = 1)
Conts.reduced.cc
plot(Conts.reduced.cc)

Pconts.reduced <- plotMultiJM(Conts.reduced)
Pconts.reduced$Surv
ggsave(save.dir.file.path("ContsReducedSurv.png", out.dir),
       width = 140, height = 90, units = "mm")
Pconts.reduced$Long
ggsave(save.dir.file.path("ContsReducedLong.png", out.dir),
       width = 140, height = 90, units = "mm")


# "Final" Multivariate ----------------------------------------------------
fin <- joint(
  long.formulas = list(
    # Continuous k=1,2
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    # Counts k = 3,4
    alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    # Binary k = 5
    hepatomegaly ~ sex + histologic + time + I(time^2) + (1 + time + I(time^2)|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian", "negbin", "genpois", "binomial"),
  disp.formulas = list(~1, ~1, ~time+age+drug+histologic, ~time+age+histologic+sex, ~1)
)
save(fin, file = save.dir.file.path("fin.RData", out.dir))

# Without hepa?
fin2 <- joint(
  long.formulas = list(
    # Continuous k=1,2
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    # Counts k = 3,4
    alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
    platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian", "negbin", "genpois"),
  disp.formulas = list(~1, ~1, ~time+age+drug+histologic, ~time+age+histologic+sex)
) 
save(fin2, file = save.dir.file.path("fin2.RData", out.dir))
# Looks better!!!

fin2.cc <- cond.ranefs(fin2, tune = .425)
png(filename = save.dir.file.path("test.png", out.dir),  # same as earlier plots really...
    units = "mm", width = 400, height = 400, res = 500)  # which is good as too big!
plot(fin2.cc)
dev.off()

Pfin2 <- plotMultiJM(fin2, 10)
Pfin2$Surv
ggsave(save.dir.file.path("fin2Surv.png", out.dir),
       width = 140, height = 90, units = "mm")
Pfin2$Long
ggsave(save.dir.file.path("fin2Long.png", out.dir),
       width = 140, height = 185, units = "mm")


# Final 3 -----------------------------------------------------------------
fin.fin <- joint(
  long.formulas = list(
    # Continuous k=1,2
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian", "gaussian")
)

fin.fin2 <- joint( # strangely(?) this is half as fast!?
  long.formulas = list(
    # Continuous k=1,2
    albumin ~ age + histologic + time + I(time^2) + (1 + time + I(time^2)|id),
    serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id)
  ),
  surv.formula = Surv(survtime, status) ~ age,
  data = PBCredx,
  family = list("gaussian", "gaussian")
)



