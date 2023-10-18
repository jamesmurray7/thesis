# DO NOT run JMbayes2 in an R project/with a custom .Rprofile.
rstudioapi::executeCommand('closeProject')
library(JMbayes2)
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)
surv.data$histologic2 <- ifelse(surv.data$histologic %in% c("3", "4"), 1, 0)
serb <- lme(fixed = serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex,
            random = ~  ns(time, knots = c(1,4))|id , data = PBCredx, method = "ML",
            control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
alb <- lme(fixed = albumin ~ time + age + histologic2,
           random = ~ time |id , data = PBCredx, method = "ML",
           control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
pro <- mixed_model(prothrombin ~ time + histologic2 + sex,
                   data = PBCredx, random = ~time|id,
                   family = Gamma(link = "log"))
M <- list(pro, serb, alb)
SM <- coxph(Surv(survtime,status)~age+sex+histologic2,data=surv.data)

J <- jm(SM, M, time_var = "time", n_chains = 3L, 
        n_iter = 10000L, n_burnin = 2000L)
save(J, file = f)

M2 <- list(serb, alb)
J2 <- jm(SM, M2, time_var = "time", n_chains = 3L, 
         n_iter = 10000L, n_burnin = 2000L)
save(J2, file = gsub("Triv", "Biv", f))

M3 <- list(serb, alb)
SM.nosex <- coxph(Surv(survtime,status)~age+histologic2,data=surv.data)
J3 <- jm(SM.nosex, M3, time_var = "time", n_chains = 3L, 
         n_iter = 10000L, n_burnin = 2000L)
save(J3, file = gsub("Triv", "Biv2", f))

# This crashed my Ubuntu earlier ~~
rm(list=ls())
rstudioapi::executeCommand('closeProject')
library(joineRML)
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)
surv.data$histologic2 <- ifelse(surv.data$histologic %in% c("3", "4"), 1, 0)

# And once more fit this with joineRML
joint.biv.joineRML <- mjoint(
  formLongFixed = list(
    "serBilir" = serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex,
    "albumin" = albumin ~ time + age + histologic2
  ),
  formLongRandom = list(
    "serBilir" = ~ 1 + ns(time, knots = c(1,4))|id,
    "albumin" = ~time|id
  ),
  data = PBCredx,
  formSurv = Surv(survtime, status) ~ age + sex + histologic2,
  timeVar = 'time', 
  control = list(
    type = 'sobol', tol.em = 1e-2, tol2 = 5e-3, convCrit = "sas"
  )
)

save(joint.biv.joineRML, file = save.dir.file.path("bivjML.RData"))

# remove sex from surv.models
joint.biv2.joineRML <- mjoint(
  formLongFixed = list(
    "serBilir" = serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex,
    "albumin" = albumin ~ time + age + histologic2
  ),
  formLongRandom = list(
    "serBilir" = ~ 1 + ns(time, knots = c(1,4))|id,
    "albumin" = ~time|id
  ),
  data = PBCredx,
  formSurv = Surv(survtime, status) ~ age + histologic2,
  timeVar = 'time', 
  control = list(
    type = 'sobol', tol.em = 1e-2, tol2 = 5e-3, convCrit = "sas"
  )
)

save(joint.biv2.joineRML, file = save.dir.file.path("biv2jML.RData"))


