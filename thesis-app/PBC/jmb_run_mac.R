# Do not run JMbayes2 in an R project.
rstudioapi::executeCommand('closeProject')
library(JMbayes2)
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)
surv.data$histologic2 <- ifelse(surv.data$histologic %in% c("3", "4"), 1, 0)
serb <- lme(fixed = serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex,
            random = ~  ns(time, knots = c(1,4))|id , data = PBCredx, method = "ML",
            control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
alb <- lme(fixed = serBilir ~ time + age + histologic2,
           random = ~ time |id , data = PBCredx, method = "ML",
           control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
pro <- mixed_model(prothrombin ~ time + histologic2 + sex,
                   data = PBCredx, random = ~time|id,
                   family = Gamma(link = "log"))
M <- list(pro, serb, alb)
SM <- coxph(Surv(survtime,status)~age+sex+histologic2,data=surv.data)

J <- jm(SM, M, time_var = "time", n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
