# Longitudinal fits for count data ----------------------------------------
rm(list=ls())
source('.Rprofile')
library(survival)
log.dir <- save.dir.file.path("Survs")
head(surv.data)
pf <- parent.frame()

file.name <- function(x) paste0("Surv_", gsub("\\~|\\+|\\s","",x), ".log")
make.and.fit <- function(svars = c("drug", "sex", "age", "histologic")){
  svars <- unname(sapply(svars, match.arg, svars)) # Checking for typos
  cc <- length(svars)
  form.rhs <- paste0("~ ", paste0(svars, collapse = " + "))
  fn <- paste0(cc, file.name(form.rhs))
  cat(sprintf("\nFitting:\n%s\n\n", form.rhs))
  form <- as.formula(paste0("Surv(survtime, status)", form.rhs), env = pf)
  ph <- coxph(form, data = surv.data, x = T)
  sink(file = save.dir.file.path(fn, log.dir))
  print(summary(ph))
  sink()
  cat(sprintf("Printed to %s\n", save.dir.file.path(fn, log.dir)))
  cat(cli::rule(line_col = .nice.orange, background_col = "mediumseagreen", line = 2))
}
# 4-variate
make.and.fit()

var3 <- combn(c("drug", "sex", "age", "histologic"), 3)
apply(var3,2,make.and.fit)
var2 <- combn(c("drug", "sex", "age", "histologic"), 2)
apply(var2,2,make.and.fit)
sapply(c("drug", "sex", "age", "histologic"),make.and.fit)

# Think a 3-variate sex/age/histologic is winner: all sig at 10% level
# C-statistic of .727.

can1 <- coxph(Surv(survtime,status)~sex+age+histologic,surv.data)
can2 <- coxph(Surv(survtime,status)~age+histologic,surv.data)
can3 <- coxph(Surv(survtime,status)~sex+histologic,surv.data)
can4 <- coxph(Surv(survtime,status)~sex+age,surv.data)

anova(can2,can1) # Implies sex not necessary
anova(can3,can1) # Implies age is necessary.
anova(can4,can1) # Implies hist. necessary.

# Okay nevermind, final model is
fin <- coxph(Surv(survtime, status)~age+histologic,surv.data)
summary(fin)
