# Dispersion models for count data ----------------------------------------
rm(list=ls())
source('.Rprofile')
truncated <- T
log.dir <- if(truncated) save.dir.file.path("logJM/truncate_at_3") else save.dir.file.path("logJM/not_truncated")
if(truncated) adni <- adni[adni$time <= 3, ]
library(gmvjoint)


# Some functions to create data -------------------------------------------
cleandata <- function(resp){
  sset <- adni[!is.na(adni[, resp]), ]
  u.rid <- unique(sset$RID)
  sset$id <- NULL
  u.df <- data.frame(RID = u.rid, id = 1:length(u.rid))
  sset2 <- merge(sset, u.df, "RID")
  sset2$status <- ifelse(sset2$survtime >= 3.1, 0, sset2$status)
  sset2$survtime <- ifelse(sset2$status == 0, 3.1, sset2$survtime)
  sset2
}

make.survdata <- function(data) data[!duplicated(data$id)] # Actually not needed

# Function to fit UNIV JM -------------------------------------------------
pf <- parent.frame()
fit.JM <- function(resp, time.spec, family, disp.formula, Surv.formula, ...){
  dd <- cleandata(resp)
  extras <- paste0(c("age_scaled", "gender"), collapse = ' + ')
  time <- switch(time.spec,
                 linear = " APOE4 * time + (1 + time|id)",
                 int = " APOE4 * time + (1|id)",
                 quad = " APOE4 * (time + I(time^2)) + (1 + time + I(time^2)|id)",
                 ns = " APOE4 * splines::ns(time, knots = c(0.5, 1.5)) + (1 + splines::ns(time, knots = c(0.5, 1.5))|id)")
  
  form <- if(nchar(extras > 1)) paste0(resp, " ~ ", extras, ' +', time) else paste0(resp, " ~ ", time)
  form <- as.formula(form)
  environment(form) <- pf
  long.formulas <- list(form); 
  disp.formulas <- lapply(1, function(i){
    a <- disp.formula  
    environment(a) <- pf
    a
  })
  
  fit <- joint(long.formulas, Surv.formula, data = dd, family = list(family), 
               disp.formulas = disp.formulas, control = list(tol.rel = 5e-3, ...))
  return(fit)
}

fit.Po.ADAS11 <- fit.JM("ADAS11", "int", "poisson", ~0, Surv(survtime, status) ~ APOE4)
summary(fit.Po.ADAS11)
plot(resid(fit.Po.ADAS11, "longit", "pearson"))
fit.GP.ADAS11 <- fit.JM("ADAS11", "ns", "genpois", ~1, Surv(survtime, status) ~ APOE4, verbose = T)
