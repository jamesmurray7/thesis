# Joint models for count data ---------------------------------------------
rm(list=ls())
source('.Rprofile')
truncated <- F
log.dir <- if(truncated) save.dir.file.path("logJM/truncate_at_3") else save.dir.file.path("logJM/not_truncated")
if(truncated) adni <- adni[adni$time <= 3, ]
# adni[adni]
library(gmvjoint)


# Some functions to create data -------------------------------------------
cleandata <- function(resp){
  sset <- adni[!is.na(adni[, resp]), ]
  u.rid <- unique(sset$RID)
  sset$id <- NULL
  u.df <- data.frame(RID = u.rid, id = 1:length(u.rid))
  sset2 <- merge(sset, u.df, "RID")
  if(truncated){
    sset2$status <- ifelse(sset2$survtime >= 3.1, 0, sset2$status)
    sset2$survtime <- ifelse(sset2$status == 0, 3.1, sset2$survtime)
  }
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


# ADAS11 ------------------------------------------------------------------
fit.Po.ADAS11.int <- fit.JM("ADAS11", "int", "poisson", ~0, Surv(survtime, status) ~ APOE4)
fit.Po.ADAS11.lin <- fit.JM("ADAS11", "linear", "poisson", ~0, Surv(survtime, status) ~ APOE4)
fit.GP.ADAS11.int <- fit.JM("ADAS11", "int", "genpois", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.ADAS11.spl <- fit.JM("ADAS11", "ns", "genpois", ~1, Surv(survtime, status) ~ APOE4)
# Poissons
summary(fit.Po.ADAS11.int)
plot(resid(fit.Po.ADAS11.int, "longit", "pearson"))
summary(fit.Po.ADAS11.lin)
plot(resid(fit.Po.ADAS11.lin, "longit", "pearson"))
# Genpois's
summary(fit.GP.ADAS11.int)
plot(resid(fit.GP.ADAS11.int, "longit", "pearson"))
summary(fit.GP.ADAS11.spl)
plot(resid(fit.GP.ADAS11.spl, "longit", "pearson")) # Fits more outliers!

# Poisson -> linear vs. genpois -> natural spline

# ADAS13 ------------------------------------------------------------------
fit.PO.ADAS13.int <- fit.JM("ADAS13", "int", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.PO.ADAS13.lin <- fit.JM("ADAS13", "linear", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.ADAS13.int <- fit.JM("ADAS13", "int", "genpois", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.ADAS13.lin <- fit.JM("ADAS13", "linear", "genpois", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.ADAS13.lin.time <- fit.JM("ADAS13", "linear", "genpois", ~time, Surv(survtime, status) ~ APOE4)

summary(fit.PO.ADAS13.int)
summary(fit.PO.ADAS13.lin)
summary(fit.GP.ADAS13.int)
summary(fit.GP.ADAS13.lin)
summary(fit.GP.ADAS13.lin.time)

# Poisson
plot(resid(fit.PO.ADAS13.int, type = "pearson"))
plot(resid(fit.PO.ADAS13.lin, type = "pearson")) # better
# Genpois's
plot(resid(fit.GP.ADAS13.int, type = "pearson"))
plot(resid(fit.GP.ADAS13.lin, type = "pearson"))
plot(resid(fit.GP.ADAS13.lin.time, type = "pearson"))

# Poisson -> linear vs. genpois -> linear + some dispersion model?

# MMSE --------------------------------------------------------------------
adni$MMSE2 <- 30-adni$MMSE
fit.Po.MM.int <- fit.JM("MMSE2", "int", "poisson", ~1, Surv(survtime, status) ~ APOE4)
summary(fit.Po.MM.int)
fit.Po.MM.lin <- fit.JM("MMSE2", "linear", "poisson", ~1, Surv(survtime, status) ~ APOE4)
summary(fit.Po.MM.lin)
fit.GP.MM.int <- fit.JM("MMSE2", "int", "genpois", ~1, Surv(survtime, status) ~ APOE4)
summary(fit.GP.MM.int)

# Poisson
plot(resid(fit.Po.MM.int, type = "pearson"))
plot(resid(fit.Po.MM.lin, type = "pearson")) # Maybe linear poisson?
# GP
plot(resid(fit.GP.MM.int, type = "pearson")) # Though residuals less problematic in GP case?


# FAQ ---------------------------------------------------------------------
fit.Po.FA.int <- fit.JM("FAQ", "int", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.Po.FA.lin <- fit.JM("FAQ", "linear", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.Po.FA.qua <- fit.JM("FAQ", "quad", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.Po.FA.spl <- fit.JM("FAQ", "ns", "poisson", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.FA.int <- fit.JM("FAQ", "int", "genpois", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.FA.lin <- fit.JM("FAQ", "linear", "genpois", ~1, Surv(survtime, status) ~ APOE4)
fit.GP.FA.qua <- fit.JM("FAQ", "quad", "genpois", ~1, Surv(survtime, status) ~ APOE4)

# Poisson
summary(fit.Po.FA.int)
summary(fit.Po.FA.lin)
summary(fit.Po.FA.qua)
summary(fit.Po.FA.spl)
AIC(fit.Po.FA.int); AIC(fit.Po.FA.lin); AIC(fit.Po.FA.qua); AIC(fit.Po.FA.spl)
plot(resid(fit.Po.FA.int, type = "pearson"))
plot(resid(fit.Po.FA.lin, type = "pearson")) # Anything but intercept model...
plot(resid(fit.Po.FA.qua, type = "pearson"))
plot(resid(fit.Po.FA.spl, type = "pearson"))

summary(fit.GP.FA.int)
summary(fit.GP.FA.lin)
summary(fit.GP.FA.qua) # Evidence to suggest quadratic fit best (overdispersed)
plot(resid(fit.GP.FA.int, type = "pearson"))
plot(resid(fit.GP.FA.lin, type = "pearson")) 
plot(resid(fit.GP.FA.qua, type = "pearson")) # Quadratic fit enforced?
