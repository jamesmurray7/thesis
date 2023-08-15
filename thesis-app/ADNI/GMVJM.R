# generalised multivariate joint models for count data --------------------
rm(list=ls())
source('.Rprofile')
truncated <- F
# log.dir <- if(truncated) save.dir.file.path("logJM/truncate_at_3") else save.dir.file.path("logJM/not_truncated")
if(truncated) adni <- adni[adni$time <= 3, ]
# adni[adni]
library(gmvjoint)

# Some functions to create data -------------------------------------------
cleandata <- function(resps){
  sset <- na.omit(adni[,c("RID","id", "APOE4", "survtime", "status", "age_scaled",
                          "gender", "time", resps)])
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

# MVs ---------------------------------------------------------------------
pf <- parent.frame()
fit.JM <- function(resps, time.specs, disps, families, Surv.formula = NULL, ...){
  if(!is.list(families)) family <- as.list(families) else family <- families
  extras <- paste0(c("age_scaled", "gender"), collapse = ' + ')
  
  long.formulas <- lapply(seq_along(resps), function(k){
    ts <- time.specs[k]
    time <- switch(ts,
                   linear = " APOE4 * time + (1 + time|id)",
                   int = " APOE4 * time + (1|id)",
                   quad = " APOE4 * (time + I(time^2)) + (1 + time + I(time^2)|id)",
                   ns = " APOE4 * splines::ns(time, knots = c(0.5, 1.5)) + (1 + splines::ns(time, knots = c(0.5, 1.5))|id)")
    form <- paste0(resps[k], " ~ ", extras, ' +', time)
    form <- as.formula(form)
    environment(form) <- pf
    form
  })
  
  disp.formulas <- as.list(disps)
  if(is.null(Surv.formula)) 
    surv.formula <- Surv(survtime, status) ~ APOE4 
  else 
    surv.formula <- Surv.formula
  
  dd <- cleandata(resps)
  
  joint(long.formulas = long.formulas, 
        surv.formula = surv.formula, data = dd, family = family, 
        disp.formulas = disp.formulas, 
        control = list(tol.rel = 5e-3, ...))
}

test1 <- fit.JM(
  resps = c("ADAS13", "FAQ"),
  time.specs = c("linear", "quad"), disps = c(~1,~1), families = c("poisson", "poisson"),
)
summary(test1)
plot(resid(test1, type = "pearson"))
plot(resid(test1, what = "surv")) # Good e.g.
cr.test1 <- cond.ranefs(test1)    
plot(cr.test1)  # Appears reasonable. plot(cr.test1, 1L)

