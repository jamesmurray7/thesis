rm(list=ls())
source(".Rprofile")
file.show("printout.txt", title = "models for reference")
log.dir <- save.dir.file.path("Disps")
resp.dirs <- dir(log.dir)
# Establish directories containing .log files
dirs.to.parse <- sapply(resp.dirs, save.dir.file.path, log.dir)
dirs.to.parse <- dirs.to.parse[!grepl("\\.log", dirs.to.parse)]


# Create .log dumps -------------------------------------------------------
#           (summaries)
for(d in dirs.to.parse){
  d.files <- dir(d, pattern = "\\.log")
  # get response
  resp <- gsub(paste0(save.dir,"\\/Disps\\/"), "", d)
  if(length(d.files) == 0){
    cat(sprintf("Nothing present in %s currently!\n", d))
    next
  }
  d.collect <- structure(matrix(NA, 0, 5),
                         dimnames=list(c(),c("AIC", "BIC", "logLik", "deviance", "df.resid") ))
  n <- 0
  for(f in d.files){
    ff <- save.dir.file.path(f, d)
    log <- readLines(ff)
    ind <- grep("AIC", log) + 1
    fixed <- gsub("\\s+", " ", trimws(log[ind]))
    if(grepl("NA", fixed)){
      cat(sprintf("\nNA values in %s, skipping...\n", f))
      next
    }
    n <- n + 1
    fixed <- setNames(as.numeric(el(strsplit(fixed,"\\s"))),
                      c("AIC", "BIC", "logLik", "deviance", "df.resid"))
    d.collect <- rbind(d.collect, f=fixed)
    row.names(d.collect)[n] <- f
    rm(fixed)
  }
  
  cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by AIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(d.collect[,"AIC"]),][1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by BIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(d.collect[,"BIC"]),][1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by logLik"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(-d.collect[,"logLik"]),][1:7,]);cat("\n")
  
  # And now sink 
  sink(file = save.dir.file.path(paste0(resp,".log"), log.dir))
  cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by AIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(d.collect[,"AIC"]),][1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by BIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(d.collect[,"BIC"]),][1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by logLik"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(d.collect[order(-d.collect[,"logLik"]),][1:7,]);cat("\n")
  sink()

  fs::file_copy(save.dir.file.path(paste0(resp,".log"), log.dir),
                paste0('./summarylogs/DISP_', resp, ".log"),
                TRUE)
  
}


# Verify via ANOVA --------------------------------------------------------

# alkaline -> ====
intonly <- glmmTMB(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
                   data = PBCredx, family = glmmTMB::nbinom2())
# Three best by BIC
m1 <- glmmTMB(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::nbinom2(), dispformula = ~time + histologic)
m2 <- glmmTMB(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::nbinom2(), dispformula = ~time + drug + histologic)
m3 <- glmmTMB(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::nbinom2(), dispformula = ~time + drug + age + histologic) # Best by AIC
# 2nd best by AIC
m4 <- glmmTMB(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::nbinom2(), dispformula = ~time + drug + age + histologic + sex) # Also best by lL.

anova(intonly, m1)
anova(m1, m2) # m2 better
anova(m2, m3) # m3 better
anova(m3, m4) # i.e. m3 `best` ->>  time + drug + age + histologic.
rm(intonly, m1, m2, m3, m4)

# Platelets -> ====
intonly <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
                   data = PBCredx, family = glmmTMB::genpois())
# Three best by BIC
m1 <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time)
m2 <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time + age)
m3 <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time + sex)
# Three best by AIC
m4 <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time + age + histologic + sex)
m5 <- glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
              data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time + age + histologic + sex + drug) # best by lL
m6 <-glmmTMB(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
             data = PBCredx, family = glmmTMB::genpois(), dispformula = ~time + age + histologic)

anova(intonly, m1)
anova(m1, m2)
anova(m2, m3)
anova(m2, m4)
anova(m4, m5)
anova(m4, m6) # i.e. m4 `best` ->> time + age + sex + histologic.

rm(intonly, m1, m2, m3, m4, m5, m6)

# Prothrombin -> ====
intonly <- glmmTMB(prothrombin ~ age + histologic + drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
                   data = PBCredx, family = Gamma(link = "log"))
# Best by AIC, BIC and lL.
m1 <- glmmTMB(prothrombin ~ age + histologic + drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
              data = PBCredx, family = Gamma(link = "log"), dispformula = ~time + drug + age + histologic + sex)
# Drop drug?
m2 <- glmmTMB(prothrombin ~ age + histologic + drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
              data = PBCredx, family = Gamma(link = "log"), dispformula = ~time + age + histologic + sex)
anova(m1, m2) #i.e. m1 `best` ->> time + drug + age + histologic + sex
rm(intonly, m1, m2)
