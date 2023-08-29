log.dir <- save.dir.file.path("Longs")
resp.dirs <- dir(log.dir)
# Establish directories containing .log files
dirs.to.parse <- sapply(resp.dirs, save.dir.file.path, log.dir)
dirs.to.parse <- dirs.to.parse[!grepl("\\.log", dirs.to.parse)]


# Create .log dumps -------------------------------------------------------
#           (summaries)
for(d in dirs.to.parse){
  d.files <- dir(d, pattern = "\\.log")
  # get response
  resp <- gsub(paste0(save.dir,"\\/Longs\\/"), "", d)
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
                paste0('./summarylogs/', resp, ".log"),
                TRUE)
  
}


# Verify via ANOVA --------------------------------------------------------
# (i think BIC ranking is producing the most sensible-sounding) `best` model
#  so use this metric first...) and then compare with best AIC/ll

# For bivs: 
# // [[1]]: sex + age; [[2]]: sex + histologic; [[3]]: age + histologic // 

# Albumin -> ====
# top 2 by BIC
m1 <- fit.glmmTMB(albumingaussian$linearbiv[[3]]$longformula, albumingaussian$linearbiv[[3]]$family)
m2 <- fit.glmmTMB(albumingaussian$linearbiv[[2]]$longformula, albumingaussian$linearbiv[[2]]$family)
# top by AIC
m3 <- fit.glmmTMB(albumingaussian$quadbiv[[3]]$longformula,albumingaussian$quadbiv[[1]]$family)
# top by lL
m4 <- fit.glmmTMB(albumingaussian$quadtriv$longformula, albumingaussian$quadtriv$family)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4)
anova(m3, m4) # i.e. m3 `best` ->> age+hist+quad time spec.

# serBilir -> ====
# (quite a lot of agreement between metrics here)
# Top by AIC/BIC
m1 <- fit.glmmTMB(serBilirgaussian$nsbiv[[2]]$longformula, gaussian)
m2 <- fit.glmmTMB(serBilirgaussian$nsbiv[[3]]$longformula, gaussian)
m3 <- fit.glmmTMB(serBilirgaussian$nstriv$longformula,gaussian)
# top by lL
m4 <- fit.glmmTMB(serBilirgaussian$nsDItriv$longformula, gaussian)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4) # i.e. m1 `best` ->> sex + histologic + cubic ns time spec.

# SGOT -> ====
# Three best by BIC
m1 <- fit.glmmTMB(SGOTgaussian$quadbiv[[3]]$longformula, gaussian)
m2 <- fit.glmmTMB(SGOTgaussian$quadtriv$longformula, gaussian) # (this is 3rd best by AIC)
m3 <- fit.glmmTMB(SGOTgaussian$linearbiv[[3]]$longformula, gaussian)
# Three best by AIC
m4 <- fit.glmmTMB(SGOTgaussian$quadDItriv$longformula, gaussian) # nb this also best by lL
m5 <- fit.glmmTMB(SGOTgaussian$quadDIbiv[[3]]$longformula, gaussian)

anova(m1, m2)
anova(m1, m3) # take m1 (best by AIC)
anova(m1, m4) # Take m4 (sex inclusion and drug interaction is sig.)
anova(m4, m5) # Take m5, sex is barely sig after drug interaction accounted for.
              # i.e. m5 `best` ->> age + histologic + drug * time^2
# spiders -> ====
# Best by BIC and AIC
m1 <- fit.glmmTMB(spidersbinomial$quadbiv[[3]]$longformula, binomial)
m2 <- fit.glmmTMB(spidersbinomial$quadbiv[[2]]$longformula, binomial)
m3 <- fit.glmmTMB(spidersbinomial$quadtriv$longformula, binomial)
# Best by lL
m4 <- fit.glmmTMB(spidersbinomial$quadDItriv$longformula, binomial)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4) # i.e. m1 `best` ->> age + histologic + quadratic time spec.

# hepatomegaly -> ====
# Three best by BIC
m1 <- fit.glmmTMB(hepatomegalybinomial$linearbiv[[2]]$longformula, binomial)
m2 <- fit.glmmTMB(hepatomegalybinomial$linearbiv[[3]]$longformula, binomial)
m3 <- fit.glmmTMB(hepatomegalybinomial$lineartriv$longformula, binomial)
# Three best by AIC
m4 <- fit.glmmTMB(hepatomegalybinomial$quadbiv[[2]]$longformula, binomial)
m5 <- fit.glmmTMB(hepatomegalybinomial$quadtriv$longformula, binomial)
m6 <- fit.glmmTMB(hepatomegalybinomial$quadDIbiv[[2]]$longformula, binomial)
# Best by lL
m7 <- fit.glmmTMB(hepatomegalybinomial$quadDItriv$longformula, binomial)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4) # take m4 (i.e. quadratic is better)
anova(m4, m5)
anova(m4, m6)
anova(m4, m7) #i.e. m4 `best` ->> sex + histologic + quadratic time spec.

# prothrombin -> ====
# BIC
m1 <- fit.glmmTMB(prothrombinGamma$quadDIbiv[[3]]$longformula, prothrombinGamma$quadDIbiv[[3]]$family) # This top by AIC and lL
m2 <- fit.glmmTMB(prothrombinGamma$linearbiv[[2]]$longformula, prothrombinGamma$linearbiv[[2]]$family) # 2nd by AIC too
m3 <- fit.glmmTMB(prothrombinGamma$linearbiv[[3]]$longformula, prothrombinGamma$linearbiv[[3]]$family) # 5th by AIC.
# AIC
m4 <- fit.glmmTMB(prothrombinGamma$lineartriv$longformula, prothrombinGamma$lineartriv$family)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4) # i.e. m1 `best` ->> age + histologic + drug * time^2.

# platelets -> ====
# BIC
m1 <- fit.glmmTMB(plateletsgenpois$nsbiv[[2]]$longformula, genpois) # Number 2 by AIC
m2 <- fit.glmmTMB(plateletsgenpois$nsbiv[[3]]$longformula, genpois)
m3 <- fit.glmmTMB(plateletsgenpois$nstriv$longformula, genpois)
# AIC
m4 <- fit.glmmTMB(plateletsgenpois$nsDIbiv[[2]]$longformula, genpois)
m5 <- fit.glmmTMB(plateletsgenpois$nsDItriv$longformula, genpois) # Number 1 by lL

anova(m1, m2) 
anova(m1, m3)
anova(m1, m4) # m4 better
anova(m4, m5) # i.e. m4 `best` ->> sex + histologic + drug * cubic splines

# alkaline -> ====
# BIC
m1 <- fit.glmmTMB(alkalinenegbin$nsbiv[[3]]$longformula, alkalinenegbin$nsbiv[[3]]$family) # first by AIC
m2 <- fit.glmmTMB(alkalinenegbin$nsbiv[[1]]$longformula, alkalinenegbin$nsbiv[[1]]$family) 
m3 <- fit.glmmTMB(alkalinenegbin$nstriv$longformula, alkalinenegbin$nstriv$family)         # second by AIC
# AIC
m4 <- fit.glmmTMB(alkalinenegbin$nsDIbiv[[3]]$longformula, alkalinenegbin$nsDIbiv[[3]]$family)
# lL
m5 <- fit.glmmTMB(alkalinenegbin$nsDItriv$longformula, alkalinenegbin$nsDItriv$family)

anova(m1, m2)
anova(m1, m3)
anova(m1, m4)
anova(m1, m5) # i.e. m1 `best` ->>  age + histologic + cubic splines
