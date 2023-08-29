log.dir <- save.dir.file.path("Longs")
resp.dirs <- dir(log.dir)
# Establish directories containing .log files
dirs.to.parse <- sapply(resp.dirs, save.dir.file.path, log.dir)

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
