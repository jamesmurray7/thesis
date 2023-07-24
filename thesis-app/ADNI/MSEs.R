# Finding Bootstrapped MSEs for B model fits ------------------------------
# (i.e. across competing families, as means for direct comparison...)

# Think we need to do a mixture of BOTH:
#  i) Fit two models FIRST, (for point-estimate?)
#  ii) Bootstrap from these repeatedly, in each case getting bootstrapped e
#      estimate.
#  iii) return a neat list.

source('.Rprofile')
truncated <- FALSE
pf <- parent.frame()
cleandata <- function(resps){
  sset <- na.omit(adni[,c("RID","id", "APOE4", "survtime", "status", "age_scaled",
                          "gender", "time", resps)])
  sc <- lapply(resps, function(r) c(scale(sset[,r])))
  sc <- do.call(cbind, sc)
  colnames(sc) <- paste0(resps, "_scaled")
  sset <- cbind(sset, sc)
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

make.formulae <- function(resp, family, time.spec = NULL, extras = NULL, 
                          surv.formula = NULL, disp.formula = NULL){
  # Set out defaults
  if(!is.null(time.spec)) time.spec <- time.spec else time.spec <- "linear"
  if(!is.null(extras)) extras <- extras else extras <- c("age_scaled", "gender")
  if(!is.null(surv.formula)) surv.formula <- surv.formula else surv.formula <- Surv(survtime, status) ~ APOE4
  if(!is.null(disp.formula)) disp.formula <- disp.formula else disp.formula <- ~1 
  if(family=="gaussianlog"){
    flag <- 1
    family <- "gaussian"
  }else if(family == "gaussianscaled"){
    flag <- 2
    family <- "gaussian"
  }else if(family == "gaussianraw"){
    flag <- 3
    family <- "gaussian"
  }else{
    family <- family; flag <- 99
  }
  family <- as.list(family)
  # Longitudinal part
  extras <- paste0(extras, collapse = ' + ')
  time <- switch(time.spec,
                 linear = " APOE4 * time + (1 + time|id)",
                 int = " APOE4 * time + (1|id)",
                 quad = " APOE4 * (time + I(time^2)) + (1 + time + I(time^2)|id)",
                 ns = " APOE4 * splines::ns(time, knots = c(0.5, 1.5)) + (1 + splines::ns(time, knots = c(0.5, 1.5))|id)")
  resp <- if(family[[1]]=="gaussian" && flag == 2) paste0(resp, "_scaled") else resp     # 2 -> scaled. (Gaussian only)
  resp <- if(family[[1]]=="gaussian" && flag == 1) paste0("log(", resp, ")") else resp   # 1 -> logged. (Gaussian only)
  
  long.formula <- paste0(resp, " ~ ", extras, ' +', time)
  long.formula <- as.formula(long.formula)
  environment(long.formula) <- pf
  long.formula <- list(long.formula)
  
  environment(disp.formula) <- pf
  disp.formula <- list(disp.formula)
  environment(surv.formula) <- pf
  
  list(long = long.formula, disp = disp.formula, family = family, Surv = surv.formula)
}

fit.univ.JM <- function(data, Flist, ...){
  x <- Flist
  return(joint(
    long.formulas = x$long, surv.formula = x$Surv,
    data = data, family = x$family, disp.formulas = x$disp,
    control = list(...)
  ))
}

# Function to resample data ----
resampledata <- function(data){
  uids <- unique(data$id)
  samps <- sample(x = uids, size = length(uids), replace = TRUE)
  newData <- setNames(lapply(1:length(uids), function(i){
    newData <- data[data$id == samps[i],]
    newData$InternalKey <- i
    newData$..old.id <- newData$id
    newData$id <- i
    newData
  }),paste0('original id ', samps)) # overkill but we dont look at this anyway.
  
  as.data.frame(do.call(rbind, newData))
}

.i <- function(x) if(!is.null(x)) return(x) else return(NULL)
boot.MSE <- function(resp, competing = list(), B = 50L){
  cat(cli::rule(line_col = 'orange', center = cli::col_blue(sprintf("Starting %s", resp))))
  cat("\n")
  parent.data <- cleandata(resp)
  # List of formulae for models (1) and (2)
  Flists <- lapply(competing, function(x){
    make.formulae(resp, x$family, .i(x$time.spec), .i(x$extras), .i(x$surv.formula), .i(x$disp.formula))
  })
  # Fit univariate models to the actual data...
  mod1 <- fit.univ.JM(parent.data, Flists[[1]], tol.rel = 5e-3)
  N <- unname(mod1$ModelInfo$nobs[1])
  MSE1 <- sum(resid(mod1)[[1]]^2)/N
  ll1 <- logLik(mod1)
  AICtab1 <- c("logLik" = c(ll1), "AIC" = attr(ll1, "AIC"), "BIC" = attr(ll1, "BIC"))
  cli::cli_alert_success("Model 1 ({Flists[[1]]$family}) fit done.")
  mod2 <- fit.univ.JM(parent.data, Flists[[2]], tol.rel = 5e-3)
  MSE2 <- sum(resid(mod2)[[1]]^2)/N
  ll2 <- logLik(mod2)
  AICtab2 <- c("logLik" = c(ll2), "AIC" = attr(ll2, "AIC"), "BIC" = attr(ll2, "BIC"))
  cli::cli_alert_success("Model 2 ({Flists[[2]]$family}) fit done.")
  cat("\n")
  cli::cli_alert_info("Starting {B} bootstrap fits...")
  pb <- utils::txtProgressBar(max = B, style = 3)
  
  AICtabs.1 <- AICtabs.2 <- structure(matrix(0, nrow = 3, ncol = B),
                                      dimnames = list(c("logLik", "AIC", "BIC"),
                                                      paste0("b",1:B)))
  MSEs <- structure(matrix(0, nrow = 2, ncol = B),
                    dimnames = list(c("mod1", "mod2"), 
                                    paste0("b",1:B)))
  for(b in 1:B){
    boot.data <- resampledata(parent.data)
    # Model 1
    mod1.b <- fit.univ.JM(boot.data, Flists[[1]], tol.rel = 1e-2)
    resid1.b <- resid(mod1.b, type = "pearson")[[1]]
    MSE1.b <- sum(resid1.b^2)/N
    ll1.b <- logLik(mod1.b)
    AICtabs.1[,b] <- c(c(ll1.b), attr(ll1.b, 'AIC'), attr(ll1.b, "BIC"))
    rm(mod1.b)
    # Model 2
    mod2.b <- fit.univ.JM(boot.data, Flists[[2]], tol.rel = 1e-2)
    resid2.b <- resid(mod2.b, type = "pearson")[[1]]
    MSE2.b <- sum(resid2.b^2)/N
    ll2.b <- logLik(mod2.b)
    AICtabs.2[,b] <- c(c(ll2.b), attr(ll2.b, 'AIC'), attr(ll2.b, "BIC"))
    rm(mod2.b)
    MSEs[,b] <- c(MSE1.b, MSE2.b)
    utils::setTxtProgressBar(pb, b)
  }
  close(pb)
  
  out <- list(
    resp = resp, B = B, mod1 = Flists[[1]], mod2 = Flists[[2]],
    MSE1 = MSE1, MSE2 = MSE2, AICtab1 = AICtab1, AICtab2 = AICtab2,
    MSEs = MSEs, AICtabs1 = AICtabs.1, AICtabs2 = AICtabs.2
  )
  class(out) <- "zxc"
  return(out)
}

bb <- boot.MSE(resp = "ADAS13", competing = list(
  one = list(family = "gaussianraw"),
  two = list(family = "poisson")
), B = 50L)

plot.zxc <- function(x, ...){
  mod1 <- x$AICtabs1; mod2 <- x$AICtabs2
  MSEs <- x$MSEs
  par(mfrow=c(1,2))
  plot(MSEs[2,] ~ MSEs[1,], pch = 20,
       col = ifelse(MSEs[2,] < MSEs[1,], "lightgreen", "lightblue"),
       ylim = c(min(MSEs), max(MSEs)), xlim = c(min(MSEs), max(MSEs)),
       main = expression(bar("MSE")))
  abline(0, 1)
  plot(mod2[1,] ~ mod1[1,], pch = 20, main = "logLik",
       col = ifelse(mod2[1,] > mod1[1,], "lightgreen", "lightblue"),
       ylim = c(min(c(mod1[1,], mod2[1,])), max(c(mod1[1,], mod2[1,]))), 
       xlim = c(min(c(mod1[1,], mod2[1,])), max(c(mod1[1,], mod2[1,]))))
  abline(0, 1)
  par(mfrow=c(1,1))
}
