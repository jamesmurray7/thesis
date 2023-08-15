t1 <- glmmTMB(ADAS13 ~ APOE4 * (time) + (1 + time|id),
              data = adni, family = poisson)
t2 <- glmmTMB(ADAS13 ~ APOE4 * (time) + (1 + time|id),
              data = adni, family = genpois())

vuong.test <- function(poiss.fit, gp.fit){
  
  # Parameters
  N <- nobs(poiss.fit)
  beta1 <- fixef(poiss.fit)$cond
  K1 <- length(beta1)
  beta2 <- fixef(gp.fit)$cond
  sigma2 <- fixef(gp.fit)$disp
  K2 <- length(beta2) + length(sigma2)
  
  # Predicted values
  mu1 <- predict(poiss.fit, type = "conditional")
  mu2 <- predict(gp.fit, type = "conditional")
  disp2 <- c(getME(gp.fit, "Xd") %*% sigma2)
  Y <- poiss.fit$frame[, glmmTMB:::responseName.default(poiss.fit)]
  
  # Likelihoods
  L1 <- dpois(Y, mu1); L2 <- HMMpa::dgenpois(Y, mu2, disp2)
  
  # All bits we need --> Test statistic?
  log.rat <- log(L1/L2)
  sum.log.rat <- sum(log.rat)
  sum.log2.rat <- sum(log.rat^2)/N
  sum.log.rat2 <- (sum(log.rat)/N)^2
  omega2 <- sum.log2.rat - sum.log.rat2
  Tstat <- sum.log.rat/(omega2 * sqrt(N)) # The test statistic
  
  
}