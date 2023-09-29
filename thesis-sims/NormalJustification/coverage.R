# Obtaining coverage on a "by one" basis, and then constructing wrapper.

# Obtaining phi, the proportion "inside" the region bounded by ellipse.
check.walks <- function(W, b.hat, Sigma){
  eig <- eigen(Sigma)
  rx <- sqrt(eig$values[1] * qchisq(0.95, 2))
  ry <- sqrt(eig$values[2] * qchisq(0.95, 2))
  rx2 <- rx^2; ry2 <- ry^2
  theta <- atan2(eig$vec[2,1], eig$vec[1,1])
  ct <- cos(theta); st <- sin(theta)
  dx <- W[,1] - b.hat[1]
  dy <- W[,2] - b.hat[2]
  checks <- (ct * dx + st * dy)^2/rx2 + (st * dx - ct * dy)^2/ry2 <= 1
  # Return the proportion
  sum(checks)/length(checks)
}

# Wrap over an object of class `Sample`.
prop.by.mi <- function(X, accept.min = 0.20, accept.max = 0.25){
  stopifnot(inherits(X, "Sample"))
  df <- do.call(rbind, lapply(1:NROW(X$Walks), function(i){
    this.df <- X$df[X$df$id==i,]
    this.b.hat <- c(X$b.hats[[i]])
    psi <- check.walks(X$Walks[[i]], this.b.hat, X$Sigmas[[i]])
    data.frame(Acc = X$Acc[i], mi = X$mi[i], psi = psi)
  }))
  df$surv <- as.numeric(X$include.survival)
  to.keep <- df$Acc >= accept.min & df$Acc <= accept.max
  df <- df[to.keep,-1]
  df[order(df$mi),]
}

quick.psi.plot <- function(psi1, psi2){
  qz <- qnorm(.975)
  
  mean1 <- with(psi1, tapply(psi, mi, mean))
  mean2 <- with(psi2, tapply(psi, mi, mean))
  med1 <- with(psi1, tapply(psi, mi, median))
  med2 <- with(psi2, tapply(psi, mi, median))
  qn1 <- do.call(rbind, with(psi1, tapply(psi, mi, quantile, c(0.25, 0.75))))
  qn2 <- do.call(rbind, with(psi2, tapply(psi, mi, quantile, c(0.25, 0.75))))
  sd1 <- with(psi1, tapply(psi, mi, sd))
  sd2 <- with(psi2, tapply(psi, mi, sd))
  
  lb1 <- mean1 - qz * sd1; ub1 <- mean1 + qz * sd1
  lb2 <- mean2 - qz * sd2; ub2 <- mean2 + qz * sd2
  
  one <- data.frame(mi = as.numeric(names(mean1)), mean = mean1, sd = sd1, lb = lb1, 
                    ub = ub1, median = med1, qn1,  surv = unique(as.logical(psi1$surv)))
  two <- data.frame(mi = as.numeric(names(mean2)), mean = mean2, sd = sd2, lb = lb2, 
                    ub = ub2, median = med2, qn2 ,surv = unique(as.logical(psi2$surv)))
  
  df <- rbind(one,two)
  
  ggplot(df, aes(x = mi, y = median, group = surv, colour = surv)) +
    # geom_hline(yintercept = c(0.95, 0.975), lty = 3, colour = "red") + 
    geom_hline(yintercept = c(0.975), lty = 3, colour = "red") + 
    geom_errorbar(aes(ymin = `X25.`, ymax = `X75.`), width = 0.2, position = position_dodge(width = 0.15)) + 
    geom_point(position = position_dodge(width = 0.15)) + 
    scale_x_continuous(breaks=1:10) +
    labs(x = expression(m[i]), y = expression(psi[i])) + 
    expand_limits(y=c(1)) + 
    theme_csda()
  
    
  
}