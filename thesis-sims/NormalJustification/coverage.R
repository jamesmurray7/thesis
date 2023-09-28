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

prop.by.mi <- function(X, accept.min = 0.20, accept.max = 0.25){
  stopifnot(inherits(X, "Sample"))
  df <- do.call(rbind, lapply(1:NROW(X$Walks), function(i){
    this.df <- X$df[X$df$id==i,]
    this.b.hat <- c(b0 = this.df[this.df$var=="b[0]", "hatb"][1], b1 = this.df[this.df$var=="b[1]", "hatb"][2])
    psi <- check.walks(X$Walks[[i]], this.b.hat, X$Sigmas[[i]])
    data.frame(Acc = X$Acc[i], mi = X$mi[i], psi = psi)
  }))
  to.keep <- df$Acc >= accept.min & df$Acc <= accept.max
  df <- df[to.keep,-1]
  df[order(df$mi),]
}

a <- prop.by.mi(X)
plot(a[,2]~a[,1])



X_ <- create.simfn(arguments = list(n = 500, ntms = 20, theta = c(-2.5,0.1)))
X_ <- dataGen(X_)
X <- Sample(X_, 2.9, T)
pp <- prop.by.mi(X)
aa <- tapply(pp$psi, pp$mi, median)
plot(as.numeric(names(aa)), aa, "l", ylim = 0:1)
abline(h = 0.95,col='red',lty=3)
