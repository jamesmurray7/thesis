#' #################################
#' Shape of integrands             #
#' -------------------             #
#' 2. Exponent in e.g. poisson     #
#'    linear predictor             #
#' #################################

rm(list=ls())
source(".Rprofile")
# Generate Gaussian data
dd <- simData(n = 500, zeta = c(0,-0.2),
              beta = c(2.0,-0.1,0.1,0.2),
              D = matrix(c(.5, .125, .125, .09), 2, 2),
              gamma = .5, family = list("negbin"))$data

jj <- joint(
  list(
    Y.1 ~ time + cont + bin + (1 + time|id)
  ),
  Surv(survtime, status) ~ bin,
  dd, list("negbin")
)

GH <- statmod::gauss.quad.prob(3, 'normal')
v <- GH$nodes; w <- GH$weights

# Take id = 1 -- this is purely for illustrative purposes
L <- jj$dmats$long
S <- jj$dmats$surv
M <- jj$ModelInfo
Samp <- Metropolis_Hastings(
  b = jj$REs[1,], Omega = jj$coeffs, Y = L$Y[[1]], X = L$X[[1]], Z = L$Z[[1]],
  W = L$W[[1]], family = jj$ModelInfo$family, Delta = jj$dmats$ph$Delta[[1]],
  S = S$S[[1]], Fi = S$Fi[[1]], l0i = S$l0i[[1]], SS = S$SS[[1]], Fu = S$Fu[[1]],
  haz = S$l0u[[1]], gamma_rep = rep(jj$coeffs$gamma, sapply(M$inds$R$b, length)), 
  beta_inds = M$inds$Cpp$beta, b_inds = M$inds$Cpp$b, K = M$K, q = M$Pcounts$q,
  burnin = 500L, N = 10000L, Sigma = gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[1,], M$Pcounts$q),
  b_dist = "t", df = 4L, tune = 3.5
)
W <- t(Samp$walks)

# Shape of E[\exp{SS^T\zeta + \sum_k\gamma_k\b{_k)}]
b <- jj$REs[1,]
Sig <- gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[1,], M$Pcounts$q)
tau <- sqrt(diag(tcrossprod(L$Z[[1]][[1]] %*% Sig, L$Z[[1]][[1]])))
I <- apply(W,1,function(x){
  exp(L$X[[1]][[1]]%*%jj$coeffs$beta + L$Z[[1]][[1]] %*% x)
})


library(ggplot2)
dens <- apply(I, 1, density)
hat <- exp(L$X[[1]][[1]]%*%jj$coeffs$beta + L$Z[[1]][[1]] %*% b)
med <- hat
mean <- hat * exp(tau^2/2)
plot.dens.with.point <- function(i){
  di <- dens[[i]]
  x <- di$x; y <- di$y
  plot(x, y ,'l', 
       xlab = expression(exp(X[i]*beta*"+"*Z[i]*hat(b)[ik])),
       ylab = 'Density')
  abline(v = med[i], col = 'blue', lty = 'dotted')
  abline(v = mean[i], col = 'red', lty = 'dotted')
  #      ylab = "Density",
  #      xlim = c(xlim.lb[i], xlim.ub[i]))
  # points(nodes[1,], c(0,0,0), pch = 'x', col = 'red')
}

plot.dens.with.point(5)

nr <- nrow(I)
ceiling(nr/2)
png(save.dir.file.path("NegbinLinpExp.png"),
    width = 140, height = 160, units = 'mm', res = 5e2)
par(mfrow = c(3,1))
for(pp in c(1, ceiling(nr/2), nr)) # start, midway, end
  plot.dens.with.point(pp)
dev.off()

