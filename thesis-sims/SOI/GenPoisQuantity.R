#' #################################
#' Shape of integrands             #
#' -------------------             #
#' 2. Quantity exclusive to GP-1   #
#'    update                       #
#' #################################

rm(list=ls())
source(".Rprofile")
# Generate GP-1 data
dd <- simData(n = 500, zeta = c(0,-0.2), sigma = list(0.2),
              beta = c(0.5, -0.2, 0.05, 0.4),
              gamma = 0.5, family = list("genpois"),
              D = matrix(c(0.5, .125, .125, .05), 2, 2),theta = c(-2.95, .1))$data

jj <- joint(
  list(
    Y.1 ~ time + cont + bin + (1 + time|id)
  ),
  Surv(survtime, status) ~ bin,
  dd, list("genpois")
)

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
  burnin = 1000L, N = 50000L, Sigma = gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[1,], M$Pcounts$q),
  b_dist = "t", df = 4L, tune = 3.5
)
Samp$AcceptanceRate
W <- t(Samp$walks)

# Shape of log(1 + exp())
b <- jj$REs[1,]
Sig <- gmvjoint:::vech2mat(attr(jj$REs, 'vcov')[1,], M$Pcounts$q)
tau <- sqrt(diag(tcrossprod(L$Z[[1]][[1]] %*% Sig, L$Z[[1]][[1]])))
I <- apply(W,1,function(x){
  log(exp(L$X[[1]][[1]]%*%jj$coeffs$beta + L$Z[[1]][[1]] %*% x) + L$Y[[1]][[1]] * (L$W[[1]][[1]]%*%jj$coeffs$sigma[[1]]))
})


dens <- apply(I, 1, density)
plot.dens.with.point <- function(i){
  di <- dens[[i]]
  x <- di$x; y <- di$y
  plot(x, y ,'l', 
       xlab = expression(log(Y[i]*phi[i]*"+"*exp(X[i]*beta*"+"*Z[i]*b[i]))),
       ylab = 'Density')
}

plot.dens.with.point(5)

nr <- nrow(I)
ceiling(nr/2)
png(save.dir.file.path("GeneralisedPoissonLogBit.png"),
    width = 140, height = 160, units = 'mm', res = 5e2)
par(mfrow = c(3,1))
for(pp in c(1, ceiling(nr/2), nr)) # start, midway, end
  plot.dens.with.point(pp)
dev.off()


