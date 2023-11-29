# Chen (2021) implementation (I think)

Chen <- function(power, 
                 Tmed, n, tbar, D, fail.rate,
                 alpha = 0.95, gamma = 0.5){
  zbeta.tilde <- qnorm(1-power)
  zalpha.tilde <- qnorm(1-alpha)
  h <- log(2)/Tmed
  h2 <- h^2
  eht <- exp(-h * tbar)
  D11 <- D[1,1]; D22 <- D[2,2]; D21 <- D[2,1]
  H <- fail.rate
  G <- D11 + D22/H * (2/h2 - eht * (tbar^2 + 2 * tbar/h + 2/h2)) + 
          D21/H * (1/h - eht * (tbar + 1/h))
  Xi <- fail.rate * n
  abs(Xi - (zalpha.tilde + zbeta.tilde)^2/(gamma^2 * G))
}

# We want to minimise the difference
Chen(0.8, 2.5, 250, 3, diag(c(1,1)), 0.1)

# Univariate data
ss <- simData(n = 250, D = matrix(c(0.25,0.125,0.125,0.09), 2, 2), beta = c(1, 0.10, 0.33, -0.50),
              gamma = 0.5, family = list("gaussian"), sigma = list(.16))
# (Exponential) rate?
SS <- survreg(Surv(survtime, status) ~ 1, data = ss$surv.data, dist = 'exponential')
rr <- 1/SS$coefficients
# Median survival time (assuming this is _either failure or censor) ---
Tmed <- median(ss$surv.data[!duplicated(ss$surv.data$id), 'survtime'])
# average follow-up time ---
tbar <- mean(ss$data$time)

optim(0.6, Chen, NULL, 
      Tmed = 2.5, n = 250, tbar = 2, D = matrix(c(0.25,0.125,0.125,0.09), 2, 2), fail.rate = rr,#1/2.578290,
      method = 'Brent', lower = 1e-5, upper = .999, control = list(fnscale = 1))

optim(0.6, Chen, NULL, 
      Tmed = Tmed, n = nrow(ss$surv.data), tbar = tbar, D = matrix(c(0.25,0.125,0.125,0.09), 2, 2), fail.rate = rr,#1/2.578290,
      method = 'Brent', lower = 1e-5, upper = .999, control = list(fnscale = 1))

# Make a wrapper for this -->
Chenwrapper.univ <- function(n = 250, D =  matrix(c(0.25,0.125,0.125,0.09), 2, 2), gamma = 0.5){
  ss <- simData(n = n, D = D, beta = c(1, 0.10, 0.33, -0.50), gamma = gamma, 
                family = list("gaussian"), sigma = list(.16))
  # (Exponential) rate?
  SS <- survreg(Surv(survtime, status) ~ 1, data = ss$surv.data, dist = 'exponential')
  rr <- 1/SS$coefficients
  # Median survival time (assuming this is _either failure or censor) ---
  Tmed <- median(ss$surv.data[!duplicated(ss$surv.data$id), 'survtime'])
  # average follow-up time ---
  tbar <- mean(ss$data$time)
  optim(0.6, Chen, NULL, 
        Tmed = Tmed, n = n, tbar = tbar, D = D, fail.rate = rr, gamma = gamma,
        method = 'Brent', lower = 1e-5, upper = .999, control = list(fnscale = 1))
}

Chenwrapper.univ(n = 5000, gamma = 0.1)
