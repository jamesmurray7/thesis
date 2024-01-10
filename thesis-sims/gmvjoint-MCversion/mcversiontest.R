rm(list=ls())
source('.Rprofile')
sourceCpp("src/testMC.cpp")
# Trivariate fit as per Sec 4.5.3
# beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
beta <- rbind(
  c(-2, 0.1, -0.1, 0.2),
  c(2, -0.1, 0.1, 0.2),
  c(1, -1, 1, -1)
)
gamma <- c(0.5, -0.5, .5)
D <- diag(c(0.25, 0.09, 0.50, 0.09, 2))
D[1,3] <- D[3,1] <- D[1,5] <- D[5,1] <- D[3,5] <- D[5,3] <- .125
family <- list('gaussian', 'poisson', "binomial")
data <- gmvjoint::simData(ntms = 10, beta = beta, D = D, n = 250,
                family = family, zeta = c(0, -0.2),
                random.formulas = list(~time,~time,~1),
                sigma = list(0.16, 0, 0), gamma = gamma, theta = c(-3,.1))$data

# Specify formulae and target families
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id),  # Poisson
  Y.3 ~ time + cont + bin + (1|id)          # Binomial
)
surv.formula <- Surv(survtime, status) ~ bin

control <- list()
disp.formulas = NULL
MCtype <- 'ordinary'; N <- 1e2

a1 <- joint(long.formulas, surv.formula, data, family, disp.formulas = NULL,
           "ordinary", N = 1e3)
a2 <- joint(long.formulas, surv.formula, data, family, disp.formulas = NULL,
            "antithetic", N = 1e3)
a3 <- joint(long.formulas, surv.formula, data, family, disp.formulas = NULL,
            "sobol", N = 1e3)

# This is Sec 4.5.3 simulation
load("/data/c0061461/THESIS/GMVJM/data/Triv.RData")

# Just take K * 100 as in joineRML...
aa1 <- joint(long.formulas, surv.formula, data= sim.sets[[1]], family, MCtype = 'sobol', N = 300,
            control=list(verbose=T))
aa2 <- joint(long.formulas, surv.formula, data= sim.sets[[1]], family, MCtype = 'ordinary', N = 300,
            control=list(verbose=T))
aa3 <- joint(long.formulas, surv.formula, data= sim.sets[[1]], family, MCtype = 'antithetic', N = 300,
            control=list(verbose=T))
