# "Noddy model"
# Just a univariate with SerBilir with splines
# and hist2 + sex in the survival sub-model
rm(list=ls())
source(".Rprofile")
library(splines)
long.formula <- list(
  serBilir ~ ns(time, knots = c(1,4)) + histologic2 + sex + (1 + ns(time, knots = c(1,4))|id)
)
surv.formula <- Surv(survtime, status) ~ age + histologic2

M <- joint(
  long.formulas = long.formula, 
  surv.formula = surv.formula, 
  data = PBCredx, family = list("gaussian")
  # control = list(tol.rel = 5e-3)
)
summary(M)


# Unpacking for benchmarking ----------------------------------------------
# dmats, coeffs and modelinfo -->
dmats <- M$dmats
sv <- dmats$surv
surv <- dmats$ph
Omega <- M$coeffs
Info <- M$ModelInfo
# REs -->
b <- M$REs
S <- attr(b, 'vcov')
Sigma <- lapply(1:Info$n, function(i) gmvjoint:::vech2mat(S[i,], Info$Pcounts$q))
b.hat <- lapply(1:Info$n, function(i) b[i,])
# GH -->
GH <- statmod::gauss.quad.prob(Info$control$gh.nodes, "normal")
w <- GH$w; v <- GH$n

# Source methods and benchmark
source("methods.R")

# This is the example of a "usual" fit
usual.fit.bench <- microbenchmark(
  `C++ entirely` = { gmvj() },
  `R wrapping C++` = { R.on.Cpp() },
  `R mimicking C++` = { R.way1() },
  times = 250L
)

df <- unpack.microbenchmark(usual.fit.bench)
library(ggplot2)
P <- microbenchmark:::autoplot.microbenchmark(usual.fit.bench)
P + theme_csda() + labs(x = "Elapsed time (milliseconds)")
P$theme[1]
ggsave("./CppR.png", width = 140, height = 60, units = "mm")
               