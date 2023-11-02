data.dir <- save.dir.file.path("data")
out.dir <- save.dir.file.path("fits")
# Negative binomial, fixed effects appear to suffer a little bit in original results
fn.nb <- create.simfn(family = list("negbin"),
                      arguments = list(D = matrix(c(0.50,0,0.,0.09),2,2),
                                       beta = c(0.5, 0.05, 0.1, -0.1), n = 500L,
                                       theta = c(-2.9, 0.1))) # about 30% failure rate
data <- fn.nb()

test <- joint(
  list(Y.1 ~ time + cont + bin + (1 + time|id)),
  surv.formula = Surv(survtime, status) ~ bin,
  data, family = list("negbin"),
  control = list(return.dmats = FALSE, tol.rel = 5e-3)
)

# Ok create N = 300
sim.sets <- createNsims(fn.nb, 300L)
DataSummary(sim.sets) # 28.37% failure rate average
out.file <- save.dir.file.path("UnivNBredux.RData", data.dir)
save(sim.sets, file =  out.file)

# Function to fit one dataset
fit.one <- function(X){
  tryCatch(joint(
    long.formulas = list(Y.1 ~ time + cont + bin + (1 + time|id)),
    surv.formula = Surv(survtime, status) ~ bin,
    data = X, family = list("negbin"),
    control = list(return.dmats = FALSE, tol.rel = 5e-3) # avoid v large lists.
  ), error = function(e) NULL)
}

# Fitting all w/ progress bar
fits <- vector("list", 300L)
pb <- utils::txtProgressBar(max = 300L, style = 3)
for(j in 1:300){
  fits[[j]] <- fit.one(sim.sets[[j]])
  utils::setTxtProgressBar(pb, j)
}
close(pb)
save(fits, file = save.dir.file.path("UnivNBredux.RData", out.dir))

# Extracting bits from ParseUnivSims.R
this.parsed <- lapply(fits, extract.from.joint)
Omega <- sapply(this.parsed, function(y) y$Omega)
SE <- sapply(this.parsed, function(y) y$SE)
lb <- sapply(this.parsed, function(y) y$lb)
ub <- sapply(this.parsed, function(y) y$ub)
et <- unname(sapply(this.parsed, function(y) y$elapsed))
tt <- unname(sapply(this.parsed, function(y) y$total))
mo <- unname(sapply(this.parsed, function(y) y$model))
out <- list(Omega = Omega, SE = SE, lb = lb, ub = ub, elapsed = et, total = tt,
            model = unique(mo))

target.mat <- create.targetmat(family = list("negbin"),
                               arguments = list(D = matrix(c(0.50,0,0.,0.09),2,2),
                                                beta = t(c(0.5, 0.05, 0.1, -0.1)), n = 500L,
                                                theta = c(-2.9, 0.1)),
                               N = 300L)
# ok cant be bothered to change underlying function, it's ignoring \beta though for some reason.
target.mat[,3] <- .5  # intercept
target.mat[,4] <- .05 # time
target.mat[,5] <- .1  # cont
target.mat[,6] <- -.1 # bin
target.mat <- t(target.mat)

# Tabulation
rn <- row.names(Omega); rn2 <- row.names(target.mat)

Omega <- Omega[-2,]; SE <- SE[-2,]; lb <- lb[-2,]; ub <- ub[-2,] # Remove offdiags
row.names(Omega) <- row.names(SE) <- row.names(lb) <- row.names(ub) <-  rn2

# Mean (SD) / (average SE)
(Emp.Mean <- rowMeans(Omega))
(Emp.SD <- apply(Omega, 1, sd))
(Avg.SE <- rowMeans(SE))

# CP
(CP <- rowSums(lb <= target.mat & ub >= target.mat) / 300)

# MSE and Bias
MSE <- rowMeans((target.mat - Omega)^2)
Bias <- rowMeans(Omega-target.mat)

x <- list(Emp.Mean = Emp.Mean, Emp.SD = Emp.SD,
          Avg.SE = Avg.SE, CP = CP, MSE = MSE, Bias = Bias,
          targets = target.mat[,1], nm ="negative binomial", mod = "negbin")


# xtable-ing
.f <- function(x) format(round(x, 3), nsmall = 3, justify = 'right')
this <- x
param <- paste0("$", names(this$Emp.Mean), ' = ', trimws(.f(this$targets)), "$")
param <- gsub("\\]", "}", gsub("D\\[","D_{", param))

MSD <- paste0(.f(this$Emp.Mean), ' (', .f(this$Emp.SD), ')')
SE <- .f(this$Avg.SE)
Bias <- .f(this$Bias)
MSE <- .f(this$MSE)
CP <- format(round(this$CP, 2), nsmall = 2)

out <- data.frame(Parameter = param, "MSD" = MSD, `SE` = SE, Bias = Bias, MSE = MSE, CP = CP,
                  row.names = NULL)
names(out)[2] <- "Mean (SD)"
# if(i > 1) out$Parameter <- NULL
tab <- out

# x-tableau
library(xtable)
align.lhs <- "cl|"
align.rhs <- paste0(rep("r", ncol(tab)-1), collapse = '')
align <- paste0(align.lhs, align.rhs)
xt <- xtable(tab, align = align)

print(xt, sanitize.text.function = identity,booktabs = FALSE,
      include.rownames = FALSE, size = 'scriptsize')

