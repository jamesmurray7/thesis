rm(list=ls())
# Example code: Cox-Snell residuals for dummy
#     Cox PH models
library(survival)
data <- survival::pbc
data <- data[!is.na(data$trt),]
data$status2 <- ifelse(data$status==2, 1, 0)
data$trt <- data$trt-1
# Fit (dummy) model
ff <- function(covar = c("sex", "trt")){
  form <- as.formula(paste0("Surv(time, status2) ~ ",
                            covar), env = parent.frame())
  ph <- coxph(form, data, x = T)
  
  # Extract necessary things 
  C <- coxph.detail(ph)
  zeta <- ph$coefficients
  S_i <- lapply(1:nrow(data), function(i) ph$x[i,,drop=F])
  
  # Estimate for baseline hazard
  l0 <- data.frame(
    time = C$time, 
    nevent = C$nevent, 
    hazard = C$hazard
  )
  
  # CS residuals
  xx <- sapply(unique(data$id), function(i){
    this <- data[data$id == i, ]
    Ti <- this$time
    l0i <- as.vector(l0[l0$time <= Ti, "hazard"])
    Si <- S_i[[i]]
    Si <- apply(Si, 2, rep, length(l0i))
    if(!"matrix"%in%class(Si)) Si <- t(Si)
    crossprod(l0i, exp(Si %*% zeta))
  })
  
  # Plot
  S <- survfit(Surv(xx, data$status2)~1)
  plot(S, mark.time = F,
       xlab = "Cox-Snell residuals",
       ylab = "Survival Probability",
       main = "Survival function of Cox-Snell residuals")
  # Overlaid unit exponential exp{-t}
  cc <- curve(exp(-x), from = 0, to = max(S$time), add = T,
        col = .nice.orange, lwd = 2)
  cc
  
  list(S=S,cs=xx,maxtime=max(S$time),covar=covar,
       cc = cc)
}

S.sex <- ff("sex"); S.trt <- ff("trt")

png(file = "./output/exampleCoxSnell.png",
    width = 140, height = 60, units = "mm", res = 1e3,
    pointsize=8)
par(mfrow=c(1,2))
plot(S.sex$S, mark.time = F,
     xlab = "Cox-Snell residuals",
     ylab = "Survival Probability",
     main = "Covariate: Sex")
curve(exp(-x), from = 0, to = S.sex$maxtmaxtime, add = T,
      col = .nice.orange, lwd = 2)
plot(S.trt$S, mark.time = F,
     xlab = "Cox-Snell residuals",
     ylab = "Survival Probability",
     main = "Covariate: Study drug")
curve(exp(-x), from = 0, to = S.trt$maxtmaxtime, add = T,
      col = .nice.orange, lwd = 2)
legend("topright", bty = "n", lty = 1, col = .nice.orange,
       legend = expression(exp(-x)))
par(mfrow=c(1,1))
dev.off()\\\\\\\\