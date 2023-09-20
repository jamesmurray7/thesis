# Example code: ROC curve for _binary_ outcome
library(survival)
data <- survival::pbc
data <- data[!is.na(data$trt),]
data$status2 <- ifelse(data$status==2, 1, 0)
ph <- coxph(Surv(time, status2) ~ trt, data, x = T)

# Obtain survival probabilities
S <- survival:::predict.coxph(ph, data, "survival")
pi <- 1 - S
C <- seq(0, 1, length.out = 101)
M <- structure(outer(pi, C, "<="),
               dimnames = list(paste0("id: ", 1:length(S)),
                               paste0("c_j = ", C)))
incident <- data$status2 == 1L
no.incident <- !incident
# True positives
TP <- colSums(M * c(incident))
# False negative
FN <- sum(incident) - TP
# False positives
FP <- colSums(M * c(no.incident))
# True negatives
TN <- sum(no.incident) - FP
# True positive rate (Sens.)
TPR <- TP/(TP + FN) 
# False positive rate (1 - spec.)
FPR <- FP/(FP + TN) 
# Youden's J statistic
J <- TPR + TN/(TN + FP) - 1
# Accuracy
Acc <- (TP + TN) / (TP + TN + FP + FN)
# Pos. Pred. Value (precision)
PPV <- TP / (TP + FP)
# Neg. Pred. Value
NPV <- TN/(TN + FN)
# F1 score
F1 <- 2 * (PPV * TPR) / (PPV + TPR)

# Make metric df
# Making a nice dataframe to report
out <- data.frame(threshold = C,
                  TP = TP, TN = TN, FP = FP, FN = FN,
                  TPR = TPR, FPR = FPR, PPV = PPV, NPV = NPV,
                  Acc = Acc, F1 = F1, J = J)
row.names(out) <- NULL

# Flip table so if multiple thresholds have 
# same TPR/FPR then we take the largest threshold
out <- out[order(out$threshold, decreasing = T), ]
# Remove duplicated TPR/FPR
out <- out[!duplicated(out[, c('TPR', 'FPR')]),]

png("output/exampleROC.png", width = 140, height = 110, units = "mm",
    res = 1e3)
# Make plot
plot(out$FPR, out$TPR, xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "ROC curve for PBC data", type = "l")
abline(0,1)
# Draw-on Youden's J
Js <- out[which.max(out$J),]
arrows(x0 = Js$FPR, x1 = Js$FPR,
       y0 = Js$FPR, y1 = Js$TPR,
       length = 0, lty = 5)
# Calculate and print AUC
sens <- out$TPR; omspec <- out$FPR;
height <- 0.5 * (sens[-1] + sens[-length(sens)])
width <- -diff(omspec)
AUC <- sum(height * width)
legend("topleft", legend = sprintf("AUC: %.3f", AUC),
       bty = "n")

dev.off()
