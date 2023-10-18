rm(list=ls())
source(".Rprofile")
library(ggplot2)
load(save.dir.file.path(x = "Joint/Multivs/biv2.RData"))

# Longit. residuals vs. fitted --------------------------------------------
R.long <- residuals(joint.biv2,'l','p')

alpha.names <- joint.biv2$ModelInfo$Resps
out <- setNames(vector("list", length(alpha.names)), alpha.names)
for(f in alpha.names){
  resp <- ifelse(f == "serBilir", "Serum bilirubin", tools::toTitleCase(f))
  
  # Make ggplots for better automatic scaling...
  R <- R.long[[f]]
  x <- attr(R, "fitted")
  out[[f]] <- data.frame(response = resp, fitted = x, pr = R, row.names = NULL)
}
out.all <- do.call(rbind, out)
out.all$response <- factor(out.all$response, c("Serum bilirubin", "Albumin"))

ggplot(out.all, aes(x = fitted, y = pr)) + 
  geom_hline(aes(yintercept = 0), lty = 5, col = 'red', lwd = 0.25) + 
  geom_point(size=.005, alpha = .2) + 
  facet_wrap(~ response, nrow = 1, scales = 'free') + 
  labs(y = "Pearson residuals", x = "Fitted") + 
  theme_csda()+ 
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    axis.title = element_text(size=5),
    strip.text = element_text(size=7,vjust=1)
  )

ggsave("FinalPearsonResiduals.png", width = 140, height = 60, units = "mm")


# Cox-Snell residuals -----------------------------------------------------
R.surv <- residuals(joint.biv2, 's')

png("FinalCoxSnellExamples.png", width = 140, height = 120, units = "mm", res = 5e2)
sf <- survfit(Surv(R.surv, unlist(joint.biv2$dmats$ph$Delta)) ~ 1)
plot(sf, mark.time = F, lwd = 0.75,
     xlab = "Cox-Snell residuals", ylab = "Survival probability")
curve(exp(-x), from = 0, to = max(sf$time), add = T, lwd = 1, col = 'steelblue')
dev.off()


# Conditional distribution of random effects ------------------------------
CR <- cond.ranefs(joint.biv2, tune = 1.1, burnin = 2000, N = 10000)
png("FinalCR.png", width = 140, height = 140, units = "mm", res = 5e2)
plot(CR)
dev.off()


# LRT with joint.biv ------------------------------------------------------
load(save.dir.file.path(x = "Joint/Multivs/biv.RData"))

anova(joint.biv2, joint.biv) # No evidence at 5% level to include sex in model.
                             # But there is some at 10% level.

# ROC
aa <- ROC(joint.biv2, PBCredx, 7, 2, control = list(b.densit = 't', df = 4, nsim=0))
