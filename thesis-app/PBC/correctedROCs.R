rm(list=ls())
source(".Rprofile")
library(ggplot2)
library(splines)
load(save.dir.file.path(x = "Joint/Multivs/biv2.RData"))
PBCredx$histologic2 <- ifelse(PBCredx$histologic%in%c("3","4"), 1, 0)
source("./correctedROCs-helpers.R")

# Window 1: [2, 3.5]
W1 <- corrected.ROC(joint.biv2, data = PBCredx, Tstart = 2, delta = 1.5)
# Window 2: [3.5, 7]
W2 <- corrected.ROC(joint.biv2, data = PBCredx, Tstart = 3.5, delta = 3.5)
# Window 3: [7, 14]
W3 <- corrected.ROC(joint.biv2, data = PBCredx, Tstart = 7, delta = 7)

Ws <- list(W1 = W1, W2 = W2, W3 = W3)
save(Ws, file = save.dir.file.path("ROCcorrected.RData"))


# Comparing with univ model (M0) ------------------------------------------
serBilir.only <- joint(
  long.formulas = list(serBilir ~  ns(time, knots = c(1, 4)) + histologic2 + 
                         sex + (1 + ns(time, knots = c(1, 4)) | id)
  ),
  surv.formula = Surv(survtime, status) ~ age + histologic2,
  data = PBCredx, family = list("gaussian"), control = list(tol.rel = 5e-3)
)
W1.0 <- corrected.ROC(serBilir.only, data = PBCredx, Tstart = 2, delta = 1.5)
W2.0 <- corrected.ROC(serBilir.only, data = PBCredx, Tstart = 3.5, delta = 3.5)
W3.0 <- corrected.ROC(serBilir.only, data = PBCredx, Tstart = 7, delta = 7)
W0s <- list(W1 = W1, W2 = W2, W3 = W3)
save(W0s, file = save.dir.file.path("ROCcorrected_M0.RData"))

# 50 replicates (only) for the extra variance 0<
W1.01 <- two.corrected.ROCs(serBilir.only, joint.biv2, data = PBCredx, Tstart = 2, delta = 1.5)
W2.01 <- two.corrected.ROCs(serBilir.only, joint.biv2, data = PBCredx, Tstart = 3.5, delta = 3.5)
W3.01 <- two.corrected.ROCs(serBilir.only, joint.biv2, data = PBCredx, Tstart = 7, delta = 7)
W.01s <- list(W1 = W1.01, W2 = W2.01, W3 = W3.01)
save(W.01s, file = save.dir.file.path("TwoROCscorrected.RData"))

# Plotting W.01s
library(ggplot2)
library(dplyr)
source("../../theme_csda.R")
YY <- lapply(W.01s, function(X){
  stopifnot(inherits(X, "two.corrected.joints"))
  # Unpack things for plotting
  M0.point.AUC <- X$M0.auc; M0.point.PE <- X$M0.PE
  M1.point.AUC <- X$M1.auc; M1.point.PE <- X$M1.PE
  AUC.c.0 <- X$M0.AUC.c; AUC.c.1 <- X$M1.AUC.c
  PE.c.0 <- X$M0.PE.c; PE.c.1 <- X$M1.PE.c
  AUC.df <- data.frame(model = c(rep("M0", X$valid.boots), rep("M1", X$valid.boots)), 
                       what = "AUC", value = c(AUC.c.0, AUC.c.1))
  AUC.df <- AUC.df[AUC.df$value < 1, ]
  PE.df <- data.frame(model = c(rep("M0", X$valid.boots), rep("M1", X$valid.boots)), 
                      what = "PE", value = c(PE.c.0, PE.c.1))
  R <- data.frame(model = "M01", what = "R", value = X$R)
  df <- rbind(AUC.df, PE.df, R)
  # Work out what window we're in
  Ts <- X$Tstart
  num <- switch(as.character(Ts),
                "2" = 1,
                "3.5" = 2,
                "7" = 3)
  df$window <- paste0("w[", num, "]")
  # Strip titles
  df %>% 
    mutate(
      what2 = case_when(
        what == "R" ~ "R[c]",
        what == "PE" ~ "widehat(PE)[c]",
        T ~ "AUC[c]"
      )
    ) %>% 
    select(-what) %>% 
    rename(what = what2)
})
df <- do.call(rbind, YY)
df$what <- factor(df$what, c("AUC[c]", "widehat(PE)[c]", "R[c]"))
df$model <- factor(df$model, c("M0", "M1", "M01"))
df$window <- factor(df$window, c("w[1]", "w[2]", "w[3]"))

df %>% 
  ggplot(aes(x = forcats::fct_rev(window), y = value, fill = model)) + 
  geom_boxplot(lwd = .15, # colour = 'grey80',
               outlier.alpha = .33, outlier.size = .5,
               fatten = 2) + 
  facet_wrap(~what, scales = 'free_x', labeller = label_parsed) + 
  scale_x_discrete("", labels = function(l) parse(text = l)) + 
  scale_y_continuous("", breaks = scales::pretty_breaks(5)) + 
  scale_fill_manual(values = c(.nice.orange, "#F03B20", "#4C98FE")) +
  theme_csda() + 
  theme(
    legend.position = 'none',
    axis.ticks.length.y = unit(0, 'mm'),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  coord_flip() 
  
ggsave("./CompareROCs.png", width = 140, height = 70, units = 'mm')

