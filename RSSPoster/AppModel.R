# The models
library(gmvjoint)
library(splines)
PBC <- na.omit(PBC[, c("id", "drug", "survtime", "status",
                       "time", "serBilir", "platelets", "alkaline")])
long.formulas <- list(
  log(serBilir) ~ drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4)) | id),
  platelets ~ drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4)) | id),
  alkaline ~ drug * time + (1 + time|id)
)
disp.formulas <- list(
  ~1,~1,~time
)
surv.formula <- Surv(survtime, status) ~ drug
family <- list("gaussian", "poisson", "negbin")

# Trivariate fit: About 2 mins
mod <- joint(long.formulas, surv.formula, PBC, family, disp.formulas, control = list(tol.rel=5e-3))
# Univariates for comparison? -- Pretty quick
mod1 <- joint(list(log(serBilir) ~ drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4)) | id)),
              surv.formula, PBC, list("gaussian"), control = list(tol.rel=5e-3)) 
mod2 <- joint(list(platelets ~ drug * ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4)) | id)),
              surv.formula, PBC, list("poisson"), control = list(tol.rel=5e-3)) 
mod3 <- joint(list(alkaline ~ drug * time + (1 + time|id)),
              surv.formula, PBC, list("negbin"), control = list(tol.rel=5e-3)) 

# Plot 1: gamma, zeta -----------------------------------------------------
qz <- qnorm(.975)
# Parse trivariate
triv <- data.frame(
  model = "Trivariate fit",
  parameter = c("gamma[bilirubin]", "gamma[platelets]", "gamma[alkaline]", "zeta[drug]"),
  estimate = c(mod$coeffs$gamma, mod$coeffs$zeta),
  SE = mod$SE[(length(mod$SE)-3):length(mod$SE)],
  row.names = NULL
)
triv$lower <- triv$estimate - qz * triv$SE
triv$upper <- triv$estimate + qz * triv$SE

univs <- do.call(rbind, lapply(list(mod1, mod2, mod3), function(x){
  resp <- x$ModelInfo$Resps
  insert <- ifelse(grepl("serBilir", resp), "bilirubin",
                   ifelse(grepl("alkaline", resp), "alkaline", "platelets"))
  S <- summary(x)$Su
  S <- as.data.frame(rbind(S[2,],S[1,]))[,c(1,2,5,6)]
  S <- cbind(model = "Univariate fit", cbind(parameter = c(paste0("gamma[", insert, "]"),
                                                       paste0("zeta[drug,", insert, "]")), S))
  row.names(S) <- NULL
  names(S) <- names(triv)
  S
}))

all.mods <- rbind(triv,univs)

# Plotting
library(ggplot2)
library(dplyr)
all.mods %>%
  dplyr::filter(parameter %in% c("gamma[bilirubin]", "gamma[platelets]", "gamma[alkaline]", "zeta[drug]")) %>% 
  mutate(parameter = factor(parameter, c("gamma[bilirubin]", "gamma[platelets]", "gamma[alkaline]", "zeta[drug]"),
                            c("gamma[bilirubin]", "gamma[platelets]", "gamma[alkaline]", "zeta[drug]"))) %>% 
  ggplot(aes(x = forcats::fct_rev(parameter), y = estimate, colour = model)) + 
  geom_hline(aes(yintercept=0), lty = 5, col = 'grey') + 
  coord_flip() +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(.25), width = .25, show.legend=T) +
  geom_point(pch = 19, position = position_dodge(.25),size = 2.1, show.legend = T) + 
  scale_color_manual(values=c(nclred, nclblue)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(6), limits = c(min(all.mods$lower), max(all.mods$upper))) +
  scale_x_discrete(labels = function(l) parse(text=l)) + 
  labs(
    x = NULL, y = "Parameter Estimate", colour = NULL
  ) + 
  theme_csda() + 

  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 21, colour = nclblue),
    legend.text = element_text(angle = 0, size = 15),
    legend.box.spacing = unit(-2, "mm"),
    legend.key.size = unit(6, "mm"),
    legend.position = c(0.845, 0.15),
    legend.key = element_blank()
  )

png("./PBCmodel.png", width = 0.45*linewidth, height = 90, units = "mm",
    res = 1500)
last_plot()
dev.off()
