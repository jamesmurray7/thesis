library(glmmTMB)

data(PBC, package = "gmvjoint")
data <- na.omit(PBC[, c("id", "survtime", "status", 
                        "albumin", "platelets", "time", "drug")])

mod1 <- glmmTMB(albumin ~ drug * time + (1 + time|id), data  = data)
mod2 <- glmmTMB(platelets ~ drug * time + (1 + time|id), data = data, family = poisson)
mod3 <- glmmTMB(platelets ~ drug * time + (1 + time|id), data = data, family = glmmTMB::genpois())

# Collate (standardised) residuals
p1 <- unname(resid(mod1, type = "pearson"))
p2 <- unname(resid(mod2, type = "pearson"))
p3 <- unname(resid(mod3, type = "pearson"))
resids <- as.data.frame(
  cbind(`Albumin (Gaussian)` = p1, `Platelet count (Poisson)` = p2,
        `Platelet count (Generalised Poisson)` = p3)
  )
resids$id <- as.numeric(as.character(data$id))
resids2 <- tidyr::pivot_longer(
  resids, -id,
  values_to = "resid"
)

# Collate fitted values
f1 <- unname(fitted(mod1, type = "response"))
f2 <- unname(fitted(mod2, type = "response"))
f3 <- unname(fitted(mod3, type = "response"))

fitteds <- as.data.frame(
  cbind(`Albumin (Gaussian)` = f1, `Platelet count (Poisson)` = f2,
        `Platelet count (Generalised Poisson)` = f3)
)
fitteds$id <- as.numeric(as.character(data$id))
fitteds2 <- tidyr::pivot_longer(
  fitteds, -id,
  values_to = "fitted"
)

df <- dplyr::left_join(resids2, fitteds2, by = c("id", "name"))

library(ggplot2)
source("../thesis-sims/theme_csda.R")
ggplot(df, aes(x = fitted, y = resid)) + 
  geom_point(size = 0.05, alpha = 0.05) + 
  facet_wrap(~name, scales = "free") + 
  geom_smooth(method = "loess", formula = y~x, se = F, colour = .nice.orange) + 
  theme_csda() + 
  scale_x_continuous(breaks = scales::pretty_breaks(6)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  labs(x = expression(hat(Y)[ij]),
       y = expression(hat(r)[ij]^{"(P)"})) + 
  theme(
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 5.5),
    axis.text = element_text(size = 4)
  )

ggsave("./output/exampleResids.png", width = 140, height = 60, units = "mm", dpi = 500)
       