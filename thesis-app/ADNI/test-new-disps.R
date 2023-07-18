library(splines)
dd <- cleandata("ADAS11")
a <- glmmTMB(
  formula = ADAS11 ~ age_scaled + gender + APOE4 * ns(time, knots = c(0.5, 1.5)) + (1 + ns(time, knots = c(0.5, 1.5))|id),
  data = dd, family = genpois(), dispformula = ~1
)

sigma <- fixef(a)$disp
mu <- predict(a, type =  "conditional")
# Var1: as per glmmTMB::family_glmmTMB; Var2 as per Z&I.
Var1 <- mu * exp(sigma)
Var2 <- (1 + sigma)^2 * mu
# Heavier underdispersion...
(sigma.new <- sqrt(Var2/mu)-1) 
(sigma.old.way <- exp(sigma/2)-1)

ii <- unname(which(table(dd$id) > 1))
dd2 <- dd[dd$id %in% ii, ]

disp <- with(dd2, tapply(ADAS11, id, var))/with(dd2, tapply(ADAS11, id, mean))
quantile(disp)

# "Dispersion factor"
(1 + sigma.old.way)^2
((1 + sigma.new)^2)[1]
exp(sigma)
