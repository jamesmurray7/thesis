# Showing longitudinal and survival parts ---------------------------------
PBC <- PBC[,c("id", "time", "survtime", "status", "serBilir", "albumin")]
PBC$sb <- log(PBC$serBilir)

Ls <- ggplot(PBC, aes(x = time, y = sb, group = id)) + 
  geom_line(lwd = 0.35, alpha = 0.25) + 
  geom_smooth(aes(group = NULL),se = F, method = "loess", 
              colour = nclred) + 
  facet_wrap(~ifelse(status==1, "Died", "Survived"),ncol=1,nrow=2) + 
  labs(x = "Follow-up time (years)", y = expression(log("Serum bilirubin")),
       title = "Longitudinal") +
  # scale_x_continuous(labels = seq(0,15,5), breaks = seq(0,15,5)) + 
  expand_limits(x=c(0,15)) + 
  theme_csda() + 
  theme(
    title = element_text(size=8.5, colour = nclblue),
    strip.text = element_text(size=7.0),
    axis.title = element_text(size=7.0, colour = "black"),
    axis.text = element_text(size = 6.5)
  )

# Make survdata
SD <- PBC[!duplicated(PBC$id), ][,c("id", "survtime", "status")]
SS <- survfit(Surv(survtime, status) ~ 1, data= SD,conf.type="plain")
se <- sqrt(cumsum(SS$n.event/(SS$n.risk*(SS$n.risk-SS$n.event))))

df <- data.frame(
  time = SS$time,
  S = SS$surv,
  Supp = pmin(1, SS$surv + 1.96 * se * SS$surv),
  Slow = SS$surv - 1.96 * se * SS$surv
)

Ss <- ggplot(df, aes(x = time, y = S)) + 
  geom_step(colour = nclred, lwd = 0.53) +
  geom_step(aes(y = Slow), lty = 5, colour = nclred, lwd = 0.35) + 
  geom_step(aes(y = Supp), lty = 5, colour = nclred, lwd = 0.35) + 
  theme_csda() + 
  expand_limits(y = c(0,1)) + 
  labs(x = "Follow-up time (years)", y = "Survival probability",
       title = "Time-to-event") + 
  theme(
    title = element_text(size=8.5, colour = nclblue),
    strip.text = element_text(size=7.0),
    axis.title = element_text(size=7.0, colour = 'black'),
    axis.text = element_text(size = 6.5)
  )

png(file = "./JM_visual.png", 
    width = linewidth, height =110, units = 'mm', res = 1000)
gridExtra::grid.arrange(Ls, Ss, nrow = 1)
dev.off()

# =========================================================================
# =========================================================================
# =========================================================================
# =========================================================================
# =========================================================================
# =========================================================================
# =========================================================================
# =========================================================================
# Other way around --------------------------------------------------------
# =========================================================================
PBC <- PBC[,c("id", "time", "survtime", "status", "serBilir", "albumin")]
PBC$sb <- log(PBC$serBilir)

Ls <- ggplot(PBC, aes(x = time, y = sb, group = id)) + 
  geom_line(lwd = 0.35, alpha = 0.25) + 
  geom_smooth(aes(group = NULL),se = F, method = "loess", 
              colour = nclred) + 
  facet_wrap(~ifelse(status==1, "Died", "Survived"),ncol=2,nrow=1) + 
  labs(y = expression(log("Serum bilirubin")), x = NULL,
       title = "Longitudinal") +
  # scale_x_continuous(labels = seq(0,15,5), breaks = seq(0,15,5)) + 
  expand_limits(x=c(0,15)) + 
  theme_csda() + 
  theme(
    title = element_text(size = 25, colour = nclblue),
    strip.text = element_text(size=21),
    axis.title = element_text(size=19, colour = "black"),
    axis.text = element_text(size = 15)
  )

# Make survdata
SD <- PBC[!duplicated(PBC$id), ][,c("id", "survtime", "status")]
SS <- survfit(Surv(survtime, status) ~ 1, data= SD,conf.type="plain")
se <- sqrt(cumsum(SS$n.event/(SS$n.risk*(SS$n.risk-SS$n.event))))

df <- data.frame(
  time = SS$time,
  S = SS$surv,
  Supp = pmin(1, SS$surv + 1.96 * se * SS$surv),
  Slow = SS$surv - 1.96 * se * SS$surv
)

Ss <- ggplot(df, aes(x = time, y = S)) + 
  geom_step(colour = nclred, lwd = 0.53) +
  geom_step(aes(y = Slow), lty = 5, colour = nclred, lwd = 0.35) + 
  geom_step(aes(y = Supp), lty = 5, colour = nclred, lwd = 0.35) + 
  theme_csda() + 
  expand_limits(y = c(0,1)) + 
  labs(x = "Follow-up time (years)", y = "Survival probability",
       title = "Time-to-event") + 
  theme(
    title = element_text(size = 25, colour = nclblue),
    strip.text = element_text(size=21),
    axis.title = element_text(size=19, colour = 'black'),
    axis.text = element_text(size = 15)
  )

png(file = "./JM_visual_otherway.png", 
    width = linewidth, height = 200, units = 'mm', res = 1000)
gridExtra::grid.arrange(Ls, Ss, nrow = 2)
dev.off()

