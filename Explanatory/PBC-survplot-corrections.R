library(gmvjoint)
data(PBC)
library(dplyr)
library(ggplot2)
source('../theme_csda.R')

# LHS: KM plot ------------------------------------------------------------
survdata <- PBC %>% distinct(id, survtime, status)
S <- survfit(Surv(survtime, status) ~ 1, survdata)
SS <- summary(S)
Sdf <- data.frame(
  t = SS$time,
  surv = SS$surv,
  lower = SS$lower, upper = SS$upper
)

P1 <- ggplot(Sdf, aes(x = t, y = surv)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper), fill = .nice.orange, alpha = .22, colour = 'grey90', lwd = .25)+
  geom_step(colour = 'black', lwd=.15) + 
  expand_limits(y = 0:1, x = 0:15.5) +
  scale_x_continuous("Follow-up time (years)", breaks=0:15) +
  labs(y = "Survival probability") + 
  theme_csda() + 
  theme(
    axis.text = element_text(size=5.5),
    axis.title = element_text(size=7),
    axis.ticks.length.x =  unit(.5,'mm')
  )

# RHS: Histogram ----------------------------------------------------------
hist(survdata %>% filter(status==1) %>% pull(survtime), breaks = 15)

survdata %>% 
  filter(status == 1) %>% 
  ggplot(aes(x = survtime)) + 
  geom_histogram(bins = 15, fill = "grey90", colour = .nice.orange,lwd=.25) +
  labs(x = "Follow-up time (years)",
       y = "Frequency") + 
  scale_x_continuous(breaks = seq(0,15,1)) + 
  theme_csda()+ 
  theme(
    axis.text = element_text(size=5.5),
    axis.title = element_text(size=7),
    axis.ticks.length.x =  unit(.5,'mm')
  ) -> P2

ggpubr::ggarrange(P1,P2,nrow=1)
ggsave("./output/PBC-KM2-correction.png", width = 140, height=60, units = "mm")
