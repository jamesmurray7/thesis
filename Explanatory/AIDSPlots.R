rm(list=ls());library(tidyverse)
source('../theme_csda.R')
aids <- joineR::aids

# Those who failed only ---------------------------------------------------
aids %>% 
  as_tibble %>% 
  filter(death == 1) %>% 
  mutate(tt = -1 * (time - obstime)) %>% 
  ggplot(aes(x = tt, y = CD4, group = id)) + 
  geom_line(lwd=.25,alpha=.25) + 
  geom_smooth(aes(group=NULL), method = 'loess', formula = y~x, colour = 'black') +
  theme_csda()

# Might be too barren of a plot, comparing all trajectories?
aids %>% distinct(id, drug) %>% count(drug)


aids %>% 
  as_tibble %>% 
  ggplot(aes(x = obstime, y = CD4, group = id)) + 
  geom_line(lwd=.25,alpha=.25) + 
  geom_smooth(aes(group = as.factor(drug)), method = 'lm', formula = y ~ x)+
  geom_smooth(aes(group = as.factor(drug)), method = 'lm', formula = y ~ poly(x^2)) +
  theme_csda()

# Dont think this worth pursuing.