rm(list=ls());library(tidyverse)
source('..//theme_csda.R')
kidneys <- joineRML::renal

kidneys$surv %>% distinct(id, failure) %>% count(failure) %>% mutate(p=prop.table(n)*100) #41% failure

lapply(kidneys, names)

# Plot the two continuous biomarkers on the same plot...

gfr <- kidneys$gfr %>% select(id, `level` = gfr, time = years, status = failure, survtime = fuyears) %>% mutate(bio = "Glomerular Filtration Rate (GFR)")
hae <- kidneys$haem %>% select(id, `level` = haematocrit, time = years, status = failure, survtime = fuyears) %>% mutate(bio = "Blood Haematocrit Level")

k <- rbind(gfr, hae)

ktib <- k %>% 
  filter(status == 1) %>% 
  filter(complete.cases(.)) %>% 
  as_tibble %>% 
  mutate(tt = -1 * (survtime - time))
  
ktib %>% 
  ggplot(aes(x = tt, y = level, group = id)) + 
  geom_vline(xintercept = 0, alpha = .25) + 
  geom_line(alpha = .25, lwd = .1) + 
  geom_smooth(aes(group = NULL), method = 'loess', formula = y ~ x, colour = 'black') + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  facet_wrap(~bio, strip.position = 'left', scales = 'free_y') +
  labs(y = NULL,
       x = 'Time (years) from kidney graft failure (0: failure time)') + 
  theme_csda() + 
  theme(strip.placement = 'outside')

ggsave('./KidneyContinuousTrajectories.png', width = 140, height = 60, units = 'mm', dpi = 1e3)  

# Binary covariate (proterinuria) -----------------------------------------

prot <- kidneys$prot %>% as_tibble %>% 
  filter(failure == 1) %>% 
  filter(complete.cases(.)) %>% 
  mutate(tt = -1 * (fuyears - years))

# Counting number of changes across IDs

change20s <- prot %>% 
  group_by(id) %>% 
  mutate(changed  = abs(c(NA, diff(proteinuria)))) %>% 
  mutate(sum.changed = sum(changed, na.rm = T)) %>% ungroup %>% 
  distinct(id, sum.changed) %>% # pull(sum.changed) %>% sort %>% tail
  filter(sum.changed > 25) %>% pull(id)

prot %>% 
  filter(id %in% change20s) %>% 
  group_by(id) %>% 
  mutate(id2 = cur_group_id(), p =factor(proteinuria, c(0, 1), c('Absent', 'Present'))) %>% 
  ggplot(aes(x = tt, y = id2, colour = p)) + 
  geom_vline(xintercept = 0, alpha = .25) +
  geom_point(size=.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(length(change20s))) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  scale_color_manual(values = c('black', 'red')) + 
  theme_csda() + 
  labs(x = 'Time (years) from kidney graft failure (0: failure time)', 
       y = 'ID', colour = 'Proteinuria') + 
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave('./KidneyBinaryTrajectories.png', width = 140, height = 60, units = 'mm', dpi = 1e3)


