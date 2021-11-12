rm(list=ls());library(tidyverse)
source('~/Documents/Bernhardt/theme_csda.R')
kidneys <- joineRML::renal

kidneys$surv %>% distinct(id, failure) %>% count(failure) %>% mutate(p=prop.table(n)*100) #41% failure

lapply(kidneys, names)
# Define a plotting function

plotter <- function(data, X, lab){
  plot.data <- data %>% 
    select(id, {{X}}, years, survtime = fuyears, failure) %>% 
    filter(failure == 1) %>% 
    filter(complete.cases(.)) %>% 
    as_tibble %>% 
    mutate(tt = -1 * (survtime - years))
  
  plot.data %>% 
    ggplot(aes(x = tt, y = {{X}}, group = id)) + 
    geom_vline(xintercept = 0, alpha = .25) + 
    geom_line(alpha = .25, lwd = .25) + 
    geom_smooth(aes(group = NULL), method = 'loess', formula = y ~ x, colour = 'black') + 
    scale_y_continuous(breaks = scales::pretty_breaks(10)) +
    scale_x_continuous(breaks = scales::pretty_breaks(10)) +
    labs(y = lab,
         x = 'Time (years) from kidney graft failure (0: failure time)') + 
    theme_csda()
}

gfrplot <- plotter(kidneys$gfr, gfr, 'Glomerular Filtration Rate (GFR)')
haeplot <- plotter(kidneys$haem, haematocrit, 'Blood Haematocrit Level')

# Some way of plotting the binary one

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
  filter(sum.changed > 20) %>% pull(id)

prot %>% 
  filter(id %in% change20s) %>% 
  group_by(id) %>% 
  mutate(id2 = cur_group_id(), p =factor(proteinuria, c(0, 1), c('Absent', 'Present'))) %>% 
  ggplot(aes(x = tt, y = id2, colour = p)) + 
  geom_vline(xintercept = 0, alpha = .25) +
  geom_point(size=.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(18)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  scale_color_manual(values = c('black', 'red')) + 
  theme_csda() + 
  labs(x = 'Time (years) from kidney graft failure (0: failure time)', 
       y = 'ID', colour = 'Proterinuria')

prot <- last_plot()


ggpubr::ggarrange(gfrplot, haeplot, prot)
