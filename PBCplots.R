library(tidyverse)
source('~/Documents/Bernhardt/theme_csda.R')
pbc <- joineRML::pbc2

# Plot of serBilir, Albumin and prothrombin time
pbc2 <- pbc %>% 
  select(id, year, years, status2, drug, serBilir, albumin, prothrombin) %>% 
  filter(complete.cases(.)) %>% 
  as_tibble %>% 
  mutate(tt = -1 * (years-year))

pbc2 %>% 
  filter(status2 == 1) %>% 
  pivot_longer(cols=serBilir:prothrombin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ 'Serum bilirubin (mg/dL)',
    biomarker == 'albumin' ~ 'Albumin (g/dL)',
    biomarker == 'prothrombin' ~ 'Prothrombin time (sec)'
  )) %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .25) + 
  geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) + 
  geom_smooth(aes(group=NULL), colour = 'red', method = 'lm', formula = y~x) + 
  geom_smooth(aes(group=NULL), colour = 'blue', method = 'lm', formula = y~x+I(x^2)) + 
  facet_wrap(~biomarker, scales = 'free', strip.position = 'left') + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))

pbc2 %>% 
  filter(status2 == 1) %>% 
  mutate(
    serBilir = log(serBilir),
    prothrombin = (.1 * prothrombin)^-4
  ) %>% 
  pivot_longer(cols=serBilir:prothrombin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ 'log(Serum~bilirubin)',
    biomarker == 'albumin' ~ 'Albumin',
    biomarker == 'prothrombin' ~ '(0.1~x~Prothrombin~time)^{-4}'
  ),
  biomarker = factor(biomarker, c('log(Serum~bilirubin)', 'Albumin', '(0.1~x~Prothrombin~time)^{-4}'))) %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .25, lwd = .25) + 
  geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) + 
  #  geom_smooth(aes(group=NULL), colour = 'red', method = 'lm', formula = y~x) + 
  #  geom_smooth(aes(group=NULL), colour = 'blue', method = 'lm', formula = y~x+I(x^2)) + 
  facet_wrap(~biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))
ggsave('~/Downloads/LOESSPBC.png', width = 140, height = 60, units = 'mm')


# Working with transformed trajectories -----------------------------------
pbc2 %>% 
  filter(status2 == 1) %>% 
  mutate(
    serBilir = log(serBilir),
    prothrombin = (.1 * prothrombin)^-4
  ) %>% 
  pivot_longer(cols=serBilir:prothrombin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ 'log(Serum~bilirubin)',
    biomarker == 'albumin' ~ 'Albumin',
    biomarker == 'prothrombin' ~ '(0.1~x~Prothrombin~time)^{-4}'
  ),
  biomarker = factor(biomarker, c('log(Serum~bilirubin)', 'Albumin', '(0.1~x~Prothrombin~time)^{-4}'))) %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .25, lwd = .25) + 
  geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) + 
  #  geom_smooth(aes(group=NULL), colour = 'red', method = 'lm', formula = y~x) + 
  #  geom_smooth(aes(group=NULL), colour = 'blue', method = 'lm', formula = y~x+I(x^2)) + 
  facet_wrap(~biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))
ggsave('~/Downloads/app1fig.png', width = 140, height = 60, units = 'mm')
tiff('~/Downloads/app1fig.tiff', width = 140, height = 60, units = 'mm', compression = 'lzw', res = 1e3)
last_plot()
dev.off()


# Stacked columns for internal poster -------------------------------------
pbc2 %>% filter(status2 == 1) %>% distinct(id) %>% nrow / 312 * 100

pbc2 %>% 
  filter(status2 == 1) %>% 
  pivot_longer(cols=serBilir:albumin, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ 'Serum bilirubin (mg/dL)',
    biomarker == 'albumin' ~ 'Albumin (g/dL)',
  )) %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .25, lwd = .33) + 
  geom_smooth(aes(group=NULL), colour = 'black', method = 'loess', formula = y~x) + 
  facet_wrap(~biomarker, scales = 'free', strip.position = 'left', nc = 1) + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1),
        axis.title = element_text(size=7))

ggsave('~/Downloads/LOESSPBC.png', width = 1015.9156 / 47 * 2.85426, height = 1.618033988749895 * ( 1015.9156 / 47 * 2.85426), units = 'mm', dpi = 1e3)
