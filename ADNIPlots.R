rm(list=ls()); library(tidyverse)
load("ADNI.RData")
source("./theme_csda.R")

# Re-scale ----------------------------------------------------------------
rescale <- function(x){
  sigma <- attr(x, "scaled:scale"); mu <- attr(x, "scaled:center")
  x * sigma + mu
}

adni <- mutate(adni, across(ADAS13:MidTemp, ~ rescale(.x)))

adni %>% distinct(id) %>% tally
adni %>% distinct(id, status) %>% colSums()

# Plot --------------------------------------------------------------------

adni2 <- adni %>% 
  filter(status == 1) %>% 
  select(PTGENDER, id, ADAS13:MidTemp, Year, survtime) %>% 
  gather("measure", "value", -id, -PTGENDER, -Year, -survtime) %>% 
  arrange(id) %>% 
  mutate(timefrom = -1 * (survtime - Year))


# Emulating Li et al. (2017) ----------------------------------------------
profplot <- adni2 %>% 
  mutate(mm = case_when(
    measure == "FAQ" ~ "FAQ score",
    measure == 'MidTemp' ~ "Middle temporal gyrus volume",
    measure == "ADAS13" ~ "ADAS13 score"
  )) %>% 
  ggplot(aes(x = timefrom, y = value, group = id)) + 
  geom_vline(xintercept = 0, alpha = .25, colour = 'black') + 
  geom_line(lty = 1, alpha = .25, lwd = .25) + 
  facet_wrap(~mm, scales = 'free_y', strip.position = 'left') + 
  theme_csda() + 
  geom_smooth(aes(group = NULL), colour = 'black',
              method = 'loess', formula = y~x) + 
  labs(x = "Time (years) from Alzheimer's conversion (0: diagnosis)",
       y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))

profplot
ggsave("~/Downloads/ADNItrajectories.png", width = 140, units = 'mm', height = 60)



