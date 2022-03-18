library(tidyverse)
source('~/Documents/Bernhardt/theme_csda.R')
pbc <- joineRML::pbc2

# mean and variances for platelet count
means <- with(pbc, tapply(platelets, id, mean))
vars <-  with(pbc, tapply(platelets, id, var))
times <- with(pbc, tapply(year, id, FUN = function(x){ length(unique(x))} ))
status <- with(pbc, tapply(status2, id, unique))
ids <- unique(pbc$id)


df <- data.frame(
  m = means[!is.na(means)],
  v = vars[!is.na],
  id = names(means[!onetp])
)

df$vrat = df$v/df$m

ggplot(df, aes(id, vrat, colour = status[!onetp] == 1)) + 
  geom_hline(aes(yintercept = 1), lty = 5, colour = 'gray50') + 
  geom_point(size = .75) + 
  labs(
    x = NULL,
    y = expression('Var(x)/mean(x)'),
    colour = 'Failed'
  ) + 
  scale_colour_manual(values = c('black', 'red')) + 
  scale_y_continuous(breaks = c(0, 1, seq(10, 100, 10)))+
  theme_csda() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) + coord_cartesian(ylim = c(0,100))

kruskal.test(df$vrat~status[!onetp])
wilcox.test(df$vrat~status[!onetp])
