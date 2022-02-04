rm(list=ls())
setwd('~/Documents/thesis-plots/')
source('theme_csda.R')
library(tidyverse)
library(joineRML)

# Changing baseline hazard ONLY -------------------------------------------

eta <- c(0,0)
gamma <- c(0,0)

thetas <- as.matrix(expand.grid(
  theta0 = c(-6, -5, -4, -3),
  theta1 = c(.10, .20, .30, .50)
))

sets <- list()
for(i in 1:nrow(thetas)){
  tt <- thetas[i, ,drop=T]
  t0 <- tt[1]; t1 <- tt[2]
  sets[[i]] <- replicate(100,
                   simData(ntms = 10, gamma.x = eta, gamma.y = gamma,, theta0 = t0, theta1 = t1)$survdat,
                   simplify = F
  )
  names(sets)[i] <- paste0('theta0 = ', t0, ', theta1 = ', t1)
}

survtimes <- lapply(sets, function(x){
  do.call(c, lapply(x, function(y) sort(unique(y[y$cens==1, 'survtime']))))
})

step1 <- do.call(c, survtimes) %>% as.data.frame() %>% 
  rownames_to_column('id') %>% 
  mutate(id = str_extract(id, 'theta0\\s\\=\\s\\-\\d\\,\\stheta1\\s\\=\\s0\\.\\d')) 

names(step1) <- c('id', 'survtime')

step2 <- step1 %>% 
  as_tibble %>% 
  separate(id, c('theta0', 'theta1'), '\\,\\s') %>% 
  mutate(
    nu = as.numeric(str_sub(theta0, -2)),
    alpha = as.numeric(str_sub(theta1, -3))
  ) %>% 
  mutate(
    nu.lbl = paste0('nu==', nu),
    alpha.lbl=paste0('alpha==', alpha)
  )

step2 %>% 
  filter(survtime != 9.1) %>% 
  ggplot(aes(x = survtime)) + 
  geom_histogram(bins = 24) + 
  facet_grid(nu.lbl~alpha.lbl, scales = 'free', labeller = label_parsed, ) + 
  theme_csda() + 
  theme(
    strip.placement = 'outside',
    strip.text = element_text(size = 11.5),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + 
  labs(x = 'Simulated failure time', y = '')

ggsave('alpha_nu_failure.png', width = 140, height = 90, units = 'mm')

# Changing c(gamma, eta) --------------------------------------------------
rm(list=ls())
# Hold baseline hazard constant at nu = -4



# Altering magnitude of vech(D) -------------------------------------------




