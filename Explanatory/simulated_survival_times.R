rm(list=ls())
source('../theme_csda.R')
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
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
                   simData(ntms = 10, gamma.x = eta, gamma.y = gamma, theta0 = t0, theta1 = t1)$survdat,
                   simplify = F
  )
  names(sets)[i] <- paste0('theta0 = ', t0, ', theta1 = ', t1)
}

survtimes <- lapply(sets, function(x){
  do.call(c, lapply(x, function(y) sort(y[y$cens==1, 'survtime'])))
})

step1 <- do.call(c, survtimes) %>% as.data.frame() %>% 
  tibble::rownames_to_column('id') %>% 
  mutate(id = stringr::str_extract(id, 'theta0\\s\\=\\s\\-\\d\\,\\stheta1\\s\\=\\s0\\.\\d')) 

names(step1) <- c('id', 'survtime')

step2 <- step1 %>% 
  as_tibble %>% 
  separate(id, c('theta0', 'theta1'), '\\,\\s') %>% 
  mutate(
    nu = as.numeric(str_sub(theta0, -2)),
    alpha = as.numeric(str_sub(theta1, -3))
  ) %>% 
  mutate(
    nu.lbl = paste0('log(nu)==', nu),
    alpha.lbl=paste0('alpha==', alpha)
  )

step2 %>% 
  filter(survtime != 9.1) %>% 
  ggplot(aes(x = survtime)) + 
  geom_histogram(aes(y =  ..density.. ), bins = 20, fill = 'white', colour = 1) + 
  geom_density(lwd = .5, colour = 'red', lty = 'dashed') +
  facet_grid(nu.lbl~alpha.lbl, scales = 'free', labeller = label_parsed, ) + 
  theme_csda() + 
  theme(
    strip.placement = 'outside',
    strip.text = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + 
  labs(x = 'Simulated failure time', y = '')

ggsave('output/alpha_nu_failure.png', width = 140, height = 90, units = 'mm')

# Changing c(gamma, eta) --------------------------------------------------
rm(list=ls())
source('theme_csda.R')
# Hold baseline hazard constant at nu = -4
theta0 <- -6; theta1 <- 0.001;
gamma.eta <- expand.grid(
  gamma = c(-1.0, -0.5, 0.0, 0.5, 1.0),
  eta1 = c(0.00),
  eta2 = c(-1.0, -0.3, 0.3, 1.0)
)
# D = diag(c(0.5^2, 0.25^2, 0.5^2, 0.25^2))

sets <- list()
for(i in 1:nrow(gamma.eta)){
  row <- gamma.eta[i,,drop=F]
  i.gamma <- c(row$gamma)
  i.zeta <- c(0, row$eta2)
  sets[[i]] <- replicate(100,
                         gmvjoint::simData(
                           n = 100, ntms = 10, fup = 9, zeta = i.zeta, gamma = i.gamma,
                           theta = c(-6, 1e-3),
                           beta = t(c(1, 0.1, 0.33, -0.5)),
                           sigma = c(.16),
                           family = list('gaussian'),
                           D = diag(rep(1,2)),
                         )$sur,
                         simplify = F
  )
  names(sets)[i] <- paste0('gamma == ',i.gamma, ', eta[2] == ', i.zeta[2])
}

survtimes <- lapply(sets, function(x){
  do.call(c, lapply(x, function(y) sort(y[y$stat==1, 'survtime'])))
})

step1 <- do.call(c, survtimes) %>% as.data.frame() 

step1$id <- do.call(c, lapply(names(survtimes), function(x) rep(x, each = length(survtimes[[x]]))))

names(step1) <- c('survtime', 'id')

step2 <- step1 %>% 
  as_tibble %>% 
  separate(id, c('gamma', 'eta'), '\\,\\s') %>% 
  mutate(
    temp.gamma = purrr::map(gamma, ~ strsplit(.x, '==')),
    temp.eta = purrr::map(eta, ~ strsplit(.x, '==')),
  ) %>% 
  unnest(c(temp.gamma, temp.eta)) %>% 
  mutate(g = purrr::map_chr(temp.gamma, el, 2),
         e = purrr::map_chr(temp.eta, el, 2)) %>% 
  select(-temp.gamma, -temp.eta) %>% 
  mutate_at(vars(g, e), as.numeric) %>% 
  arrange(g, e) %>% 
  mutate(
    gamma.lbl = forcats::fct_inorder(paste0('gamma == ', g)),
    eta.lbl = forcats::fct_inorder(paste0('zeta[2] == ', e))
  )

step2 %>% 
  filter(survtime != 9.1) %>% 
  ggplot(aes(x = survtime)) + 
  geom_histogram(aes(y = ..density..), bins = 20, colour = 1, fill = 'white') + 
  geom_density(lwd = 0.5, colour = 'red', lty = 'dashed') + 
  facet_grid(eta.lbl~gamma.lbl, scales = 'free', labeller = label_parsed) + 
  theme_csda() + 
  theme(
    strip.placement = 'outside',
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + 
  labs(x = 'Simulated failure time', y = '')

ggsave('gamma_zeta_failure2.png', width = 140, height = 90, units = 'mm')


# Altering magnitude of vech(D) -------------------------------------------
rm(list=ls())
source('theme_csda.R')
theta0 <- -5; theta1 <- 1e-3
gamma <- c(0.5, 0)
eta <- c(0.0, 0.0)

int.slope <- expand.grid(
  intercept = c(3, 1, .5, .5^2),
  slope = c(1, .5, .25, .25^2)
)

sets <- list()
for(i in 1:nrow(int.slope)){
  D <- diag(c(int.slope[i, 1], int.slope[i, 2], .01, .01))
  sets[[i]] <- replicate(100,
                         simData(ntms = 10, gamma.x = eta, gamma.y = gamma, theta0 = theta0, theta1 = theta1,
                                 D = D)$survdat,
                         simplify = F
  )
  names(sets)[i] <- paste0('D11 == ', int.slope[i, 1], ', D22 == ', int.slope[i, 2])
}

survtimes <- lapply(sets, function(x){
  do.call(c, lapply(x, function(y) sort(y[y$cens==1, 'survtime'])))
})

step1 <- do.call(c, survtimes) %>% as.data.frame() 

step1$id <- do.call(c, lapply(names(survtimes), function(x) rep(x, each = length(survtimes[[x]]))))

names(step1) <- c('survtime','id')

step2 <- step1 %>% 
  as_tibble %>% 
  separate(id, c('D11', 'D22'), '\\,\\s') %>% 
  mutate(
    temp.D11 = map(D11, ~ strsplit(.x, '==')),
    temp.D22 = map(D22, ~ strsplit(.x, '==')),
  ) %>% 
  unnest(c(temp.D11, temp.D22)) %>% 
  mutate(intercept = map_chr(temp.D11, el, 2),
         slope = map_chr(temp.D22, el, 2)) %>% 
  select(-temp.D11, -temp.D22) %>% 
  mutate_at(vars(intercept, slope), as.numeric) %>% 
  arrange(intercept, slope) %>% 
  mutate(
    intercept.lbl = fct_inorder(paste0('D[11] == ', intercept)),
    slope.lbl = fct_rev(fct_inorder(paste0('D[22] == ', slope)))
  )

step2 %>% 
  filter(survtime != 9.1) %>% 
  ggplot(aes(x = survtime)) + 
  geom_histogram(aes(y = ..density..), bins = 20, colour = 1, fill = 'white') + 
  geom_density(lwd = 0.5, colour = 'red', lty = 'dashed') + 
  facet_grid(slope.lbl~intercept.lbl, scales = 'free', labeller = label_parsed) + 
  theme_csda() + 
  theme(
    strip.placement = 'outside',
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + 
  labs(x = 'Simulated failure time', y = '')

ggsave('int_slope_failure.png', width = 140, height = 90, units = 'mm')
