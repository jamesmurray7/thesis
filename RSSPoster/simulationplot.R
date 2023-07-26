# RIPPED FROM gmvjoint/paper-sims/Simulation/parse.R
rm(list=ls())
source(".Rprofile")
load("/data/c0061461/GLMM_Paper_Sims/Revision2/fits_regspaceN200_new.RData")
vech <- gmvjoint:::vech
# Restate true parameter values
N <- 200
D.true <- diag(c(0.25, 0.09, 0.50, 0.10, 2.00))
D.true[1,3] <- D.true[3,1] <- D.true[1,5] <- D.true[5,1] <- D.true[3,5] <- D.true[5,3] <- 0.25
beta.true <- c(2, -0.1, 0.1, -0.2, 2, -0.1, 0.1, -0.2, 1, -1, 1, -1)
sigma.true <- .16
gamma.true <- c(.5,-.5,.5)
zeta.true <- -.2
target <- c(vech(D.true), beta.true, sigma.true, gamma.true, zeta.true)
target.mat <- t(apply(t(matrix(target)),2,rep,N))
# Parsing fits.
qz <- qnorm(.975)
# Function to extract estimates
extract.estimates <- function(L) sapply(L, function(x){
  setNames(c(vech(x$coeffs$D), x$coeffs$beta, unlist(x$coeffs$sigma)[1], x$coeffs$gamma, x$coeffs$zeta),
           names(x$SE))
})
# Function to extract SEs
extract.SE <- function(L) sapply(L, function(x){
  x$SE
})

ests <- lapply(fits, extract.estimates)
SEs <- lapply(fits, extract.SE)
times <- lapply(fits, function(x) sapply(x, function(y) y$elapsed))

# Empirical Means, SDs and average estimated SE.
emp.mean <- lapply(ests, rowMeans)
emp.sd <- lapply(ests, function(x) apply(x, 1, sd))
avg.se <- lapply(SEs, rowMeans)

# 95% coverage probability
CPs <- Map(function(estimate, se){
  lb <- estimate - qz * se; ub <- estimate + qz * se
  rowSums(lb <= target.mat & ub >= target.mat)/N
}, estimate = ests, se = SEs)

MSEs <- Map(function(estimate){
  rowMeans((target.mat - estimate)^2)
}, estimate = ests)

Biases <- Map(function(estimate){
  rowMeans(estimate-target.mat)
}, estimate = ests)

Elapseds <- Map(function(TIME){
  EM <- TIME[1,]
  PP <- TIME[2,]
  EM + PP
}, TIME = times)

Totals <- Map(function(TIME){
  TIME[3,]
}, TIME = times)

iter.per.second <- Map(function(TIME){
  TIME[4,]/TIME[1,]
}, TIME = times)

df.Elapsed <- expand.grid(r = c(5, 10, 15),
                          omega = c(0.1, 0.3, 0.5))
df.Elapsed <- cbind(df.Elapsed, t(apply(df.Elapsed,1,function(x){
  r <- x[1]; fail <- paste0(100*x[2], "%")
  .lookup <- paste0("n = 250, mi = ", r, ", failure = ", fail)
  qn <- quantile(Elapseds[[.lookup]])
  c(qn[2], qn[3], qn[4])
})))

df.Totals <- expand.grid(r = c(5, 10, 15),
                          omega = c(0.1, 0.3, 0.5))
df.Totals <- cbind(df.Totals, t(apply(df.Totals,1,function(x){
  r <- x[1]; fail <- paste0(100*x[2], "%")
  .lookup <- paste0("n = 250, mi = ", r, ", failure = ", fail)
  qn <- quantile(Totals[[.lookup]])
  c(qn[2], qn[3], qn[4])
})))

sapply(iter.per.second, mean)

# Plotting survival parameters --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
source('../theme_csda.R')

all.ests <- do.call(rbind, lapply(seq_along(ests), function(x){
  current <- names(ests)[x]
  mi <- stringr::str_extract(current, 'mi\\s\\=\\s\\d?\\d') %>% 
        gsub("mi\\s\\=\\s",'',.) %>% as.integer
  rate <- stringr::str_extract(current, "\\d\\d\\%$") %>% 
        gsub("\\%", "", .) %>% as.integer
  
  survs <- as.data.frame(t(ests[[x]][(nrow(ests[[x]])-3):nrow(ests[[x]]),]))
  info <- data.frame(r = rep(mi, nrow(survs)),
                     omega = rep(rate/100, nrow(survs)))
  cbind(info, survs)
}))


# plot of estimates for gamma + zeta --------------------------------------

all.ests %>% 
  pivot_longer(`gamma_1`:`zeta_bin`) %>% 
  mutate(name = case_when(
    name == "gamma_1" ~ "gamma[1]",
    name == "gamma_2" ~ "gamma[2]",
    name == "gamma_3" ~ "gamma[3]",
    T ~ "zeta"
  )) %>% 
  mutate(target = case_when(
    name == "gamma[1]" ~ .5,
    name == "gamma[2]" ~ -.5,
    name == "gamma[3]" ~ .5,
    T ~ -.2
  )) %>% 
  mutate(r = factor(r, levels = c(5,10,15)),
         omega = factor(omega, levels = c(.1, .3, .5), labels = c("10%", "30%", "50%"))) %>% 
  ggplot(aes(x = omega, y = value, fill = r)) + 
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey5', alpha = .5) + 
  geom_boxplot(outlier.alpha = .33, lwd = 0.35, fatten = 2,
               outlier.size = .50) + 
  facet_wrap(~name, scales = 'free_y', labeller = label_parsed) + 
  labs(fill = "Maximal profile length", 
       #x = "Approximate failure rate",
       x = NULL,
       y = "Estimate") + 
  theme_csda() + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme(
    legend.position = "none",
    strip.text = element_text(size = 21, colour = nclred),
    axis.title = element_text(size = 19, colour = 'black'),
    axis.text = element_text(size = 15)
  ) -> q

G <- ggplotGrob(q)
strip.titles <- which(sapply(G$grobs, function(x) grepl("strip", x$name))) # zeta must be last one?
sapply(strip.titles, function(i) G$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label) # no, its second for some reason
# => We want to change the second one!
G$grobs[[strip.titles[2]]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "black"
gam.zet.plot <- ggpubr::as_ggplot(G)

# Plot of elapsed times ---------------------------------------------------

elapsed.df <- do.call(rbind, lapply(seq_along(Elapseds), function(i){
  current <- names(Elapseds)[i]
  mi <- stringr::str_extract(current, 'mi\\s\\=\\s\\d?\\d') %>% 
    gsub("mi\\s\\=\\s",'',.) %>% as.integer
  rate <- stringr::str_extract(current, "\\d\\d\\%$") %>% 
    gsub("\\%", "", .) %>% as.integer
  cat(sprintf("current: %s; r: %d, omega: %d\n", current, mi, rate))
  elapsed <- Elapseds[[i]]
  data.frame(elapsed = elapsed, r = factor(mi, levels = c(5,10,15)), omega = factor(rate, c("10", "30", "50"),
                                                                                    c("10%", "30%", "50%")))
}))

elapsed.df %>% 
  ggplot(aes(x = omega, y = elapsed, fill = r)) + 
  geom_boxplot(outlier.alpha = .33, lwd = 0.35, fatten = 2,
               outlier.size = .50) + 
  theme_csda() + 
  labs(y = "Elapsed time (s)",
       x = NULL,
       fill = "Length of follow-up") + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme(
    legend.position = "right",
    strip.text = element_text(size = 21, colour = nclred),
    axis.title = element_text(size = 19, colour = 'black'),
    axis.text = element_text(size = 15),
    legend.title = element_text(angle = 270, size = 19),
    legend.text = element_text(angle = 270, hjust = 0.5,vjust=0, size=  15),
    legend.key.height = unit(13, "mm"),
    legend.box.spacing = unit(1,'mm')
  ) -> elapsed.plot 

gridExtra::grid.arrange(gam.zet.plot, elapsed.plot, nrow = 1,
                        bottom = grid::grid.text("Approximate failure rate", gp = grid::gpar(cex=1.5)))

png("./simulation.png", width = linewidth2, height = 150, units = "mm", res = 1500)
gridExtra::grid.arrange(gam.zet.plot, elapsed.plot, nrow = 1,
                        bottom = grid::grid.text("Approximate failure rate", gp = grid::gpar(cex=1.5)))
dev.off()
