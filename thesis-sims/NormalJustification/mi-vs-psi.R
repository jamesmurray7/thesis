.r <- function(){rm(list=ls());source(".Rprofile")}


# mi vs psi ---------------------------------------------------------------
.r()
gauss2 <- get.psi.mi2("gaussian")
poiss <- get.psi.mi2("poisson")
negbin <- get.psi.mi2("negbin")
Gammas <- get.psi.mi2("Gamma")
bin <- get.psi.mi2("binomial")
genp <- get.psi.mi2("genpois")

fams <- as.data.frame(rbind(gauss[[1]], poiss[[1]], negbin[[1]], Gammas[[1]], bin[[1]], genp[[1]]))
library(dplyr)
library(ggplot2)
fams2 <- fams %>% 
  mutate_at("family", forcats::fct_inorder) %>% 
  group_by(family, mi) %>% 
  summarise(mn = mean(psi), sd = sd(psi),
            med = median(psi), X25 = quantile(psi, .25), X75 = quantile(psi, .75),
            .groups = "keep") %>% 
  ungroup

fams2

ggplot(fams2, aes(x = as.factor(mi), y = med, colour = family)) + 
  geom_point() 

ggplot(fams, aes(x = as.factor(mi), y = psi, colour = family, group = family)) + 
  geom_smooth(aes()) + 
  facet_wrap(~surv, scales='free_y')
  

ggplot(fams, aes(x = as.factor(mi), y = psi, colour = family)) + 
  geom_boxplot()

sim <- create.simfn(list("binomial"), arguments = list(n=1000, D = diag(c(2,.5)), random.formulas = list(~time)))
(X_ <- dataGen(sim))
S <- Sample(X_, 5.25, T, T, T, NMC = 5000)
pbmi <- prop.by.mi(S,.2,.25)

ggplot(pbmi, aes(x = as.factor(mi), y = psi, group = surv)) + 
  geom_smooth(aes())

ggplot(pbmi, aes(x = as.factor(mi), y = psi)) + 
  geom_boxplot()
