.r <- function(){rm(list=ls());source(".Rprofile")}


# mi vs psi ---------------------------------------------------------------
.r()
gauss2 <- get.psi.mi2("gaussian")
poiss <- get.psi.mi2("poisson")
negbin <- get.psi.mi2("negbin")
Gammas <- get.psi.mi2("Gamma")
bin <- get.psi.mi2("binomial", tune = c(3.25, 5.00))

e <- T
while(e){
  genp <- try(get.psi.mi2("genpois"), silent = TRUE)
  e <- inherits(genp, "try-error")
}


fams <- as.data.frame(rbind(gauss2[[1]], poiss[[1]], negbin[[1]], Gammas[[1]], bin[[1]], genp[[1]]))
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


# This plot ---------------------------------------------------------------
fams %>%
  mutate(famname2 = purrr::map_chr(family, family.dir.name)) %>% 
  mutate_at("famname2", ~gsub("\\-Normal\\-Justification$", "", .x)) %>%
  mutate_at("famname2", ~gsub("isedpoi", "ised Poi", .x)) %>% 
  mutate_at("famname2", ~gsub("ivebin", "ive bin", .x)) %>% 
  ggplot(aes(x = as.factor(mi), y = psi, fill = as.factor(surv))) + 
    geom_boxplot(outlier.alpha = 0.33, lwd = 0.25, 
                 fatten = 2, outlier.size = 0.50) + 
    facet_wrap(~famname2, scales = 'free_y') + 
    scale_fill_manual(values = c("#FFEDA0", "#F03B20")) + 
    labs(x = expression(m[i]), y = expression(psi["  "*i])) + 
    theme_csda() + 
    theme(
      legend.position = "none",
      strip.text = element_text(vjust = 1)
    )
# Save...
ggsave(save.dir.file.path("psi-mi-families.png"), width = 140, height = 90, units = "mm")


# Now creating ellipse per family -----------------------------------------
(f <- .valid.families)
xx <- c("gauss2", "poiss", "bin", "genp", "negbin", "Gammas")
for(i in seq_along(f)){
  out <- paste0(save.dir.file.path(family.dir.name(f[i])), "/ellipse-",f[i],".png")
  # print(out)
  to.eval <- paste0("plot(", xx[i], "[[2]])")
  # print(to.eval)
  png(file = out, width = 140, height = 110, units = "mm", res = 500)
  eval(parse(text=to.eval))
  dev.off()
  cat(sprintf("----%s done; saved in %s----\n", f[i], out))
}

