.r <- function(){rm(list=ls());gc();source(".Rprofile")}
.xyz <- function(f) file.path(save.dir.file.path(family.dir.name(f)), paste0(f,"-psi-new.RData"))

# Gaussian ----------------------------------------------------------------
.r()
sim <- create.simfn(list("gaussian"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 3, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("gaussian"))

# Poisson -----------------------------------------------------------------
.r()
sim <- create.simfn(list("poisson"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 3, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("poisson"))


# Binomial ----------------------------------------------------------------
.r()
sim <- create.simfn(list("binomial"), arguments = list(n = 500L, theta = c(-3,.1), beta = t(c(2.0,-0.1,0.1,0.2)),
                                                       D = matrix(c(.5, .125, .125, .09), 2, 2)))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 2, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("binomial"))

# Gamma -------------------------------------------------------------------
.r()
sim <- create.simfn(list("Gamma"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 3, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("Gamma"))


# Negbin ------------------------------------------------------------------
.r()
sim <- create.simfn(list("negbin"), arguments = list(n = 500L, theta = c(-3,.1)))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 3, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("negbin"))

# Genpois -----------------------------------------------------------------
.r()
sim <- create.simfn(list("genpois"), arguments = list(n = 500L, theta = c(-3,.1),
                                                      D = matrix(c(0.5, .125, .125, .05), 2, 2),
                                                      sigma = list(0.2), beta = t(c(.5,-.2,.05,.4))))
X_ <- dataGen(sim)
(S <- Sample(X_, TUNE = 3, T, T, T, "t", 4, 1000L, 5000L))
psis <- prop.by.mi(S, .2, .3)
save(psis, file = .xyz("genpois"))

# Plotting all ------------------------------------------------------------
dirs <- fs::dir_ls(save.dir, type = 'directory')
dirs <- dirs[grepl("Justification$", dirs)]
out <- vector("list", length(dirs))
p <- 1
for(d in dirs){
  f <- gsub("\\/data\\/c0061461\\/THESIS\\/", "", d)
  f <- stringr::str_extract(f, '\\w+')
  cat(sprintf("\nDirectory: %s;\nFamily: %s", d, f))
  
  if(f == "Generalisedpoisson") f <- "Generalised Poisson"
  if(f == "Negativebinomial") f <- "Negative binomial"
  
  to.load <- dir(d, pattern = 'new\\.RData', full.names = T)
  assign("xxx", get(load(to.load)))
  
  xxx$f <- f
  
  out[[p]] <- xxx
  cat("\n")
  rm(xxx)
  p <- p + 1
  cat(cli::rule(col = 'red'))
}

out <- do.call(rbind, out)

library(ggplot2)
ggplot(out, aes(x = mi, y = psi, group = mi)) + 
  facet_wrap(~f,scales='free') + 
  geom_boxplot(outlier.alpha = 0.33, lwd = 0.25, 
               fatten = 2, outlier.size = 0.50, fill = .nice.orange, colour = 'black') + 
  theme_csda() + 
  scale_x_continuous(expression(m[i]), breaks = 1:10) + 
  scale_y_continuous(expression(psi[~~i])) + 
  theme(
    strip.text = element_text(size=7,vjust=1)
  )

ggsave(save.dir.file.path("psi-vs-mi-new.png"), width = 140, height = 90, units = "mm")


# ########################################
# THE BELOW IS OLD                       #
# ########################################

# mi vs psi ---------------------------------------------------------------
.r()
gauss <- get.psi.mi2("gaussian",file.name = "test.RData")
poiss <- get.psi.mi2("poisson", file.name = "test.RData")
negbin <- get.psi.mi2("negbin", file.name = "test.RData")
Gammas <- get.psi.mi2("Gamma", file.name = "test.RData")
bin <- get.psi.mi2("binomial", file.name = "test.RData", tune = c(2.0,2.0))

e <- T
while(e){
  genp <- try(get.psi.mi2("genpois"), silent = TRUE)
  e <- inherits(genp, "try-error")
}


fams <- as.data.frame(rbind(gauss[[1]], poiss[[1]], negbin[[1]], Gammas[[1]], bin[[1]])) #, genp[[1]]))
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

