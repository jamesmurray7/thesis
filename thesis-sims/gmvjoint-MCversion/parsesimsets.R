rm(list=ls())
source(".Rprofile")
# Load in appropriate one ----
# load(save.dir.file.path("fitsMC4-nosf.RData"))
# load(save.dir.file.path("fitsMC4.RData"))
# load(save.dir.file.path("fitsGH.RData"))
load(save.dir.file.path("fitsGH-n500"))
fits <- GHfits
load(save.dir.file.path("fitsMC4-nosf-n500.RData"))

Antifits <- lapply(MCfits, el, 2)
Quasfits <- lapply(MCfits, el, 3)

# Objects are called: fits {GH}; Antifits {Antithetic MC}; Quasfits {Sobol sequence}

# Making target vector and matrix
D <- diag(c(.25,.09,.50,.09,2))
D[1,3] <- D[3,1] <- D[1,5] <- D[5,1] <- D[3,5] <- D[5,3] <- .125
targets <- setNames(
  c(gmvjoint:::vech(D), 
    c(-2,.1,-.1,.2),
    c(2,-.1,.1,.2),
    c(1,-1,1,-1),
    .16,
    .5, -.5, .5,
    -0.2
  ),
  names(fits[[1]]$SE)
)
targets.mat <- t(apply(t(targets), 2, rep, NROW(fits)))

# Function to extract
extract.from.joint <- function(x){
  stopifnot(inherits(x, "joint"))
  qz <- qnorm(.975)
  K <- 3
  co <- x$coeffs
  Omega <- setNames(c(gmvjoint:::vech(co$D), co$beta, co$sigma[[1]], co$gamma, co$zeta),
                    names(x$SE))
  SE <- x$SE
  lb <- Omega - qz * SE; ub <- Omega + qz * SE
  et <- x$elapsed.time
  return(list(
    Omega = Omega, SE = SE, lb = lb, ub = ub, 
    elapsed = et[1] + et[2], total = et[3], iters = et[4]
  ))
}

# Function to parse and make things to tabulate.
parser <- function(X){
  X.parsed <- lapply(X, extract.from.joint)
  
  Omega <- sapply(X.parsed, '[[', "Omega")
  SE <- sapply(X.parsed, '[[', "SE")
  lb <- sapply(X.parsed, '[[', "lb")
  ub <- sapply(X.parsed, '[[', "ub")
  et <- sapply(X.parsed, '[[', "elapsed")
  tot<- sapply(X.parsed, '[[', "total")
  iters <- sapply(X.parsed, '[[', "iters") 
  
  Mean <- rowMeans(Omega)
  SD <- apply(Omega, 1, sd)
  MeanSE <- rowMeans(SE)
  Bias <- rowMeans(Omega-targets.mat)
  MSE <- rowMeans((targets.mat-Omega)^2)
  CP <- rowMeans(lb < targets.mat & ub > targets.mat)
  
  return(list(Mean = Mean, SD = SD, MeanSE = MeanSE,
         Bias = Bias, MSE = MSE, CP = CP,
         et = unname(et), 
         tot = unname(tot), iters = unname(iters)))
}

GHparsed <- parser(fits)
ATparsed <- parser(Antifits)
QIparsed <- parser(Quasfits)

# Now get into xtable-ready format
.f <- function(x, d = 3) format(round(x, d), nsmall = d)
.qn <- function(x) quantile(x, c(.5, .25, .75))
to.tab <- function(X, method = "GH"){
  params <- names(X$Mean)
  # The table
  msd <- paste0(.f(X$Mean), ' (', .f(X$SD), ')')
  mSE <- .f(X$MeanSE)
  bias <- .f(X$Bias); MSE <- .f(X$MSE); CP <- .f(X$CP, 2)
  tab <- data.frame(
    Parameter = params,
    "Mean (SD)" = msd, `SE` = mSE, Bias = bias, MSE = MSE, CP = CP,
    row.names = NULL
  )
  # Elapsed time info
  et <- .qn(X$et); tot <- .qn(X$tot); iters <- .qn(X$iters)
  cap <- paste0(
    sprintf("Joint models fit by %s had median [IQR] elapsed time for model convergence and calculation standard errors ",
            switch(method,
                   GH = "Gauss-Hermite quadrature",
                   Antithetic = "Antithetic Monte Carlo",
                   Quasi = "Quasi Monte Carlo via the Sobol sequence")),
    sprintf("%.2f [%.2f, %.2f] seconds and total computation time %.2f [%.2f, %.2f], with median number of iterations %.0f",
            et[1], et[2], et[3], tot[1], tot[2], tot[3], median(X$iters))
  )
    
  list(tab = tab,
       cap = cap)
}

GHtab <- to.tab(GHparsed, "GH")
ATtab <- to.tab(ATparsed, "Antithetic")
QItab <- to.tab(QIparsed, "Quasi")

tab <- cbind(GHtab$tab, ATtab$tab[,-1], QItab$tab[,-1])
p <- tab$Parameter
p <- gsub("\\[", "_{", p)
p <- gsub("]$", "}", p)
p <- gsub("zeta_bin", "zeta", p)
p <- gsub("\\(Intercept\\)$", "0}", p)
p <- gsub("time$", "1}", p)
p <- gsub("cont$", "1}", p)
p <- gsub("bin$", "1}", p)
p <- gsub("Y\\.1\\_", "beta_{1", p)
p <- gsub("Y\\.2\\_", "beta_{2", p)
p <- gsub("Y\\.3\\_", "beta_{3", p)
p <- paste0("$\\", p, "=", .f(targets), "$")
tab$Parameter <- p

library(xtable)
xt <- xtable(tab, caption = paste0(GHtab$cap, ATtab$cap, QItab$cap))
print(xt, include.rownames = F,
      sanitize.text.function = identity)

# Make the usual plot -----------------------------------------------------
fits <- list(
  "Gauss-Hermite quadrature" = GHfits,
  "Antithetic Monte Carlo" = Antifits,
  "Quasi Monte Carlo" = Quasfits
)

fits2 <- do.call(rbind, lapply(seq_along(fits), function(i){
  df <- data.frame(method = names(fits)[i], t(sapply(fits[[i]], function(x) c(x$coeffs$gamma, x$coeffs$zeta))))
  tidyr::pivot_longer(df, cols = "gamma_1":"zeta_bin")
}))

library(dplyr)
library(ggplot2)

fits2 %>% 
  mutate(
    name = case_when(
      grepl("\\_bin$", name) ~ gsub("\\_bin$", "", name),
      T ~ paste0(gsub("\\_", '[', name), "]")
    ),
    target = case_when(
      name == "zeta" ~ -0.2,
      name == "gamma[2]" ~ -0.5,
      T ~ 0.5
    )
  ) %>% 
  mutate_at("method", forcats::fct_inorder) %>% 
  ggplot(aes(x = method, y = value, fill = method)) + 
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey3') + 
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50) + 
  facet_wrap(~name, scales = 'free_y', labeller = label_parsed) +
  theme_csda() + 
  labs(x = "E-step method", y = "Estimate") + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  theme(
    legend.position = 'none',
    strip.placement = 'outside',
    axis.text.x = element_text(size = 4.5)
  )

ggsave("/data/c0061461/THESIS/gmvjoint-MC/MCEM_gamzet.png",
       device = grDevices::png,
       type = 'cairo', width = 140, height = 90, units = 'mm', res = 250L)






