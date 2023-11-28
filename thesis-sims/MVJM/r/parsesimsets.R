N <- 100
# r -----------------------------------------------------------------------
load.dir <- save.dir.file.path('r/fits')
fits.file.names <- dir(load.dir, pattern = '\\.RData')
# Change order slightly
fits.file.names <- fits.file.names[c(3,4,1,2)]  # 3, 5, 10, 15
.loader <- function(filename){
  out <- assign('data', get(load(save.dir.file.path(filename, load.dir))))
  out
}

# Get the estimates, SEs, 2.5% and 97.5% estimates
parsed <- setNames(lapply(seq_along(fits.file.names), function(i){
  this.fits <- .loader(fits.file.names[i])
  # This is a list of length N, each containing the four parsed sets of data: 
  this.parsed <- lapply(this.fits, extract.from.joint)
  Omega <- sapply(this.parsed, function(y) y$Omega)
  SE <- sapply(this.parsed, function(y) y$SE)
  lb <- sapply(this.parsed, function(y) y$lb)
  ub <- sapply(this.parsed, function(y) y$ub)
  et <- unname(sapply(this.parsed, function(y) y$elapsed))
  tt <- unname(sapply(this.parsed, function(y) y$total))
  it <- unname(sapply(this.parsed, function(y) y$iter))
  list(Omega = Omega, SE = SE, lb = lb, ub = ub, elapsed = et, total = tt, iter = it)
}), gsub(".RData", '', fits.file.names))

# Get the target matrix:
target.mat <- t(create.targetmat())

# Get the table components print out.
tab.components <- lapply(parsed, function(x){
  Omega <- x$Omega; SE <- x$SE; lb <- x$lb; ub <- x$ub
  rn <- row.names(Omega); rn2 <- row.names(target.mat)
  
  # Work out what components in D to keep
  Dnames.tmat <- gsub("\\\\", '', rn2)
  Dnames.tmat <- Dnames.tmat[grepl('^D', Dnames.tmat)]
  rn.to.keep <- c(match(Dnames.tmat, rn), (match("Y.1_(Intercept)", rn):length(rn)))
  
  Omega <- Omega[rn.to.keep,]; SE <- SE[rn.to.keep,]; lb <- lb[rn.to.keep,]; ub <- ub[rn.to.keep,]
  row.names(Omega) <- row.names(SE) <- row.names(lb) <- row.names(ub) <-  rn2
  
  # Mean (SD) / (average SE)
  Emp.Mean <- rowMeans(Omega)
  Emp.SD <- apply(Omega, 1, sd)
  Avg.SE <- rowMeans(SE)
  
  # CP
  CP <- rowSums(lb < target.mat & ub > target.mat) / N
  
  # MSE and Bias
  MSE <- rowMeans((target.mat - Omega)^2)
  Bias <- rowMeans(Omega-target.mat)
  
  list(Emp.Mean = Emp.Mean, Emp.SD = Emp.SD,
       Avg.SE = Avg.SE, CP = CP, MSE = MSE, Bias = Bias,
       targets = target.mat[,1])
})

# And create the thing we want to pass to xtable.
.f <- function(x) format(round(x, 3), nsmall = 3, justify = 'right')
nm <- paste0("$", gsub("^r", "r=", names(parsed)), "$")
to.xtab <- do.call(cbind, lapply(seq_along(tab.components), function(i){
  this <- tab.components[[i]]
  param <- paste0("$", names(this$Emp.Mean), ' = ', trimws(.f(this$targets)), "$")
  param <- gsub("\\]", "}", gsub("D\\[","D_{", param))
  
  MSD <- paste0(.f(this$Emp.Mean), ' (', .f(this$Emp.SD), ')')
  SE <- .f(this$Avg.SE)
  Bias <- .f(this$Bias)
  MSE <- .f(this$MSE)
  CP <- format(round(this$CP, 2), nsmall = 2)
  
  out <- data.frame(Parameter = param, "MSD" = MSD, `SE` = SE, Bias = Bias, MSE = MSE, CP = CP,
                    row.names = NULL)
  names(out)[2] <- "Mean (SD)"
  if(i > 1) out$Parameter <- NULL
  out
}))

# Caption with usual stuff + computation times.
caption <- paste(sapply(seq_along(parsed), function(x){
  y <- parsed[[x]]
  et <- y$elapsed; tot <- y$tot
  et <- et[et<800]; tot <- tot[tot<800]
  qn.et <- quantile(et, probs = c(.50, .25, .75)); qn.tot <- quantile(tot, probs = c(.50, .25, .75))
  paste0("Median [IQR] elapsed time taken for approximate EM algorithm to converge",
         " and standard error calculation for ", nm[x], " was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.et[1], qn.et[2], qn.et[3]),
         " and total computation time was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.tot[1], qn.tot[2], qn.tot[3]))
}), collapse = '; ')
caption <- paste("Parameter estimates for differing sample sizes $\\r\\in\\lbr3,5,10,15\\rbr$. `Mean (SD)' denotes",
                 "the average estimated value with standard deviation of the parameter estimate. `SE' denotes the",
                 "mean standard error calculated at each model fit. Coverage probabilities are calculated from",
                 "$\\hbO\\pm1.96\\mathrm{SE}(\\hbO)$", caption, '\\!.')
align.lhs <- "cl|"
align.rhs <- paste0(rep("r", ncol(to.xtab)-1), collapse='')
align <- paste0(align.lhs, align.rhs)
library(xtable)
xt <- xtable(to.xtab, caption = caption,
             align = align)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- paste0("&", paste0("\\multicolumn{5}{c}{", nm, "}", collapse = '&'), '\\\\')

print(xt, sanitize.text.function = identity,
      add.to.row = addtorow,
      booktabs = FALSE,
      include.rownames = FALSE, size = 'scriptsize')


# Figure for survival parameters and comp.time ----------------------------
library(dplyr)
library(ggplot2)
gam.zet <- do.call(rbind,lapply(seq_along(parsed), function(x){
  ests <- parsed[[x]]$Omega; et <- parsed[[x]]$elapsed
  rn <- which(grepl("gamma", row.names(ests)))
  rn <- c(rn, max(rn)+1)
  df <- as.data.frame(t(ests[rn,]))
  df %>% tidyr::pivot_longer(everything()) %>% 
    mutate(r = readr::parse_number(nm[x]), param = case_when(
      grepl("gamma", name) ~ paste0(gsub("_", '[~', name),']'),
      name == "zeta_bin" ~ 'zeta'
    ),
    target = case_when(
      name == 'gamma_1' | name == "gamma_3" ~ 0.5,
      name == 'gamma_2' ~ -0.5,
      name == 'zeta_bin' ~ -0.2
    ))
}))

gam.zet %>% 
  mutate(rlab = forcats::fct_inorder(gsub("n=","",r))) %>% 
  ggplot(aes(x = rlab, y = value, fill = rlab)) + 
  geom_hline(aes(yintercept = target), lty=5, colour='grey3')+
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50) + 
  facet_wrap(~param, scales = 'free_y', labeller = label_parsed) +
  theme_csda() + 
  labs(x = expression("Length of follow-up period, "*r), y = "Estimate") + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  theme(
    legend.position = 'none',
    strip.text.x=element_text(margin=margin(b=3,t=5)),
    # strip.text = element_text(family="Times New Roman"),
    strip.placement = 'outside'
  )

ggsave(save.dir.file.path("r_gammazeta.png", load.dir),
       width = 148, height = 110, units = 'mm')


# Bias plot ---------------------------------------------------------------
all.bias <- do.call(rbind, lapply(seq_along(parsed), function(y){
  x <- parsed[[y]]
  Omega <- x$Omega; SE <- x$SE; lb <- x$lb; ub <- x$ub
  rn <- row.names(Omega); rn2 <- row.names(target.mat)
  
  # Work out what components in D to keep
  Dnames.tmat <- gsub("\\\\", '', rn2)
  Dnames.tmat <- Dnames.tmat[grepl('^D', Dnames.tmat)]
  rn.to.keep <- c(match(Dnames.tmat, rn), (match("Y.1_(Intercept)", rn):length(rn)))
  
  Omega <- Omega[rn.to.keep,]; SE <- SE[rn.to.keep,]; lb <- lb[rn.to.keep,]; ub <- ub[rn.to.keep,]
  row.names(Omega) <- row.names(SE) <- row.names(lb) <- row.names(ub) <-  rn2
  
  Bias <- Omega-target.mat
  Bias <- as.data.frame(t(Bias))
  
  Bias %>%
    tidyr::pivot_longer(everything()) %>% mutate(r = nm[y])
  
}))

step <- all.bias %>% 
  mutate(
    grp = case_when(
      grepl("^\\\\D", name) ~ "vech(D)",
      grepl("^\\\\beta|\\\\sigma", name) ~ "Longitudinal",
      T ~ "Survival"
    ),
    r = stringr::str_remove_all(r, "\\$"),
    r = forcats::fct_inorder(r)
  ) %>% 
  mutate(grp = forcats::fct_inorder(grp)) %>% 
  mutate(
    name = stringr::str_remove_all(name, "\\\\"),
    name2 = case_when(
      grepl("^D", name) ~ stringr::str_remove_all(name, "\\,"),
      grepl("^beta", name) ~ stringr::str_replace(stringr::str_replace(name, "\\}", "]"),
                                                  "\\_\\{", "["),
      grepl("sigma", name) ~ paste0("sigma[", stringr::str_extract(name, "\\d$"),"]^2"),
      grepl("gamma", name) ~ paste0("gamma[", stringr::str_extract(name, "\\d$"),"]"),
      T ~ name
    )
  )

plots <- setNames(vector("list", 3L), levels(step$grp))

for(j in levels(step$grp)){
  sset <- step %>% filter(grp == j) %>% mutate(name2 = forcats::fct_inorder(name2))
  sset$name2 <- forcats::fct_rev(sset$name2)
  sset.params <- sapply(levels(sset$name2), function(i) parse(text = as.character(i)))
  plots[[j]] <- sset %>% 
    ggplot(aes(y = name2, x = sqrt(abs(value)), fill = r)) + 
    geom_hline(aes(yintercept = 0), lty = 1, col = 'grey50') + 
    geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
                 outlier.size = .50) + 
    scale_y_discrete(labels = sset.params) + 
    scale_fill_brewer(palette = 'YlOrRd') + 
    facet_wrap(~grp, scales = 'free', strip.position = "right") +
    labs(y = '', x = expression(sqrt("|Bias|")), fill = "") + 
    # coord_flip() + 
    theme_csda() + 
    theme(
      # strip.text = element_blank()
      legend.position = ifelse(j == "Survival", "bottom", "none"),
      strip.placement = 'outside'
    )
}

png(file = save.dir.file.path("r_biases.png", load.dir),
    units = "mm", width = 148, height = 160, res = 300)
gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], 
                        layout_matrix = rbind(
                          c(1,1,2,2),
                          c(1,1,2,2),
                          c(3,3,2,2),
                          c(3,3,2,2)
                        ))
dev.off()
