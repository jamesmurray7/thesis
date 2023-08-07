N <- 100
# n -----------------------------------------------------------------------
load.dir <- save.dir.file.path('K/fits')
fits.file.names <- dir(load.dir, pattern = '\\.RData')
# Change order slightly
# fits.file.names <- fits.file.names[c(1,3,4,5,2)]  # {1,2,3,5,10}
Ks <- as.numeric(stringr::str_extract(fits.file.names, '\\d?\\d'))
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
  iter <- unname(sapply(this.parsed, function(y) y$iter))
  list(Omega = Omega, SE = SE, lb = lb, ub = ub, elapsed = et, total = tt,
       iter = iter)
}), gsub(".RData", '', fits.file.names))

# for each k in K
Ktab.function <- function(k){
  ind <- paste0("K", k)
  D.true.diag <-  rep(c(.25, .06), k)
  target.mat <- t(create.targetmat(K = k, 
                                   arguments = list(D = diag(D.true.diag, ncol = 2*k, nrow = 2*k))))
  
  x <- parsed[[ind]]
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
}

Ktabs <- setNames(lapply(Ks, Ktab.function), paste0("K", Ks))

# And create the thing we want to pass to xtable.
.f <- function(x) format(round(x, 3), nsmall = 3, justify = 'right')
to.xtab <- lapply(seq_along(Ktabs), function(i){
  this <- Ktabs[[i]]
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
  # if(i > 1) out$Parameter <- NULL
  out
})

# Caption with usual stuff + computation times.
nm <-  paste0("$K=", Ks,"$")
caption <- lapply(seq_along(parsed), function(x){
  y <- parsed[[x]]
  et <- y$elapsed; tot <- y$tot
  et <- et[et<800]; tot <- tot[tot<800]
  
  qn.et <- quantile(et, probs = c(.50, .25, .75)); qn.tot <- quantile(tot, probs = c(.50, .25, .75))
  paste0("Parameter estimates for ", nm[x], ". Median [IQR] elapsed time taken for approximate EM algorithm to converge",
         " and standard error calculation was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.et[1], qn.et[2], qn.et[3]),
         " and total computation time was ", sprintf("%.3f [%.3f, %.3f] seconds.", qn.tot[1], qn.tot[2], qn.tot[3]))
})

library(xtable)
xt <- lapply(seq_along(to.xtab), function(x){
  align.lhs <- "cl|"
  align.rhs <- paste0(rep("r", ncol(to.xtab[[x]])-1), collapse = '')
  align <- paste0(align.lhs, align.rhs)
  this <- to.xtab[[x]]
  xt <- xtable(to.xtab[[x]], caption = caption[[x]], align = align)
})

for(i in seq_along(xt)){
  xx <- readline(sprintf("Press Enter for K=%d printout", Ks[i]))
  print(xt[[i]], sanitize.text.function = identity,
                              booktabs = FALSE,
                              include.rownames = FALSE, size = 'scriptsize')
}

# Figure of elapsed time progression --------------------------------------
library(dplyr)
library(ggplot2)
times <- do.call(rbind, lapply(seq_along(parsed), function(i){
  Kn <- paste0(Ks[i])
  et <- parsed[[i]]$elapsed; tot <- parsed[[i]]$tot
  data.frame(K = Kn, elapsed = et, total = tot)
}))

times <- times %>% 
  tidyr::pivot_longer(elapsed:total) %>% 
  mutate(name = ifelse(name == "elapsed", "EM + Standard error", "Total computation"))

medians <- times %>% group_by(K, name) %>% summarise(med = median(value))

times %>% 
  ggplot(aes(x = K, y = value, fill = K)) + 
  geom_vline(xintercept = c(3.5,4.5), lty = 3, alpha = .25)+
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, position = position_dodge(width = 0.9)) + 
  # geom_line(data = medians, aes(x = K, y = med, group = name), lty = 'solid', colour = 'black') + 
  scale_y_log10("Elapsed time (seconds)") + 
  facet_wrap(~name, scales = 'free_y') + 
  scale_x_discrete(expression("Number of longitudinal responses,"~K)) + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme_csda() + 
  theme(legend.position = 'none')

ggsave(save.dir.file.path("K_times.png", load.dir),
       width = 148, height = 110, units = 'mm')


# iterations per second ---------------------------------------------------
times <- do.call(rbind, lapply(seq_along(parsed), function(i){
  Kn <- paste0(Ks[i])
  et <- parsed[[i]]$iter / parsed[[i]]$elapsed
  data.frame(K = Kn, iter_ps = et)
}))

times %>% 
  ggplot(aes(x = K, y = iter_ps, fill = K)) + 
  geom_vline(xintercept = c(3.5,4.5), lty = 3, alpha = .25)+
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, position = position_dodge(width = 0.9)) + 
  # geom_line(data = medians, aes(x = K, y = med, group = name), lty = 'solid', colour = 'black') + 
  scale_y_log10("Iterations per second") +
  scale_x_discrete(expression("Number of longitudinal responses,"~K)) + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme_csda() + 
  theme(legend.position = 'none')

ggsave(save.dir.file.path("K_iter_per_second.png", load.dir),
       width = 148, height = 90, units = 'mm')
