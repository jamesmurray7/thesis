rm(list=ls())
source(".Rprofile")
load.dir <- save.dir.file.path("fits")
N <- 100
fits.file.names <- dir(load.dir, pattern = "\\.RData") # Should already be in order
file.names.noext <- gsub("\\.RData$", "", fits.file.names)

loaded <- setNames(lapply(seq_along(fits.file.names), function(i){
  this <- fits.file.names[i]; this.nm <- file.names.noext[i]
  K <- as.numeric(gsub("^K", "", gsub("\\_.*$", "", this.nm)))
  method <- gsub("^K\\d\\_", "", this.nm)
  cat(sprintf("Loading saved fit object for K: %d; method: %s\n", K, method))
  if(method == "gmv"){
    assign("fits", get(load(save.dir.file.path(this, load.dir))))
    out <- lapply(fits, extract.from.joint)
  }else{
    assign("step", get(load(save.dir.file.path(gsub("jML","gmv", this), load.dir))))
    step <- extract.from.joint(step[[1]])
    assign("fits", get(load(save.dir.file.path(this, load.dir))))
    out <- lapply(fits, extract.from.jML, names(step$SE))
  }
  out
}), file.names.noext)

parsed <- lapply(loaded, function(x){
  this.fits <- x
  # This a list of length N, either full of joineRML or gmvjoint fits.
  # Regardless, they *should* have the same names 
  Omega <- sapply(this.fits, '[[', "Omega")
  SE <- sapply(this.fits, '[[', "SE")
  lb <- sapply(this.fits, '[[', "lb")
  ub <- sapply(this.fits, '[[', "ub")
  et <- unname(sapply(this.fits, '[[', "elapsed"))
  tt <- unname(sapply(this.fits, '[[', "total"))
  list(
    Omega = Omega, SE = SE, lb = lb, ub = ub, et = et, tt = tt
  )
})

Ktab.function <- function(k, method = c("gmv", "jML")){
  D.true.diag <-  rep(c(.25, .06), k)
  target.mat <- t(create.targetmat(K = k, 
                                   arguments = list(D = diag(D.true.diag, ncol = 2*k, nrow = 2*k))))
  x <- parsed[[paste0("K", k, "_", method)]]
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
  
  return(list(Emp.Mean = Emp.Mean, Emp.SD = Emp.SD,
              Avg.SE = Avg.SE, CP = CP, MSE = MSE, Bias = Bias,
              targets = target.mat[,1], method = switch(method,
                                                        jML = "joineRML",
                                                        gmv = "gmvjoint"), K = k))
  
}

Ktabs <- setNames(lapply(seq_along(parsed), function(x){
  this.nm <- names(parsed)[x]
  K <- as.numeric(gsub("^K", "", gsub("\\_.*$", "", this.nm)))
  method <- gsub("^K\\d\\_", "", this.nm)
  Ktab.function(K, method)
}),names(parsed))

.f <- function(x) format(round(x, 3), nsmall = 3, justify = "right")
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
                    row.names = NULL, method = this$method, stringsAsFactors = F)
  names(out)[2] <- "Mean (SD)"
  out
  
})


# Forget above; just plot stuff of interest? ------------------------------

# Elapsed time ------------------------------------------------------------
library(dplyr)
library(ggplot2)
source("../theme_csda.R")

nm <- names(parsed)
elapsed <- do.call(rbind, lapply(seq_along(parsed), function(i){
  this <- parsed[[i]]; this.nm <- nm[i]
  K <- as.numeric(gsub("^K", "", gsub("\\_.*$", "", this.nm)))
  method <- gsub("^K\\d\\_", "", this.nm)
  et <- this$et
  data.frame(method = ifelse(method == "gmv", "Approximate EM", "joineRML"), K = K, elapsed = et,
             row.names = NULL, stringsAsFactors = FALSE)
}))

elapsed %>% 
  ggplot(aes(x = as.factor(K), y = elapsed, fill = method)) + 
  geom_vline(xintercept = c(1.5,2.5), lty = 3, alpha = .25)+
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, position = position_dodge(width = 0.9)) + 
  scale_y_log10("Elapsed time (seconds)") + 
  facet_wrap(~method, scales = 'free_y') + 
  scale_x_discrete(expression("Number of longitudinal responses,"~K)) + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme_csda() + 
  theme(
    legend.position = 'none',
    strip.placement = "outside",
    strip.text = element_text(vjust=1)
  )

ggsave(save.dir.file.path("comparison_elapsed.png", load.dir),
       width = 140, height = 65, units = "mm")


# Parameter estimates -----------------------------------------------------
# 1. Survival parameters;
# 2. Fixed effects;
# 3. vech(D).

gam.zet <- do.call(rbind, lapply(seq_along(parsed), function(i){
  this <- parsed[[i]]; this.nm <- nm[i]
  K <- as.numeric(gsub("^K", "", gsub("\\_.*$", "", this.nm)))
  method <- gsub("^K\\d\\_", "", this.nm)
  est <- this$Omega
  targets <- t(create.targetmat(K = K))
  rn <- row.names(est); rnt <- row.names(targets)
  rn.to.keep <- which(grepl("gamma\\_1", rn)):length(rn)
  rnt.to.keep <- which(sapply(rn[rn.to.keep], grepl, rnt)[,1]):length(rnt)
  
  est <- as.data.frame(t(est[rn.to.keep,]))
  tar <- as.data.frame(t(targets[rnt.to.keep,]))
  inds <- 1:nrow(est)
  est <- cbind(inds, est); tar <- cbind(inds, tar)
  lb <- cbind(inds, as.data.frame(t(this$lb[rn.to.keep,])))
  ub <- cbind(inds, as.data.frame(t(this$ub[rn.to.keep,])))
  est <- est %>% tidyr::pivot_longer(-inds, names_to = "parameter", values_to = "estimate")
  lb <- lb %>% tidyr::pivot_longer(-inds, names_to = "parameter", values_to = "lower")
  ub <- ub %>% tidyr::pivot_longer(-inds, names_to = "parameter", values_to = "upper") 
  tar <- tar %>% tidyr::pivot_longer(-inds, names_to = "parameter", values_to = "target") %>% 
    mutate(parameter = ifelse(parameter == "\\zeta", "zeta_bin", gsub("\\\\", "", parameter)))
  out <- left_join(est, lb, by = c("inds", "parameter")) %>% 
    left_join(., ub, by = c("inds", "parameter")) %>% 
    left_join(., tar,by = c("inds", "parameter")) %>% 
    mutate(
      parameter = case_when(
        grepl("gamma", parameter) ~ paste0(gsub("_", '[', parameter),']'),
        parameter == "zeta_bin" ~ "zeta"
      ),
      method = ifelse(method == "gmv", "Approximate EM", "joineRML"),
      K = forcats::fct_inorder(paste0("K = ", K))
    )
  
  out
}))

gam.zet$parameter <- factor(gam.zet$parameter, levels = unique(sort(gam.zet$parameter)))

unq.targets <- gam.zet %>% distinct(parameter, target)
means <- gam.zet %>% 
  group_by(K, method, parameter) %>% 
  summarise(m = mean(estimate),
            l = mean(lower), u = mean(upper),
            .groups = "keep")

gam.zet2 <- left_join(gam.zet, means, by = c("K", "method", "parameter"))

gam.zet2 %>% 
  ggplot(aes(x = parameter, y = estimate, fill = forcats::fct_inorder(method))) + 
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, position = position_dodge(width = 1.1)) + 
  geom_point(aes(y = target), show.legend = FALSE, pch = 18, colour = "cyan4", size = 1) + 
  geom_point(aes(y = l, colour = forcats::fct_inorder(method)), position = position_dodge(width = 1.1), pch = 4,
             show.legend = FALSE, colour = rgb(0, 0, 0, .005), size = 0.66) + 
  geom_point(aes(y = u, colour = forcats::fct_inorder(method)), position = position_dodge(width = 1.1), pch = 4,
             show.legend = FALSE, colour = rgb(0, 0, 0, .005), size = 0.66) + 
  facet_wrap(~K, scales = "free_x") + 
  scale_x_discrete(labels = function(l) parse(text=l)) + 
  scale_fill_brewer(palette = "YlOrRd") + 
  labs(x = NULL, y = "Estimate", fill = NULL) + 
  theme_csda() + 
  theme(
    legend.position = "bottom"
  )

ggsave(save.dir.file.path("Kgammazeta.png", load.dir),
       width = 140, height = 90, units = "mm")

