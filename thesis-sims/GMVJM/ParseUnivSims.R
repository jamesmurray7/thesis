rm(list=ls())
source(".Rprofile")

in.dir <- save.dir.file.path("fits")
(files <- dir(in.dir, pattern = "\\.RData"))
univs <- files[which(grepl("Univ", files))]
# Hardcoding the three I want
univs <- univs[c(2,3,5)]
nms <- gsub("\\.RData$|^Univ|\\_Surgery", "", univs)
# Load list of fitted joint objects
.loader <- function(f){
  assign("out", get(load(save.dir.file.path(f, in.dir))))
  out
}

parsed <- setNames(lapply(univs, function(x){
  ff <- .loader(x)
  this.parsed <- try(lapply(ff, extract.from.joint), silent = T)
  out <- NULL
  if(!inherits(this.parsed, "try-error")){
    Omega <- sapply(this.parsed, function(y) y$Omega)
    SE <- sapply(this.parsed, function(y) y$SE)
    lb <- sapply(this.parsed, function(y) y$lb)
    ub <- sapply(this.parsed, function(y) y$ub)
    et <- unname(sapply(this.parsed, function(y) y$elapsed))
    tt <- unname(sapply(this.parsed, function(y) y$total))
    mo <- unname(sapply(this.parsed, function(y) y$model))
    out <- list(Omega = Omega, SE = SE, lb = lb, ub = ub, elapsed = et, total = tt,
                model = unique(mo))
  }
  out
}), nms)

# Figure (3x gamma + zetaa) -----------------------------------------------
library(dplyr)
library(ggplot2)

i <- 1
gam.zet.sig <- do.call(rbind, lapply(parsed, function(x){
  gamzetsig <- x$Omega[grepl("gamma|zeta|sigma|shape|phi", rownames(x$Omega)),]
  sig.shape.ind <- grep("sigma|shape|phi", rownames(gamzetsig))
  gam.ind <- grep("gamma", rownames(gamzetsig))
  zet.ind <- grep("zeta", rownames(gamzetsig))
  mm <- x$model
  if(is.null(mm)) return(NULL)
  # Univariate -> exploit this
  gam <- data.frame(Parameter = "gamma", estimate = unname(gamzetsig[gam.ind,,drop=T]),
                    target = 0.5,
                    stringsAsFactors = F, row.names = NULL)
  zet <- data.frame(Parameter = "zeta", estimate = unname(gamzetsig[zet.ind,,drop=T]),
                    target = -0.2,
                    stringsAsFactors = F, row.names = NULL)
  sig <- data.frame(Parameter = "sigma", estimate = unname(gamzetsig[sig.shape.ind,,drop=T]),
                    target = ifelse(mm == "Gamma", 2, ifelse(mm == "negbin", 1, -0.3)),
                    stringsAsFactors = F, row.names = NULL)
  out <- rbind(gam,zet,sig)
  out$i <- i; out$model <- mm
  rownames(out) <- NULL
  i <<- i + 1
  return(out)
}))

# Plot sigma separately!
gam.zet.sig %>% 
  filter(Parameter != "sigma") %>% 
  mutate(
    model = case_when(
      model == "genpois" ~ "Generalised Poisson",
      model == "negbin" ~ "Negative binomial",
      T ~ model
    )
  ) %>% 
  mutate(
    model.f = forcats::fct_inorder(model),
  ) %>% 
  ggplot(aes(x = model.f, y = estimate)) +
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey3') +
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, fill = .nice.orange) + 
  facet_wrap(~Parameter, scales = 'free', labeller = label_parsed, nr = 1) +
  theme_csda() + 
  labs(x = NULL, y = "Estimate") + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  theme(
    legend.position = 'none', 
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 5)
  )

ggsave(save.dir.file.path("gamzet.png", in.dir),
       width = 140, height = 60, units = "mm")

# sigma, with strip titles `hacked`.
gam.zet.sig %>% 
  filter(Parameter == "sigma") %>% 
  mutate(
    model = case_when(
      model == "genpois" ~ "Generalised Poisson",
      model == "negbin" ~ "Negative binomial",
      T ~ model
    )
  ) %>% 
  mutate(
    model.f = forcats::fct_inorder(model),
  ) %>% 
  ggplot(aes(x = model.f, y = estimate)) +
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey3') +
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, fill = .nice.orange) + 
  facet_wrap(~model, scales = 'free', nr=1, # Overwrite the panelling with \sigma
             labeller = function(l) parse(text="sigma"))+
  theme_csda() + 
  labs(x = NULL, y = "Estimate") + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  theme(
    legend.position = 'none', 
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    strip.text = element_text(size = 10)
  )

ggsave(save.dir.file.path("sigma.png", in.dir),
       width = 140, height = 60, units = "mm")

# Tabulation --------------------------------------------------------------
N=300L

# Pre-tabulating...
tabs <- setNames(lapply(seq_along(parsed), function(i){
  x <- parsed[[i]]
  if(is.null(x)) return(NULL)
  mm <- x$model
  if(mm != "genpois")
    target.mat <- t(create.targetmat(family = mm, N = N))
  else
    target.mat <- t(create.targetmat(family = mm, N = N,
                                     arguments = list(
                                       sigma = list(-0.3),
                                       random.formulas = list(~1), D = matrix(.30,1,1))
                                     ))
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
       targets = target.mat[,1], nm = nms[i], mod = mm)
  
}), nms)


# To xtable-ready object --------------------------------------------------
.f <- function(x) format(round(x, 3), nsmall = 3, justify = 'right')
to.xtab <- lapply(tabs, function(x){
  if(is.null(x)) return(NULL)
  this <- x
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

# Caption with usual stuff + computation times ----------------------------
caption <- lapply(parsed, function(x){
  if(is.null(x)) return(NULL)
  mm <- x$model
  fam <- switch(mm,
                Gamma = "Gamma",
                genpois = "Generalised Poisson",
                negbin = "Negative binomial")
  et <- x$elapsed; tot <- x$tot
  et <- et[et<800]; tot <- tot[tot<800]
  
  qn.et <- quantile(et, probs = c(.50, .25, .75)); qn.tot <- quantile(tot, probs = c(.50, .25, .75))
  list(
    caption = paste0("Parameter estimates for univariate ", fam, " joint models. Median [IQR] elapsed time taken for approximate EM algorithm to converge",
                     " and standard error calculation was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.et[1], qn.et[2], qn.et[3]),
                     " and total computation time was ", sprintf("%.3f [%.3f, %.3f] seconds.", qn.tot[1], qn.tot[2], qn.tot[3])),
    fam = fam
  )
  
})

# xtables -----------------------------------------------------------------
library(xtable)
xts <- lapply(seq_along(tabs), function(i){
  tab <- to.xtab[[i]]; cap <- caption[[i]]
  if(is.null(tab)) return(NULL)
  fam <- cap$fam; cap <- cap$caption
  align.lhs <- "cl|"
  align.rhs <- paste0(rep("r", ncol(tab)-1), collapse = '')
  align <- paste0(align.lhs, align.rhs)
  xt <- xtable(tab, caption = cap, align = align)
})

invisible(
  sapply(xts, print, sanitize.text.function = identity,
         booktabs = FALSE,
         include.rownames = FALSE, size = 'scriptsize')
)


# Elapsed time table ------------------------------------------------------
make.qn <- function(y){
  stopifnot(length(y)>1)
  qn <- quantile(y, probs = c(.5,.25,.75), names = F)
  qn <- .f(qn)
  paste0(qn[1], " [", qn[2], ", ", qn[3], "]")
}

Ga.el <- make.qn(parsed$Ga$elapsed); Ga.tt <- make.qn(parsed$Ga$total)
NB.el <- make.qn(parsed$NB$elapsed); NB.tt <- make.qn(parsed$NB$total)
GP.el <- make.qn(parsed$GP$elapsed); GP.tt <- make.qn(parsed$GP$total)

tab <- structure(rbind(
  cbind(Ga.el, Ga.tt),
  cbind(NB.el, NB.tt),
  cbind(GP.el, GP.tt)
), dimnames = list(c("Gamma", "Negative binomial", "Generalised Poisson"),
                   c("EM algorithm", "Total computation time")))
xtable(tab)
