rm(list=ls())
source(".Rprofile")

in.dir <- save.dir.file.path("fits")
(files <- dir(in.dir))
triv <- files[which(grepl("Triv", files))]

# Load list of fitted joint objects
load(save.dir.file.path(triv, in.dir))

parsed <- lapply(fits, extract.from.joint)

# Figure (3x gamma + zetaa) -----------------------------------------------
library(dplyr)
library(ggplot2)


i <- 1
gam.zet <- do.call(rbind, lapply(parsed, function(x){
  gamzet <- x$Omega[grepl("gamma|zeta", names(x$Omega))]
  mm <- x$model
  out <- data.frame(Parameter = names(gamzet), estimate = gamzet, 
                    model = mm, i = i,
                    row.names = NULL)
  i <<- i + 1
  return(out)
}))

gam.zet %>% 
  mutate(
    paramlab = case_when(
      grepl("gamma", Parameter) ~ paste0("gamma[", stringr::str_extract(Parameter, "\\d$"),"]"),
      T ~ "zeta"
    ),
    target = rep(c(.5,-.5,.5,-.2), 300) # Fixed
  ) %>% 
  mutate(
    paramlab = forcats::fct_inorder(paramlab)
  ) %>% 
  ggplot(aes(x = paramlab, y = estimate)) +
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey3') +
  geom_boxplot(outlier.alpha = .33, lwd = 0.25, fatten = 2,
               outlier.size = .50, fill = .nice.orange) + 
  facet_wrap(~paramlab, scales = 'free', labeller = label_parsed, nr = 1) +
  theme_csda() + 
  labs(x = NULL, y = "Estimate") + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  theme(
    legend.position = 'none', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10)
  )

ggsave(save.dir.file.path("Triv_surv.png", in.dir),
       width = 140, height = 60, units = "mm")


# Tabulation --------------------------------------------------------------
N=300L
target.mat <- t(create.targetmat(N=N))

Omega <- sapply(parsed, '[[', "Omega")
SE <- sapply(parsed, '[[', "SE")
lb <- sapply(parsed, '[[', "lb")
ub <- sapply(parsed, '[[', "ub")
et <- unname(sapply(parsed, '[[', "elapsed"))
tt <- unname(sapply(parsed, '[[', "total"))

rn <- row.names(Omega); rn2 <- row.names(target.mat)
rn;rn2

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

tab.components <- list(
  Emp.Mean = Emp.Mean, Emp.SD = Emp.SD,
  Avg.SE = Avg.SE, CP = CP, MSE = MSE, Bias = Bias,
  targets = target.mat[,1]
)


# To xtable-ready object --------------------------------------------------
.f <- function(x) format(round(x, 3), nsmall = 3, justify = 'right')

this <- tab.components
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

# Caption -----------------------------------------------------------------
qn.et <- quantile(et, probs = c(.50, .25, .75)); qn.tot <- quantile(tt, probs = c(.50, .25, .75))
paste0("Median [IQR] elapsed time taken for approximate EM algorithm to converge",
       " and standard error calculation for was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.et[1], qn.et[2], qn.et[3]),
       " and total computation time was ", sprintf("%.3f [%.3f, %.3f] seconds", qn.tot[1], qn.tot[2], qn.tot[3])) -> caption


# xtable-out --------------------------------------------------------------
to.xtab <- out
align.lhs <- "cl|"
align.rhs <- paste0(rep("r", ncol(to.xtab)-1), collapse='')
align <- paste0(align.lhs, align.rhs)
library(xtable)
xt <- xtable(to.xtab, caption = caption,
             align = align)

print(xt, sanitize.text.function = identity,
      booktabs = FALSE,
      include.rownames = FALSE, size = 'scriptsize')
