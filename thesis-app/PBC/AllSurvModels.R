# All survival models -----------------------------------------------------
rm(list=ls())
qz <- qnorm(.975)
source('.Rprofile')
library(survival)
log.dir <- save.dir.file.path("Survs2")
head(surv.data)
pf <- parent.frame()

# Split histologic as 1-2 vs 3-4
surv.data$histologic2 <- ifelse(surv.data$histologic %in% c("3", "4"), 1, 0)
surv.data$histologic3 <- ifelse(surv.data$histologic %in% c("3", "4"), "3-4", as.character(surv.data$histologic))
surv.data$histologic3 <- factor(surv.data$histologic3, c("1", "2", "3-4"))

# file.name <- function(x) paste0("Surv_", gsub("\\~|\\+|\\s","",x), ".log")
make.and.fit <- function(svars = c("drug", "sex", "age", "histologic2")){
  svars <- unname(sapply(svars, match.arg, svars)) # Checking for typos
  cc <- length(svars)
  form.rhs <- paste0("~ ", paste0(svars, collapse = " + "))
  # fn <- paste0(cc, file.name(form.rhs))
  cat(sprintf("\nFitting:\n%s\n\n", form.rhs))
  form <- as.formula(paste0("Surv(survtime, status)", form.rhs), env = pf)
  ph <- coxph(form, data = surv.data, x = T)
  aic <- AIC(ph)
  conc <- summary(ph)$conc
  C <- conc[1]; C.lb <- conc[1] - qz * conc[2]; C.ub <- conc[1] + qz * conc[2]
  name <- toupper(substr(sort(svars),1,1))
  list(name = paste0(name, collapse = "+"), aic = aic,
       C = C, C.lb = C.lb, C.ub = C.ub)
}
# 4-variates
var4 <- make.and.fit()
# 3-variates
var3 <- combn(c("drug", "sex", "age", "histologic2"), 3)
var3s <- apply(var3,2,make.and.fit)
# 2-variates
var2 <- combn(c("drug", "sex", "age", "histologic2"), 2)
var2s <- apply(var2,2,make.and.fit)
# univariates
univs <- lapply(c("drug", "sex", "age", "histologic2"),make.and.fit)

# Bring all together
S <- c(list(var4), var3s, var2s, univs)

S.df <- do.call(rbind, lapply(S, function(x){
  data.frame(name = x$name, AIC = x$aic, 
             C = x$C, C.lb = x$C.lb, C.ub = x$C.ub)
}))

S.df$name <- forcats::fct_inorder(S.df$name)

# Add labels to 3 smallest AICs
S.df$AIC.lab <- ifelse(S.df$AIC < 1410, as.character(round(S.df$AIC, 1)), "")

library(ggplot2)
P1 <- ggplot(S.df, aes(x = name, y = AIC, label = AIC.lab)) + 
  geom_point(size=.66) +
  geom_text(nudge_y = 1.6, size = 1) + 
  labs(y = "AIC", x = "") + 
  theme_csda() + 
  theme(
    axis.text.x = element_text(angle = 360-45, vjust = -1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.title = element_text(size = 5.5),
  )
  

P2 <- ggplot(S.df, aes(x = name, y = C)) + 
  geom_point(size=.66) +
  geom_errorbar(aes(ymin = C.lb, ymax = C.ub), width = .15, lwd = .25) + 
  labs(y = "Harrell's C-index", x = "") + 
  theme_csda()  + 
  theme(
    axis.text.x = element_text(angle = 360-45, vjust = -1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.title = element_text(size = 5.5),
  )


ggpubr::ggarrange(P1, P2)
ggsave("./allsurvmodels.png", width = 140, height = 60, units = "mm")
