rm(list=ls())
source('.Rprofile')
library(splines)
log.dir <- save.dir.file.path("Joint/Univs")

load(save.dir.file.path("AllUnivs.RData",log.dir))


# Plotting all survival sub-models ----------------------------------------
nj <- names(joints)
ests <- do.call(rbind, lapply(seq_along(joints), function(i){
  j <- nj[i]
  J <- joints[[i]]
  S <- as.data.frame(summary(J)$Surv)
  S$response <- j
  S$parameter <- row.names(S)
  row.names(S) <- NULL
  S
}))

ests$response <- ifelse(ests$response == "serBilir", "Serum bilirubin", tools::toTitleCase(ests$response))
ests$response <- ifelse(ests$response == "SGOT", "AST", ests$response)

library(dplyr)
library(ggplot2)
ests <- ests %>% 
  mutate(
    parameter2 = case_when(
      grepl("age", parameter) ~ 'zeta[~1]',
      grepl("sex", parameter) ~ "zeta[~2]",
      grepl("histologic2", parameter) ~ "zeta[~3]",
      T ~ "gamma"
    )
  )

ggplot(ests, aes(x = parameter2, y = Estimate)) + 
  geom_hline(aes(yintercept=0), lwd = .25, lty = 3) +
  geom_point(size=.5) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, lwd = .25) + 
  scale_x_discrete("Parameter", labels = function(l) parse(text=l)) +
  scale_y_continuous(breaks=scales::pretty_breaks(6)) + 
  facet_wrap(~response, scales = "free_y", nrow = 2) + 
  theme_csda() + 
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=4),
    axis.title = element_text(size=6),
    strip.text = element_text(size=6,vjust=1)
  )

ggsave("./UnivSurvModels.png", width = 140, height = 90, units = "mm")


# Elapsed times -----------------------------------------------------------
# Not sure why bothered with this ! 
# ets <- do.call(rbind, lapply(joints, function(x){
#   resp <- x$ModelInfo$Resps
#   time <- x$elapsed.time
#   time <- c(time[1] + time[2], time[3])
#   et <- as.data.frame(time)
#   et$measure <- row.names(et); row.names(et) <- NULL
#   et$response <- resp
#   et
# }))
# 

# Producing _all_ residual plots ------------------------------------------
# Make ggplots for better automatic scaling...
# Ensure outputted to screen in alphabetical order...
# Recall SGOT is 'AST'
alpha.names <- c("albumin", "alkaline", "SGOT", "hepatomegaly", "platelets", "prothrombin", "serBilir")
out <- setNames(vector("list", length(alpha.names)), alpha.names)
for(f in alpha.names){
  resp <- ifelse(f == "serBilir", "Serum bilirubin", tools::toTitleCase(f))
  resp <- ifelse(resp == "SGOT", "AST", resp)
  
  # Make ggplots for better automatic scaling...
  R <- residuals(joints[[f]], "l", "p")
  x <- attr(R[[1]], "fitted")
  out[[f]] <- data.frame(response = resp, fitted = x, pr = R[[1]], row.names = NULL)
}
out.all <- do.call(rbind, out)

out.all %>% 
  ggplot(aes(x = fitted, y = pr)) + 
  geom_hline(aes(yintercept = 0), lty = 5, col = 'red', lwd = 0.25) + 
  geom_point(size=.005, alpha = .2) + 
  facet_wrap(~ response, nrow = 2, scales = 'free') + 
  labs(y = "Pearson residuals", x = "Fitted") + 
  theme_csda()+ 
  theme(
    axis.text.x = element_text(size=3),
    axis.text.y = element_text(size=3),
    axis.title = element_text(size=4),
    strip.text = element_text(size=4,vjust=1)
  )

ggsave("UnivPearsonResiduals.png", width = 140, height = 60, units = "mm")


# Survival residuals ------------------------------------------------------

for(f in alpha.names){
  resp <- ifelse(f == "serBilir", "Serum bilirubin", tools::toTitleCase(f))
  resp <- ifelse(resp == "SGOT", "AST", resp)
  
  # Make ggplots for better automatic scaling...
  R <- residuals(joints[[f]], "s")
  plot(R)
  legend("topright", legend = resp, bty = "n")
}






