# Setup -------------------------------------------------------------------
rm(list=ls())
source('.Rprofile')
library(dplyr)
library(ggplot2)
log.dir <- save.dir.file.path("Joint/Multivs")
(to.load <- dir(log.dir, pattern = "\\.RData"))
to.load <- to.load[to.load!="OneBigModel.RData"]

# Establish "clean names" ----
names.lookup <- data.frame(file.name = gsub("\\.RData$", "", to.load),
                           clean = c("Binary", "Blood clotting and flow",
                                     "Continuous", "Counts", "Liver enzymes",
                                     "Liver health and function"))
names.lookup$clean <- factor(names.lookup$clean, 
                             c("Continuous", "Counts", "Binary",
                               "Blood clotting and flow", "Liver enzymes",
                               "Liver health and function"))


resp.lookup <- data.frame(
  alpha.names = c("albumin", "alkaline", "SGOT", "hepatomegaly", "platelets", "prothrombin", "serBilir"),
  gamma.subscript = c("Alb", "Alk", "AST", "Hep", "Pla", "Pro", "Bil")
)


# Load individual GMVJMs, get ready to plot -------------------------------
out <- vector("list", nrow(names.lookup))
p <- 1
for(f in to.load){
  assign("X", get(load(save.dir.file.path(f, log.dir))))
  cat(f,"\n")
  print(X)
  f.model <- names.lookup[names.lookup$file.name == gsub("\\.RData$", "", f), "clean"]
  S <- as.data.frame(summary(X)$Surv)
  S$model <- f.model
  S$parameter <- row.names(S); row.names(S) <- NULL
  S$gs <- stringr::str_extract(S$parameter, "\\_.*$") %>% gsub("\\_","",.) %>%
    match(., resp.lookup$alpha.names) %>% resp.lookup[., 2]
  S <- S %>% 
    mutate(
      parameter2 = case_when(
        grepl("age", parameter) ~ 'zeta[~1]',
        grepl("sex", parameter) ~ "zeta[~2]",
        grepl("histologic2", parameter) ~ "zeta[~3]",
        T ~ paste0("gamma[~", gs, ']')
      )
    )
  out[[p]] <- S; p <- p + 1
}

out.all <- do.call(rbind, out)

out.all %>% 
  ggplot(aes(x = parameter2, y = Estimate)) + 
  geom_hline(aes(yintercept=0), lwd = .25, lty = 3) +
  geom_point(size=.5) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, lwd = .25) + 
  scale_x_discrete("Parameter", labels = function(l) parse(text=l)) +
  scale_y_continuous(breaks=scales::pretty_breaks(6)) + 
  facet_wrap(~model, scales = "free", nrow = 2) + 
  theme_csda() + 
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=4),
    axis.title = element_text(size=6),
    strip.text = element_text(size=6,vjust=1)
  )








