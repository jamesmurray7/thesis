# Setup -------------------------------------------------------------------
rm(list=ls())
source('.Rprofile')
library(dplyr)
library(ggplot2)
log.dir <- save.dir.file.path("Joint/Multivs")
(to.load <- dir(log.dir, pattern = "\\.RData"))
to.load <- to.load[!to.load%in%c("OneBigModel.RData", "biv.RData", "triv.RData")]
load(save.dir.file.path("../Univs/AllUnivs.RData", log.dir)) # For comparing?

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

out.all %>% 
  filter(model %in% c("Blood clotting and flow", 
                      "Liver enzymes", "Liver health and function")) %>% 
  ggplot(aes(x = parameter2, y = Estimate)) + 
  geom_hline(aes(yintercept=0), lwd = .25, lty = 3) +
  geom_point(size=.5) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, lwd = .25) + 
  scale_x_discrete("Parameter", labels = function(l) parse(text=l)) +
  scale_y_continuous(breaks=scales::pretty_breaks(6)) + 
  facet_wrap(~model, scales = "free", nrow = 1) + 
  theme_csda() + 
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=4),
    axis.title = element_text(size=6),
    strip.text = element_text(size=6,vjust=1)
  )
  
ggsave("IntermediateMultivs.png",
       width = 140, height = 60, units = "mm")

# Printing `xtables` for appendix -----------------------------------------
require(xtable)

xtable(joint.bloods, 
       label = "tab:appendix-PBC-multivs-bloods",
       booktabs = FALSE, size = "footnotesize", max.row = 10)
xtable(joint.enzymes,
       label = "tab:appendix-PBC-multivs-enzymes",
       booktabs = FALSE, size = "footnotesize", max.row = 14)
xtable(joint.health, 
       label = "tab:appendix-PBC-multivs-health",
       booktabs = FALSE, size = "footnotesize", max.row = 7)

# Outputting cov/cor matrices ---------------------------------------------
D.R <- function(x){
  stopifnot(inherits(x, "joint"))
  D <- x$coeffs$D
  row.names(D) <- colnames(D) <- NULL
  R <- cov2cor(D)
  inds <- x$ModelInfo$inds$R$b
  rn <- lapply(seq_along(inds), function(i){
    ii <- 0:(length(inds[[i]]) - 1)
    paste0("$\\D_{", i, ",", ii, "}$")
  })
  D[upper.tri(D, F)] <- R[upper.tri(R, F)]
  D <- apply(D, 2, function(a) format(round(a, 4), nsmall = 4, justify = "right"))
  diag(D) <- paste0("\\cbb{",diag(D),"}") # Shading on diagonal
  row.names(D) <- colnames(D) <- unlist(rn)
  xt <- xtable(D)
  print(xt, sanitize.text.function = identity, sanitize.rownames.function = identity,
        sanitize.colnames.function = identity)
  invisible(D)
}

# Printing
D.R(joint.bloods)
D.R(joint.enzymes)
D.R(joint.health)


