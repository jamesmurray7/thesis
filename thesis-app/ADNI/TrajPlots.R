library(ggplot2)
plot.dir <- save.dir.file.path('plots')
filename <- function(x) save.dir.file.path(paste0("ADNILongProfile_", gsub("\\.", "_", x), ".png"), plot.dir)

# Plotting longitudinal trajectories --------------------------------------
# NOT trying to fit a good model -> just show y~x linear fit

nm.adni <- names(adni)
biomarkers <- nm.adni[which(nm.adni == "AV45"):which(nm.adni=="ICV")] # `block` in data
biomarkers <- biomarkers[!grepl("^Ecog", biomarkers) & biomarkers != "AV45" & biomarkers != "MOCA"] # Remove non-populated ones

plotfn <- function(b){ # b string containing name of biomarker.
  # Proportion MISSING (%)
  prop <- 100 * sum(is.na(adni[, b]))/nrow(adni)
  # Number (%) ids with at least one measurement present
  prop.id <- 100 * length(unique(adni[!is.na(adni[,b]), 'id']))/length(unique(adni$id))
  sset <- adni[, c("id", "time", "APOE4", "status", b, "survtime")]
  names(sset)[5] <- 'var'
  sset$A <- factor(ifelse(sset$APOE4 == 1, "Present", "Absent"), c("Present", "Absent"))
  sset$status <- ifelse(sset$survtime > 3, 0, sset$status)
  sset$S <- factor(ifelse(sset$status == 1, "Conversion to AD", "Censored"), c("Conversion to AD", "Censored"))
  
  p <- ggplot(data = sset[sset$time<=3,], aes(x = time, y = var, group = id)) + 
  geom_line(alpha = 0.10, lwd = 0.25) + 
  geom_smooth(aes(group = NULL, colour = A),
              method = 'lm', formula = y~x) + 
    geom_smooth(aes(group = NULL, colour = A),
                method = 'lm', formula = y~splines::ns(x, knots = c(0.5, 1.5)),
                lty = 3) + 
  facet_wrap(~S) +
  scale_color_brewer(palette = 'Dark2') + 
  labs(y = b, x = "Year", colour = expression("APOE "*epsilon[4]*":"),
       caption = sprintf("%.2f%% missing; %.2f%% of subjects have one measurement.",
                         prop, prop.id)) + 
  theme_csda() + 
    theme(plot.caption = element_text(size = 7),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  
  suppressWarnings(print(p))
  
  ggsave(filename(b), width = 148, height = 110, units = 'mm')
  
  cat(sprintf("Done %s.\n",b))
  invisible(b)
}

sapply(biomarkers, plotfn)

