rm(list=ls())
source(".Rprofile")
# file.show("printout.txt", title = "models for reference")
log.dir <- save.dir.file.path("Disps")
resp.dirs <- dir(log.dir)
# Establish directories containing .log files
dirs.to.parse <- sapply(resp.dirs, save.dir.file.path, log.dir)
dirs.to.parse <- dirs.to.parse[!grepl("\\.log|OLD", dirs.to.parse)]


# Create .log dumps -------------------------------------------------------
#           (summaries)
# And storing 
out <- setNames(vector("list", length(dirs.to.parse)), names(dirs.to.parse))
for(d in dirs.to.parse){
  d.files <- dir(d, pattern = "\\.log")
  # get response
  resp <- gsub(paste0(save.dir,"\\/Disps\\/"), "", d)
  if(length(d.files) == 0){
    cat(sprintf("Nothing present in %s currently!\n", d))
    next
  }
  d.collect <- structure(matrix(NA, 0, 5),
                         dimnames=list(c(),c("AIC", "BIC", "logLik", "deviance", "df.resid") ))
  n <- 0
  for(f in d.files){
    ff <- save.dir.file.path(f, d)
    log <- readLines(ff)
    ind <- grep("AIC", log) + 1
    fixed <- gsub("\\s+", " ", trimws(log[ind]))
    if(grepl("NA", fixed)){
      cat(sprintf("\nNA values in %s, skipping...\n", f))
      next
    }
    n <- n + 1
    fixed <- setNames(as.numeric(el(strsplit(fixed,"\\s"))),
                      c("AIC", "BIC", "logLik", "deviance", "df.resid"))
    d.collect <- rbind(d.collect, f=fixed)
    row.names(d.collect)[n] <- f
    rm(fixed)
  }
  
  cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by AIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), AIC));cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by BIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), BIC));cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by logLik"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), logLik));cat("\n")
  
  # And now sink 
  sink(file = save.dir.file.path(paste0(resp,".log"), log.dir))
  cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by AIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), AIC));cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by BIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), BIC));cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by logLik"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), logLik));cat("\n")
  sink()

  fs::file_copy(save.dir.file.path(paste0(resp,".log"), log.dir),
                paste0('./summarylogs/DISP_', resp, ".log"),
                TRUE)
  
  # "Rotate" around and store
  d.collect <- as.data.frame(d.collect)
  d.collect$name <- sub("^\\d.*?\\_","",row.names(d.collect))
  d.collect$name <- sub("\\.log$", "", d.collect$name)
  row.names(d.collect) <- NULL
  d.collect <- dplyr::arrange(d.collect, BIC)
  d.out <- tidyr::pivot_longer(data = d.collect[1:10,], AIC:df.resid, names_to = "criterion", values_to = "value")
  d.out$response <- resp
  out[[resp]] <- d.out
}

# Comparing to intercept only fits.
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)
prt <- glmmTMB(prothrombin ~ histologic2 + sex + time + (1 + time|id),
               data = PBCredx, family = Gamma(link="log"))
alk <- glmmTMB(alkaline ~ age + histologic2 + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
               data = PBCredx, family = glmmTMB::nbinom2())
pla <- glmmTMB(platelets ~ age + histologic2 + ns(time, knots = c(1, 4)) + (1 + ns(time, knots = c(1, 4))|id),
               data = PBCredx, family = glmmTMB::genpois())

# Manually insert these...
out$alkaline <- rbind(out$alkaline, data.frame(name = "1", criterion = "BIC", value = BIC(alk), response = "alkaline"))
out$platelets <- rbind(out$platelets, data.frame(name = "1", criterion = "BIC", value = BIC(pla), response = "platelets"))
out$prothrombin <- rbind(out$prothrombin, data.frame(name = "1", criterion = "BIC", value = BIC(prt), response = "prothrombin"))

out.all <- do.call(rbind, out)
  
# BIC only
bics <- out.all[out.all$criterion=="BIC", ] 
bics <- bics[!is.na(bics$value),] # no idea where these come from!
#  Tidy a bit
bics$response <- tools::toTitleCase(bics$response)
bics$name1 <- substr(toupper(bics$name),1,1)
bics$name1 <- factor(bics$name1, c("1", "A", "D", "H", "S", "T"))
# For drawing ring around min BIC
bics$min.bic <- ifelse(bics$value %in% tapply(bics$value, bics$response, min), bics$value, NA)

library(ggplot2)
ggplot(bics, 
       aes(x = name1, y = value)) + 
  # Draw ring around smallest
  geom_point(aes(y = min.bic), size = .75, shape = 0, colour = "black") + 
  geom_point(size=.5, colour = .nice.orange) + 
  scale_shape_manual(values = c(19, 4)) + 
  facet_wrap(~ response, scales = "free_y", nrow = 1) + 
  theme_csda() + 
  # scale_colour_brewer(palette = "YlOrRd") + 
  # Shift YlOrRd because hard to see!
  labs(x = "", y = "BIC") + 
  theme(
    axis.text.x = element_text(size=5.5),
    axis.text.y = element_text(size=4),
    strip.text = element_text(size=6,vjust=1)
  )

ggsave("./allDispmodels.png", width = 140, height = 45, units = "mm")
