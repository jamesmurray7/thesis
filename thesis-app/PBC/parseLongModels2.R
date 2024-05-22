source(".Rprofile")
log.dir <- save.dir.file.path("Longs")
resp.dirs <- dir(log.dir)
# Establish directories containing .log files
dirs.to.parse <- sapply(resp.dirs, save.dir.file.path, log.dir)
dirs.to.parse <- dirs.to.parse[!grepl("\\.log|OLD", dirs.to.parse)]


# Loop through and populate list of responses -----------------------------
out <- setNames(vector("list", length(dirs.to.parse)), names(dirs.to.parse))
for(d in dirs.to.parse){
  d.files <- dir(d, pattern = "\\.log")
  # get response
  resp <- gsub(paste0(save.dir,"\\/Longs\\/"), "", d)
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
      #cat(sprintf("\nNA values in %s, skipping...", f))
      next
    }
    n <- n + 1
    fixed <- setNames(as.numeric(el(strsplit(fixed,"\\s"))),
                      c("AIC", "BIC", "logLik", "deviance", "df.resid"))
    d.collect <- rbind(d.collect, f=fixed)
    row.names(d.collect)[n] <- f
    rm(fixed)
  }
  
  # Same printout as in parseLongModels(1).R
  cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by AIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), AIC)[1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by BIC"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), BIC)[1:7,]);cat("\n")
  cat(cli::rule(left = paste0(d, " -- Ranked by logLik"), background_col = "black",
                line = 2, col = .nice.orange))
  cat("\n")
  print(dplyr::arrange(as.data.frame(d.collect), logLik)[1:7,]);cat("\n")
  
  # "Rotate" around and store
  d.collect <- as.data.frame(d.collect)
  d.collect$name <- sub("^\\d.*?\\_","",row.names(d.collect))
  d.collect$name <- sub("\\.log$", "", d.collect$name)
  row.names(d.collect) <- NULL
  d.collect <- dplyr::arrange(d.collect, BIC)
  out[[resp]] <- tidyr::pivot_longer(data = d.collect[1:10,], AIC:df.resid, names_to = "criterion", values_to = "value")
  
}

# Repair names / stuff for plotting
SAH <- function(x) paste(toupper(substr(sort(el(strsplit(x, "_"))),1,1)),collapse="+")
out.all <- do.call(rbind, lapply(seq_along(out), function(i){
  resp.i <- names(out)[i]
  x <- out[[i]]
  # Time specification
  x$time.spec <- gsub("\\_", "", stringr::str_extract(x$name, "^.*?\\_"))
  # Drug interaction present?
  x$druginteract <- grepl("DI", x$name)
  # Covariates used
  x$covs <- unname(sapply(sub("DI\\_","",sub("^.*?\\_","",x$name)), SAH))
  
  data.frame(response = resp.i, time = x$time.spec, 
             drugInt = x$druginteract, covariates = x$covs, 
             criterion = x$criterion, value = x$value)
}))

library(ggplot2)
out.all$response <- ifelse(out.all$response == "serBilir", "Serum bilirubin", tools::toTitleCase(out.all$response))
out.all$response <- ifelse(out.all$response == "SGOT", "AST", out.all$response)

out.all$time <- sapply(out.all$time, switch,
                       int = "Linear + intercept only",
                       linear = "Linear",
                       ns = "Natural cubic splines",
                       quad = "Quadratic")
out.all$time <- factor(out.all$time, 
                       c("Linear + intercept only", "Linear", "Quadratic", "Natural cubic splines"))

# BIC only
bics <- out.all[out.all$criterion=="BIC", ] 
# For drawing ring around min BIC
bics$min.bic <- ifelse(bics$value %in% tapply(bics$value, bics$response, min), bics$value, NA)


ggplot(bics[bics$response!="Spiders",], 
       aes(x = covariates, y = value, colour = as.factor(time),
           shape = drugInt)) + 
  # Draw ring around smallest
  geom_point(aes(y = min.bic), size = .75, shape = 0, colour = "black") + 
  geom_point(size=.5) + 
  scale_shape_manual(values = c(19, 4)) + 
  facet_wrap(~ response, scales = "free", nrow = 2) + 
  theme_csda() + 
  guides(shape = F) + 
  # scale_colour_brewer(palette = "YlOrRd") + 
  # Shift YlOrRd because hard to see!
  # scale_colour_manual(values = c("#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")) + 
  scale_color_brewer(palette = 'Dark2') +
  labs(x = "", y = "BIC", colour = "Time specification:") + 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size=6),
    legend.text = element_text(size=5.5),
    axis.text.x = element_text(size=5.5),
    axis.text.y = element_text(size=4),
    strip.text = element_text(size=6,vjust=1)
  )

ggsave("./allLongmodels-corrections.png", width = 140, height = 90, units = "mm")
