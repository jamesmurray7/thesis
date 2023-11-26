library(xtable)
load(save.dir.file.path("ROC_windows.RData"))

lapply(ROCs, function(R){
  cat(cli::rule())
  cat(sprintf("\nWindow: [%.1f, %.1f]\n", R$Tstart, R$Tstart + R$delta))
  cap <- paste0(sprintf("Window: [%.1f, %.1f]. ", R$Tstart, R$Tstart + R$delta),
                "`Acc': Accuracy. ",
                "Probabilistic thresholds c_j which are repeated are ommitted and the greater presented.")
  xt <- xtable(R$metrics, , digits = c(1,2,0,0,0,0,2,2,2,2,2,2,2)
               caption = cap)
  xtable:::print.xtable(xt, 
                        include.rownames = FALSE,
                        sanitize.text.function = identity,
                        sanitize.colnames.function = identity,
                        size = 'tiny')
  cat(cli::rule())
  cat("\n")
})
