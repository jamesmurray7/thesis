N <- 100
# K -----------------------------------------------------------------------
out <- save.dir.file.path('data')
filename <- function(x) save.dir.file.path(paste0("K", x, '.RData'), out)
choices <- c(3,5,7)
for(i in choices){
  cli::cli_alert_info("Creating data for K = {i}.")
  Ddiag <- rep(c(.25, .06), i)
  fn <- create.simfn(arguments = list(theta = c(-2.9, 0.1), 
                                      n = 250, D = diag(Ddiag, ncol = 2*i, nrow = 2*i)), K = i)
  sim.sets <- createNsims(fn, N)
  save(sim.sets, file = filename(i))
  cli::cli_alert_success("Done, saved in:\n{filename(i)}.")
  DataSummary(sim.sets)
  cli::cli_rule()
  rm(sim.sets)
}
gc()


