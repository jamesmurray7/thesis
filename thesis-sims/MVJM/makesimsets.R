N <- 100
# n -----------------------------------------------------------------------
out <- save.dir.file.path('n/data')
filename <- function(x) save.dir.file.path(paste0("n", x, '.RData'), out)
choices <- c(100, 250, 500, 1000)
for(i in choices){
  cli::cli_alert_info("Creating data for n = {i}.")
  fn <- create.simfn(arguments = list(n = i, theta = c(-3, 0.1)))
  sim.sets <- createNsims(fn, N)
  save(sim.sets, file = filename(i))
  cli::cli_alert_success("Done, saved in:\n{filename(i)}.")
  DataSummary(sim.sets)
  cli::cli_rule()
  rm(sim.sets)
}
gc()