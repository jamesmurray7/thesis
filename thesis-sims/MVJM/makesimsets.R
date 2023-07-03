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


# K -----------------------------------------------------------------------
out <- save.dir.file.path('K/data')
filename <- function(x) save.dir.file.path(paste0("K", x, '.RData'), out)
choices <- c(1,2,3,5,7)
for(i in choices){
  cli::cli_alert_info("Creating data for K = {i}.")
  Ddiag <- rep(c(.25, .06), i)
  fn <- create.simfn(arguments = list(theta = c(-2.7, 0.1), 
                                      n = 500, D = diag(Ddiag, ncol = 2*i, nrow = 2*i)), K = i)
  sim.sets <- createNsims(fn, N)
  save(sim.sets, file = filename(i))
  cli::cli_alert_success("Done, saved in:\n{filename(i)}.")
  DataSummary(sim.sets)
  cli::cli_rule()
  rm(sim.sets)
}
gc()

# r -----------------------------------------------------------------------
out <- save.dir.file.path("r/data")
filename <- function(x) save.dir.file.path(paste0("r", x ,".RData"), out)
choices <- c(3,5,10,15)
for(i in choices){
  cli::cli_alert_info("Creating data for r = {i}.")
  fn <- create.simfn(arguments = list(n = 500, ntms = i,
                                      theta = c(-2.7, 0.0)))
  sim.sets <- createNsims(fn, N)
  save(sim.sets, file = filename(i))
  cli::cli_alert_success("Done, saved in:\n{filename(i)}.")
  DataSummary(sim.sets)
  cli::rule(col = "tomato2")
  rm(sim.sets)  
}


# omega -------------------------------------------------------------------
out <- save.dir.file.path('omega/data')
filename <- function(x) save.dir.file.path(paste0("om", x*100 ,".RData"), out)
choices <- c(0.1, 0.3, 0.5, 0.7)
for(i in choices){
  cli::cli_alert_info("Creating data for omega = {i}.")
  theta.i <- if(i == 0.1){
    c(-4.30, 0.1)
  }else if(i == 0.3){
    c(-3.00, 0.1)
  }else if(i == 0.5){
    c(-2.15, 0.1)
  }else if(i == 0.7){
    c(-1.40, 0.1)
  }else{
    cli::cli_abort("{i} not matched!")
  }
  
  cli::cli_inform("This corresponds to a theta = {theta.i}.")
  fn <- create.simfn(arguments = list(n = 500, ntms = 10,
                                      theta = theta.i))
  sim.sets <- createNsims(fn, N)
  save(sim.sets, file = filename(i))
  cli::cli_alert_success("Done, saved in:\n{filename(i)}.")
  DataSummary(sim.sets)
  cli::rule(col = "tomato2")
  rm(sim.sets)  
}
