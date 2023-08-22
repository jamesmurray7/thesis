randomfail.sim <- function(){
  theta <- c(sample(c(-3, -3.5, -2.75), 1), .1)
  
  fn <- create.simfn(arguments = list(n = 500L, theta = theta,
                                      random.formulas = list(~time, 
                                                             ~time,
                                                             ~1)))
  return(.quietly(fn))
}

#sims <- replicate(200L, randomfail.sim(), FALSE)
