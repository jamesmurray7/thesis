library(gmvjoint)
thesis.dir <- '/data/c0061461/THESIS'

# Unchanging stuff --------------------------------------------------------
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1+time|id),
  Y.2 ~ time + cont + bin + (1+time|id),
  Y.3 ~ time + cont + bin + (1+time|id)
)
surv.formula <- Surv(survtime, status) ~ bin

# Define function to simulate
.sim <- function(n = 250, ntms = 5, fup = 5,
                 beta = "normal-usual", D = "paper",
                 theta = "low", 
                 sigma = c(.16, .16, .16), 
                 zeta = c(0, -0.2), gamma = c(.5, -.5, .5), 
                 REs = TRUE){
  if(!class(beta)%in%"character"){ 
    beta <- beta
  }else{
    beta <- do.call(rbind, replicate(3, c(2, -.1, .1, -.2), simplify = F))
    sgn <- rbind(rep(1, 4), rep(-1, 4), rep(1, 4))
    beta <- beta * sgn
  }
  
  if(!class(D)%in%"character"){
    D <- D
  }else{
    D <- diag(c(0.25, 0.09, 0.25, 0.09, 0.25, 0.09))
    for(i in c(1,3,5))
      for(j in c(1,3,5))
        if(i != j) D[i,j] <- D[j,i] <- .125
  }
  
  if(theta == "low")
    theta <- c(-4.3, .1)
  else if(theta == "med")
    theta <- c(-2.9, .1)
  else if(theta == "high")
    theta <- c(-2.15, .1)
  else if(theta == "very-high")
    theta <- c(-1.5, .1)
  else
    stop("theta specified wrong.\n")
  
  simData(n = n, ntms = ntms, fup = fup,
          family = as.list(rep("gaussian", 3)),
          sigma = sigma, beta = beta, D = D,
          gamma = gamma, zeta = zeta, theta = theta, 
          return.ranefs = REs)
}

# Setting up "Standard": n = 250, ntms = {5, 10, 15}, 
N <- 100
generate <- function(save.dir = 'Gaussian-Standard'){
  
  if(!dir.exists(file.path(thesis.dir, save.dir))) dir.create(file.path(thesis.dir, save.dir))
  
  to.sim <- as.data.frame(expand.grid(r = seq(5,15,5),
                          omega = c("low", "med", "high", "very-high")))
  
  data <- setNames(apply(to.sim, 1, function(i){
    r = as.numeric(i[1]); theta <- i[2]
    replicate(N,
              .sim(ntms = r, theta = theta)$data,
              simplify = FALSE)
  }),
  apply(to.sim, 1, function(i) paste0("r = ", i[1], ", omega = ", i[2])))
  
  saveRDS(data, file = file.path(thesis.dir, save.dir, "data.RDS"))
  cat(sprintf("saved in %s.\n", file.path(thesis.dir, save.dir, "data.RDS")))
  
}
generate()
