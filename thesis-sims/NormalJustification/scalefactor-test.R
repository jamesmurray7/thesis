#' ##############################################
#' Dummy program to try and find a scale-factor #
#' to place on Sigma.hat to reduce psi_i        #  
#' ##############################################
wrap <- function(f = "gaussian"){
  obj <- function(a, sim){
    df <- do.call(c, lapply(1:NROW(sim$Walks), function(i){
      this.df <- sim$df[sim$df$id==i,]
      this.b.hat <- c(sim$b.hats[[i]])
      psi <- check.walks(sim$Walks[[i]], this.b.hat, sim$Sigmas[[i]] * a)
      psi
    }))
    abs(mean(df - 0.95))
  }
  
  sim <- create.simfn(list(f), arguments = list(n = 100L, theta = c(-3,.1)))
  X_ <- dataGen(sim)
  S <- Sample(X_, return.walks = T, TUNE = 3, force.intslope = T, b.dist = "t", df = 4,burnin = 1000, NMC = 5000L)
  
  optim(1, obj, NULL, sim = S, method = 'Brent', lower = 1e-2, upper = 2,
        control = list(abstol = 1e-3, reltol = 1e-2))$par
}

# This takes ages
r.gaussian <- replicate(100, wrap("gaussian"))
dd <- data.frame(i = 1:100L, a = r.gaussian)
qq <- quantile(dd$a)
ggplot(dd, aes(x = i, y = a)) + 
  geom_hline(yintercept = c(qq[2], qq[4]), lty = 5, alpha = .5, colour = .nice.orange, lwd = .25) + 
  geom_hline(yintercept = qq[3], lty = 5, alpha = .5, colour = 'steelblue', lwd = .25) + 
  geom_point(size = .25) + 
  scale_y_continuous(expression(a), breaks = scales::pretty_breaks(6)) + 
  scale_x_continuous("",breaks=c()) + 
  theme_csda()
ggsave(save.dir.file.path("sf.png"), width = 140, height = 60, units = "mm")
write.table(qq, file = save.dir.file.path("sf.txt"))

# Quick look -> These look okay and about 0.75, too.
wrap("poisson")
wrap("Gamma")