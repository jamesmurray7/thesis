#' ###############
#' Idea for bivariate plot
#' ###############

sim <- create.simfn()
X_ <- dataGen(sim)
WW <- Sample(X_, 5.5, T) # Ensure the walks are returned, and MH is tuned well-enough.

# Take first id for now
test <- WW$df[WW$df$id==1,]
wide <- as.data.frame(cbind("b0"=test[test$var=="b[0]",c(b0 = "condx")],
              "b1"=test[test$var=="b[1]",c(b1 = "condx")]))

library(ggplot2)
ggplot(test,aes(x = ApproxBias)) + 
  geom_boxplot() +  # Not sure what this _actually_ shows.
  facet_wrap(~var)

b0 <- test[test$var=="b[0]", "condx"]
b1 <- test[test$var=="b[1]", "condx"]
rx <- range(b0); ry <- range(b1)
b.hat <- c(test[test$var=="b[0]", "hatb"][1], test[test$var=="b[1]", "hatb"][1])
S <- WW$Sigmas[[1]]

# Theoretical with {b.hat, Sigma.hat} -->
ran <- mvtnorm::rmvnorm(1e3, b.hat, S)
ran <- as.data.frame(ran)
ggplot(ran, aes(x = V1, y = V2)) +
  geom_point(data=as.data.frame(WW$Walks[[1]]), aes(x = V1, y = V2), alpha = .15) + 
  geom_density_2d() + # I'm not sure what these contour lines actually show?
  theme_csda()        # Docu says density estimation, but again what is this?

# Filled one looks quite cool, a bit useless though.
ggplot(ran, aes(x = V1, y = V2)) +
  geom_density_2d_filled() +
  geom_point(data=as.data.frame(WW$Walks[[1]]), aes(x = V1, y = V2), alpha = .15) + 
  theme_csda()

# Decent random scatter, some correlation (as expected?)
plot(WW$Walks[[1]], pch = 20)

# I think this is a decent little plot -->
# colours from https://www.colorhexa.com/3cb371
ggplot(as.data.frame(WW$Walks[[1]]), aes(x = V1, y = V2)) + 
  geom_point(alpha = .15, size = 0.01) + 
  geom_point(data = data.frame(x=b.hat[1], y=b.hat[2]), aes(x = x, y = y),
             colour = "#b33c7e", size = 0.3, pch = 18) + 
  stat_ellipse(data = ran, type = "norm", colour = "#b3713c", level = 0.90, lwd = 0.25) + 
  stat_ellipse(data = ran, type = "norm", colour = "#3cb371", level = 0.95, lwd = 0.25) + 
  stat_ellipse(data = ran, type = "norm", colour = "#713cb3", level = 0.99, lwd = 0.25) + 
  labs(x = expression(b[0]), y = expression(b[1])) + 
  theme_csda() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 4),
    axis.text = element_text(size = 2, hjust=0, vjust=0.5),
    axis.ticks = element_blank()
  )
ggsave("~/Downloads/temp.png", width = 60, height = 60, units = "mm")  
