# ###########################################################################
# Createbivcontours.R -------------
#
# The idea as of 22/09/23:
# For some simulated random walks, we want to randomly sample n subjects 
# and produce n plots in same manner as in bivcontouridea.
# As of now I'm not sure if this replaces e.g. the plots from the GLMM 
# paper! -------> Will ask Pete on 28/09/23.
#
# Thesis page dims to keep in mind: 
# Text width:  149.99825mm; 
# Text height: 221.99739mm;
# Might make more sense landscape?
#
# ###########################################################################
library(ggplot2)

# Function to create _one_ plot for _one_ subject.
create.biv.contour <- function(X, id, trimmed = FALSE, axes = TRUE, N = 1e3){
  stopifnot(inherits(X, "Sample"))
  
  # Get everything we need associated with this id
  this.Walks <- X$Walks[[id]]
  if(trimmed) this.Walks <- trim.Walks(this.Walks)
  this.Walks.df <- as.data.frame(this.Walks)
  this.Acc <- X$Acc[id]; this.mi <- X$mi[id]
  this.df <- X$df[X$df$id==id,]
  this.b.hat <- c(b0 = this.df[this.df$var=="b[0]", "hatb"][1], b1 = this.df[this.df$var=="b[1]", "hatb"][2])
  this.b.true <- X$true.b[id,]
  this.Sigma <- X$Sigmas[[id]]
  
  # Make plot title (probably change in future)
  tit <- bquote("id: "*.(id)*", "*m[i]==.(this.mi)*".")
  
  # Simulate based on b.hat and Sigma.hat
  ran <- as.data.frame(mvtnorm::rmvnorm(N, this.b.hat, this.Sigma))
  names(this.Walks.df) <- names(ran)
  
  # Make the plot -->
  P <- ggplot(this.Walks.df, aes(x = b0, y = b1)) + 
    stat_ellipse(data = ran, type = "norm", colour = "#b3713c", level = 0.90, lwd = 0.25) + 
    stat_ellipse(data = ran, type = "norm", colour = "#3cb371", level = 0.95, lwd = 0.25) + 
    stat_ellipse(data = ran, type = "norm", colour = "#713cb3", level = 0.99, lwd = 0.25) + 
    geom_point(alpha = 0.15, size = 0.01) + 
    geom_point(data = as.data.frame(t(this.b.hat)), aes(x = b0, y = b1),
               colour = "#b33c7e", size = 0.3, pch = 18) + 
    labs(x = expression(b[0]), y = expression(b[1]), title = tit) + 
    theme_csda() +
    theme(
      panel.grid = element_blank()
    )
  
  if(axes){
    P <- P + theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 4),
      axis.text = element_text(size = 2, hjust=0, vjust = 0.5),
      axis.ticks = element_blank()
    )
  }else{
    P <- P + theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  }
  
  P
}
