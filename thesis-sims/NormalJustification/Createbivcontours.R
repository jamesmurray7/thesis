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
create.biv.contour <- function(X, id = NULL, trimmed = FALSE, # strange results when trimmed - weird!
                               axes = c("none", "left", "bottomleft", "bottom", "all", "nonebiggertitle"), N = 1e3){
  stopifnot(inherits(X, "Sample"))
  if(is.null(id)){
    cat("\n--> id not supplied, choosing a random subject...\n\n")
    id <- sample(1:NROW(X), 1)
  }
  axes <- match.arg(axes)
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
  # tit <- bquote("id: "*.(id)*", "*m[i]==.(this.mi)*".")
  tit <- bquote(m[i]==.(this.mi))
  
  # Simulate based on b.hat and Sigma.hat
  ran <- as.data.frame(mvtnorm::rmvnorm(N, this.b.hat, this.Sigma))
  names(this.Walks.df) <- names(ran)
  
  # Make the plot -->
  P <- ggplot(this.Walks.df, aes(x = b0, y = b1)) + 
    geom_point(alpha = 0.15, size = 0.01) + 
    stat_ellipse(data = ran, type = "norm", colour = "#b3713c", level = 0.90, lwd = 0.25) + 
    stat_ellipse(data = ran, type = "norm", colour = "#3cb371", level = 0.95, lwd = 0.25) + 
    stat_ellipse(data = ran, type = "norm", colour = "#713cb3", level = 0.99, lwd = 0.25) + 
    geom_point(data = as.data.frame(t(this.b.hat)), aes(x = b0, y = b1),
               colour = "#b33c7e", size = 0.3, pch = 18) + 
    labs(x = expression(b[0]), y = expression(b[1]), title = tit) + 
    theme_csda() +
    theme(
      panel.grid = element_blank(),
      title = element_text(size = 4, vjust = 0)
    )
  
  if(axes == "all"){
    P <- P + theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 4),
      axis.text = element_text(size = 2, hjust=0, vjust = 0.5),
      axis.ticks = element_blank()
    )
  }else if (axes == "left"){ # For e.g. the left-hand `wall`.
    P <- P + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 4),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  }else if (axes == "bottomleft"){
    P <- P + theme(
      axis.title.x = element_text(size = 4),
      axis.title.y = element_text(size = 4),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  }else if(axes == "bottom"){
    P <- P + theme(
      axis.title.x = element_text(size = 4),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  }else if(axes == "nonebiggertitle"){
    P <- P + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      title = element_text(size = 6, vjust = -1)
    )
  }else{ # NONE
    P <- P + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  }
  
  P
}

# text.Width      text.Height
.w <- 221.99739; .h <- 149.99825 # This is in mm
# Function to plot `staircase`-type plot
plot.staircase <- function(X, num.to.plot = 30, num.per.row = 10, file.name = NULL){
  stopifnot(inherits(X, "Sample"))
  nrow <- num.to.plot%/%num.per.row
  ncol <- num.per.row # Slightly redundant, oh well.
  if(is.null(file.name)) # Make a junk filename if none provided
    file.name <- paste0(paste(sample(LETTERS, 4), collapse = "", sep=""),
                        "_",format(Sys.time(), "%T"),".png")
  
  mis <- X$mi
  df <- data.frame(id = 1:length(mis), mis = mis, Acc = X$Acc)
  df <- df[df$Acc >= 0.20, ]       # Safeguarding against unruly fits (mainly occuring under genpois).
  df2 <- df[order(df$mis),]
  df2 <- df2[1:num.to.plot,]
  ids <- df2$id
  
  # List of plots
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x))
  file.name <- save.dir.file.path(file.name, 
                                  save.dir.file.path(family.dir.name(X$family)))
  
  # Make the plot
  png(file.name, width = floor(.w) - 1, height = floor(.h) - 1, units = "mm", res = 500)
  gridExtra::grid.arrange(grobs = Ps, nrow = nrow, ncol = ncol,
                          bottom = grid::textGrob(expression(b[0]), vjust = 0,
                                                  gp = grid::gpar(cex=.75)),
                          left = grid::textGrob(expression(b[1]), rot = 90, vjust = 1,
                                                gp = grid::gpar(cex=.75)))
  dev.off()
  cat(sprintf("\nSaved in %s.\n", file.name))
}



# "Zoomed" plot for main.tex ----------------------------------------------
# pseudo-version of plot.staircase, but only choosing six profiles.
zoom.for.main <- function(X, file.name = NULL){
  stopifnot(inherits(X, "Sample"))
  if(is.null(file.name)) # Make a junk filename if none provided
    file.name <- paste0(paste(sample(LETTERS, 4), collapse = "", sep=""),
                        "_ZOOM_",format(Sys.time(), "%T"),".png")
  mis <- X$mi
  df <- data.frame(id = 1:length(mis), mis = mis, Acc = X$Acc)
  # Ensure sensbile acceptance rates are considered.
  df <- df[df$Acc >= 0.20, ]       

  # Randomly sample two ids from each tertile of available follow-up ----
  df <- df[order(df$mis),]
  unique.mis <- unique(df$mis)
  # Sort by tertiles and take three from each?
  grps <- cut(unique.mis, quantile(unique.mis, probs=c(0,.33,.66,1)), 
              include.lowest = T, labels = F)
  # Merge this on...
  rhs <- data.frame(mis = unique.mis, grps = grps)
  df2 <- merge(df, rhs, "mis")
  
  # These are the sampled ids
  ids <- unlist(with(df2, tapply(id, grps, sample, 2, F)))
  
  # Always take at least one m[i]=1 if possible ----
  one.flag <- df$mis==1
  ids.with.one <- df2[one.flag, "id"]
  
  # Get m_i for each chosen ids
  chosen.mis <- sapply(ids, function(ii) df2[df2$id==ii,"mis"])
  ids <- ids[order(chosen.mis)]
  # Replace first subject with m_i==1 if available.
  if(all(chosen.mis>1) && length(ids.with.one) >= 1)
    ids[1] <- ids.with.one[1]
  
  # Generate list of plots ----
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x, axes = "nonebiggertitle"))
  file.name <- save.dir.file.path(file.name, 
                                  save.dir.file.path(family.dir.name(X$family)))
    
  # Make the plot ----
  png(file.name, width = 140, height = 110, units = "mm", res = 1e3)
  gridExtra::grid.arrange(grobs = Ps, nrow = 3, ncol = 2, # always plotting 6 alone
                          bottom = grid::textGrob(expression(b[0]), vjust = 0,
                                                  gp = grid::gpar(cex=.75)),
                          left = grid::textGrob(expression(b[1]), rot = 90, vjust = 1,
                                                gp = grid::gpar(cex=.75)))
  dev.off()
  cat(sprintf("\nSaved in %s.\n", file.name))
}


