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
                               axes = c("none", "left", "bottomleft", "bottom", "all", "nonebiggertitle"), N = 1e3,
                               level = 0.95, df = 2){
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
  this.b.hat <- setNames(c(X$b.hats[[id]]), c("b0", "b1"))
  this.b.true <- X$true.b[id,]
  this.Sigma <- X$Sigmas[[id]]
  
  # Make plot title (probably change in future)
  # tit <- bquote("id: "*.(id)*", "*m[i]==.(this.mi)*".")
  
  # Simulate based on b.hat and Sigma.hat
  ran <- as.data.frame(mvtnorm::rmvnorm(N, this.b.hat, this.Sigma))
  # test <- car::dataEllipse(ran[,1],ran[,2], draw = F, segments = 150)
  # test <- test$`0.95`
  
  # Trying things...
  # https://stackoverflow.com/questions/76175820/function-to-calculate-whether-a-point-is-within-ellipse-using-stat-ellipse
  # eig <- eigen(this.Sigma)
  eig <- eigen(cov(ran))
  rx <- sqrt(eig$values[1] * qchisq(level, df))
  ry <- sqrt(eig$values[2] * qchisq(level, df))
  rx2 <- rx^2; ry2 <- ry^2
  theta <- atan2(eig$vec[2,1], eig$vectors[1,1])
  ct <- cos(theta); st <- sin(theta)
  cm <- colMeans(ran)
  dx <- this.Walks[,1] - cm[1]
  dy <- this.Walks[,2] - cm[2]
  checks <- (ct * dx + st * dy)^2/rx2 + (st * dx - ct * dy)^2/ry2 <= 1
  pc.checks <- round(sum(checks)/length(checks), 3)
  tit <- bquote(m[i]==.(this.mi)*","~psi[i]==.(pc.checks))
  # Plot -->
  # png("~/Downloads/temp.png", width = 140, height = 90, units = "mm", res = 1e3)
  # plot(test[,2]~test[,1], type="l",
  #      ylab = expression(b[1]), xlab = expression(b[0]),
  #      ylim = range(this.Walks[,2]), xlim = range(this.Walks[,1]))
  # # abline(v=this.b.hat[1],h=this.b.hat[2])
  # points(this.Walks[checks,1], this.Walks[checks,2], pch = 19, col = 'mediumseagreen', cex = .25)
  # points(this.Walks[!checks,1], this.Walks[!checks,2], pch = 19, col = 'red', cex = .25)
  # points(this.b.hat[2] ~ this.b.hat[1], col = 'blue', pch = 4, cex = .5)
  # legend("topright", bty = "n", legend = sprintf("%.2f%%", 100 * mean(checks)))
  # dev.off()
  
  names(this.Walks.df) <- names(ran)
  
  # Make the plot -->
  P <- ggplot(this.Walks.df, aes(x = b0, y = b1)) + 
    geom_point(alpha = 0.15, size = 0.01) + 
    # stat_ellipse(data = ran, type = "norm", colour = "#b3713c", level = 0.90, lwd = 0.25) + 
    # stat_ellipse(data = ran, type = "norm", colour = "#713cb3", level = 0.99, lwd = 0.25) + 
    stat_ellipse(data = ran, type = "norm", colour = "#3cb371", level = 0.95, lwd = 0.25) + 
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
plot.staircase <- function(X, num.to.plot = 30, num.per.row = 10, file.name = NULL, 
                           N = 1e3, level = 0.95, df = 2){
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
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x, N = N))
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

compare.theor.distn <- function(X1, X2){
  stopifnot(inherits(X1, "Sample") & inherits(X2, "Sample"))
  if(sum(X1$include.survival, X2$include.survival) != 1)
    stop("Need one Sample with include.survival = T, not zero/both!")
  
  b.hats1 <- X1$b.hats; b.hats2 <- X2$b.hats
  Sigma.hats1 <- X1$Sigmas
  eig.S1 <- lapply(Sigma.hats1, eigen)
  theta.S1 <- sapply(eig.S1, function(x){
    xx <- atan2(x$vec[2,1], x$vectors[1,1])
    if(xx < 0) xx <- xx + 2 * pi
    xx
  })
  semi.maj1 <- sapply(eig.S1, function(x) sqrt(x$values[1] * qchisq(0.95, 2)))
  semi.min1 <- sapply(eig.S1, function(x) sqrt(x$values[2] * qchisq(0.95, 2)))
  Sigma.hats2 <- X2$Sigmas
  eig.S2 <- lapply(Sigma.hats2, eigen)
  theta.S2 <- sapply(eig.S2, function(x){
    xx <- atan2(x$vec[2,1], x$vectors[1,1])
    if(xx < 0) xx <- xx + 2 * pi
    xx
  })
  semi.maj2 <- sapply(eig.S2, function(x) sqrt(x$values[1] * qchisq(0.95, 2)))
  semi.min2 <- sapply(eig.S2, function(x) sqrt(x$values[2] * qchisq(0.95, 2)))
  # Keep monitoring acceptance probs!
  A1 <- X1$Acc; A2 <- X2$Acc
  
  if(X2$include.survival){
    loc.diff <- Map(function(b2,b1) b2 - b1, b2 = b.hats2, b1 = b.hats1)
    ang.diff <- theta.S2 - theta.S1
    maj.diff <- semi.maj2 - semi.maj1
    min.diff <- semi.min2 - semi.min1
    Acc.diff <- A2 - A1
  }else{
    loc.diff <- Map(function(b1,b2) b1 - b2, b1 = b.hats1, b2 = b.hats2)
    ang.diff <- theta.S1 - theta.S2
    maj.diff <- semi.maj1 - semi.maj2
    min.diff <- semi.min1 - semi.min2
    Acc.diff <- A1 - A2
  }
  
  structure(list(loc.diff = loc.diff,
                 ang.diff = ang.diff,
                 maj.diff = maj.diff,
                 min.diff = min.diff,
                 Acc.diff = Acc.diff,
                 Acc = rbind(A1, A2),
                 mi = X1$mi,
                 which.surv = if(X1$include.survival) "X1" else "X2"),
            class = "dist.compare")
}

print.dist.compare <- function(x){
  stopifnot(inherits(x, "dist.compare"))
  # Strictly between 20-30%
  acc <- apply(x$Acc,2,function(y){
    y1 <- round(y[1],2); y2 <- round(y[2],2)
    c1 <- y1 >= 0.20 & y1 <= 0.30
    c2 <- y2 >= 0.20 & y2 <= 0.30
    c1 & c2
  })
  cat(sprintf("\n----\n%s is the survival density\n----\n", x$which.surv))
  cat(sprintf("\n%.2f%% instances with two usable samples\n", 100 * mean(acc)))
  cat("\nSummary of b.hat differences:\n ")
  print(summary(setNames(as.data.frame(do.call(rbind, x$loc.diff)), c("b0", "b1"))[acc,]))
  cat("\nSummary of angle differences:\n")
  print(summary(x$ang.diff[acc]))
  cat("\nSummary of semi-major axes differences:\n")
  print(summary(x$maj.diff[acc]))
  cat("\nSummary of semi-minor axes differences:\n")
  print(summary(x$min.diff[acc]))
}

plot.dist.compare <- function(x){
  stopifnot(inherits(x, "dist.compare"))
  acc <- apply(x$Acc, 2, function(y){
    y1 <- round(y[1],2); y2 <- round(y[2],2)
    c1 <- y1 >= 0.20 & y1 <= 0.30
    c2 <- y2 >= 0.20 & y2 <= 0.30
    c1 & c2
  })
  loc.df <- setNames(as.data.frame(do.call(rbind, x$loc.diff)), c("b0", "b1"))[acc,]
  loc.df$mi <- x$mi[acc]
  P1 <- ggplot(loc.df, aes(x = b0, y = b1, colour = mi)) + 
    geom_point(size=.5) + 
    labs(x = expression(b[0]^{"S"}-b[0]^{"NS"}),
         y = expression(b[1]^{"S"}-b[1]^{"NS"}),
         colour = expression(m[i])) + 
    scale_colour_gradient(low = .nice.orange, high = "red3",
                          breaks = c(2,4,6,8,10)) + 
    theme_csda()+
    theme(
      legend.position = "none",
      axis.text = element_text(size = 5)
    )
  
  r.df <- data.frame(mi = x$mi, ry = x$min.diff, rx = x$maj.diff)[acc,]
  P2 <- ggplot(r.df, aes(x = rx, y = ry, colour = mi)) + 
    geom_point(size=.5) + 
    labs(x = expression(r[x]^{"S"}-r[x]^{"NS"}),
         y = expression(r[y]^{"S"}-r[y]^{"NS"}),
         colour = expression(m[i])) + 
    scale_color_gradient(low = .nice.orange, high = "red3",
                         breaks = c(2,4,6,8,10)) + 
    theme_csda() + 
    theme(
      legend.key.width = unit(3, "mm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      axis.text = element_text(size = 5)
    )
  
  ang.df <- data.frame(mi = x$mi, theta = x$ang.diff)[acc,]
  P3 <- ggplot(ang.df, aes(x = mi, y = theta, fill = mi, group = mi)) + 
    geom_hline(yintercept = 0, lty = 5, alpha = 0.33) +
    geom_boxplot(outlier.alpha = 0.33, lwd = 0.25, 
                 fatten = 2, outlier.size = 0.50) + 
    labs(x = expression(m[i]),
         y = expression(vartheta^{"S"}-vartheta^{"NS"})) + 
    scale_fill_gradient(low = .nice.orange, high = "red3") + 
    scale_x_continuous(breaks = 1:10) + 
    theme_csda() + 
    theme(
      legend.position = 'none',
    )
  
  gridExtra::grid.arrange(P1, P2, P3,
                          layout_matrix = rbind(
                            c(1,1,2,2),
                            c(1,1,2,2),
                            c(3,3,3,3),
                            c(3,3,3,3)
                          ))
}













# OLD CODE ----------------------------------------------------------------
# "Zoomed" plot for main.tex ----------------------------------------------
# EDIT -> don't actually like these, so scrapping this!!
# pseudo-version of plot.staircase, but only choosing six profiles.
zoom.for.main <- function(X, file.name = NULL, N = 1e4){
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
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x, axes = "nonebiggertitle", N = N))
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


