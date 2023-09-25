# Test program to produce e.g. an a4r png image to be placed in thesis appendices
rm(list=ls())
gc()
source(".Rprofile")
# Sticking with Gaussian
sim <- create.simfn()
X_ <- dataGen(sim)
X <- Sample(X_, 5.5, T) # Ensure the walks are returned, and MH is tuned well-enough.

# Continue as if landscape
w <- 221.99739; h <- 149.99825 # This is in mm
(area <- w*h)                  # mm for whole page?
(area.p.s <- area/X_$args$n)   # mm (area) for each subject?
sqrt(area.p.s)                 # 18.24807 w/h if square(!) -> tiny, maybe just do 30-50 subjects?

# Simply make a grid, as I think this resizes automatically?
# gridExtra::grid.arrange fills in ACROSS rows first!
plot.random.ids <- function(X, num.to.plot = 30, num.per.row = 10,
                             file.name = NULL){
  stopifnot(inherits(X, "Sample"))
  nrow <- num.to.plot%/%num.per.row
  ncol <- num.per.row # Slightly redundant, oh well.
  if(is.null(file.name)) # Make a junk filename if none provided
    file.name <- paste0(paste(sample(LETTERS, 4), collapse = "", sep=""),
                        "_",format(Sys.time(), "%T"),".png")
  
  # Random sample
  ids <- sample(1:NROW(X$Walks), num.to.plot)
  # Order these by their profile length to make plot a bit nicer.
  ids <- ids[order(X$mi[ids])]
  
  # List of plots
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x))
  file.name <- save.dir.file.path(file.name, 
                                  save.dir.file.path(family.dir.name(X$family)))
  # Make the plot
  png(file.name, width = floor(w)-1, height = floor(h)-1, units = "mm", res = 500)
  gridExtra::grid.arrange(grobs = Ps, nrow = nrow, ncol = ncol,
                          bottom = grid::textGrob(expression(b[0]), vjust = 0),
                          left = grid::textGrob(expression(b[1]), rot = 90, vjust = 1))
  dev.off()
  cat(sprintf("\nSaved in %s.\n", file.name))
}
print.random.ids(X)

# This has _increasing_ m[i] guaranteed.
plot.staircase <- function(X, num.to.plot = 30, num.per.row = 10, file.name = NULL){
  stopifnot(inherits(X, "Sample"))
  nrow <- num.to.plot%/%num.per.row
  ncol <- num.per.row # Slightly redundant, oh well.
  if(is.null(file.name)) # Make a junk filename if none provided
    file.name <- paste0(paste(sample(LETTERS, 4), collapse = "", sep=""),
                        "_",format(Sys.time(), "%T"),".png")
  
  mis <- X$mi
  df <- data.frame(id = 1:length(mis), mis = mis)
  df2 <- df[order(df$mis),]
  df2 <- df2[1:num.to.plot,]
  ids <- df2$id
  
  # List of plots
  Ps <- lapply(ids, function(x) create.biv.contour(X, id = x))
  file.name <- save.dir.file.path(file.name, 
                                  save.dir.file.path(family.dir.name(X$family)))
  
  # Make the plot
  png(file.name, width = floor(w)-1, height = floor(h)-1, units = "mm", res = 500)
  gridExtra::grid.arrange(grobs = Ps, nrow = nrow, ncol = ncol,
                          bottom = grid::textGrob(expression(b[0]), vjust = 0,
                                                  gp = grid::gpar(cex=.75)),
                          left = grid::textGrob(expression(b[1]), rot = 90, vjust = 1,
                                                gp = grid::gpar(cex=.75)))
  dev.off()
  cat(sprintf("\nSaved in %s.\n", file.name))
}
