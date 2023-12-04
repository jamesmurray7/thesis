cat("\f")

library(gmvjoint)
library(Rcpp)
library(RcppArmadillo)
options(
	prompt = glue::glue_col("{green SOI> }"),
	continue = "  "
)

save.dir <- '/data/c0061461/THESIS/SOI'
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

.quietly <- function(fun){
  capture.output(a <- fun(), file = nullfile())
  return(a)
}

.startmessage <- function(){
  cat("\f")
  message("Shapes of integrands (SOI) - Thesis plots")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  sourceCpp("../NormalJustification/mh.cpp")
  cat("\n")
}

source('../theme_csda.R')

.startmessage()


