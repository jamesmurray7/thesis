cat("\f")

options(
	prompt = glue::glue_col("{green gmvjoint-MC> }"),
	continue = "  "
)

save.dir <- '/data/c0061461/THESIS/gmvjoint-MC'
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

.startmessage <- function(){
  cat("\f")
  message("gmvjoint - MC version")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  cat("\n")
}

source('R/jointMC.R')
source('R/EMUpdateMC.R')
source('R/inits.R')
source('R/makeDataMatrices.R')
source('R/misc.R')
source('R/survmod.R')
source('R/loglik.R')
source('R/vcovMC.R')
source('../theme_csda.R')

library(glmmTMB)
library(survival)
library(Rcpp)
library(RcppArmadillo)

.startmessage()


