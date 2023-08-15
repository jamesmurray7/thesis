cat("\f")

library(gmvjoint)
options(
	prompt = glue::glue_col("{green MVJM> }"),
	continue = "  "
)

save.dir <- '/data/c0061461/THESIS/MVJM'
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

# Shortcuts for making D/gamma/beta
.makeD <- function(K){
  D <- diag(rep(c(0.25, 0.09), K))
  off.inds <- as.data.frame(expand.grid(1:(2*K), 1:(2*K)))
  off.inds <- off.inds[off.inds$Var1 %% 2 != 0 & off.inds$Var2 %% 2 != 0 & off.inds$Var1 != off.inds$Var2,]
  D[cbind(off.inds$Var1, off.inds$Var2)] <- .125
  D
}

.makegamma <- function(K) sapply(1:K, function(k) ifelse(k %% 2 != 0, 1, -1) * .5)
.makebeta <- function(K) t(sapply(1:K, function(k) ifelse(k %% 2 != 0, 1, -1) * c(2, -0.1, 0.1, -0.2)))

.quietly <- function(fun){
  capture.output(a <- fun(), file = nullfile())
  return(a)
}

.startmessage <- function(){
  cat("\f")
  message("MVJM - Thesis plots")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  cat("\n")
}

source('stockSimfn.R')
source('../theme_csda.R')

.startmessage()


