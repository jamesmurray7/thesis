cat("\f")

library(gmvjoint)
options(
	prompt = glue::glue_col("{green GMVJM> }"),
	continue = "  "
)

save.dir <- '/data/c0061461/THESIS/GMVJM'
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

# Shortcuts for making D/gamma/beta for three de-facto families
.makeD <- function(family){
  if(!"list"%in%class(family)) stop("zzz")
  Ds <- lapply(family, function(f){
    switch(f,
           gaussian = diag(c(.25, .09)),
           poisson = diag(c(.50,.09)),
           binomial = matrix(2, 1, 1))
  })  
  D <- gmvjoint:::bDiag(Ds)
  off.inds <- as.data.frame(expand.grid(1:nrow(D), 1:nrow(D)))
  off.inds <- off.inds[off.inds$Var1%%2!=0 & off.inds$Var2%%2!=0&off.inds$Var1!=off.inds$Var2,]
  D[cbind(off.inds$Var1,off.inds$Var2)] <- .125
  if(gmvjoint:::is.not.SPD(D)) stop("not SPD")
  return(D)
}

.makegamma <- function(K) sapply(1:K, function(k) ifelse(k %% 2 != 0, 1, -1) * .5)
.makebeta <- function(K){ # Assumes it goes Gauss -> Poisson -> Binomial (ad inf.)
  
})

.quietly <- function(fun){
  capture.output(a <- fun(), file = nullfile())
  return(a)
}

.startmessage <- function(){
  cat("\f")
  message("GMVJM - Thesis plots")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  cat("\n")
}

source('stockSimfn.R')
source('../theme_csda.R')

.startmessage()


