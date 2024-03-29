cat("\f")

library(gmvjoint)
options(
	prompt = glue::glue_col("{green Quadrature> }"),
	continue = "  "
)

save.dir <- '/data/c0061461/THESIS/quadnodes'
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

# Shortcuts for making D/gamma/beta for all families
.makeD <- function(family){
  if(!"list"%in%class(family)){
    if(length(family) == 1L)
      family <- list(family)
    else
      stop("zzz")
  }
  Ds <- lapply(family, function(f){
    if(f == "gaussian"){
      return(diag(c(.25, .09)))
    }else if(f == "Gamma" | f == "genpois"){
      return(diag(c(.20, .05)))
    }else if(f %in% c("poisson", "negbin")){
      return(diag(c(.50, .09)))
    }else{
      return(matrix(2, 1, 1))
    }
  })
  D <- gmvjoint:::bDiag(Ds)
  # Populate covariance between random intercepts.
  off.inds <- which(D==.25|D==.20|D==.50|D==2, arr.ind = T)[,1]
  off.inds <- as.data.frame(expand.grid(off.inds, off.inds))
  off.inds <- off.inds[off.inds$Var1!=off.inds$Var2,]
  D[cbind(off.inds$Var1,off.inds$Var2)] <- .125
  if(gmvjoint:::is.not.SPD(D)) stop("not SPD")
  return(D)
}

.makegamma <- function(K) sapply(1:K, function(k) ifelse(k %% 2 != 0, 1, -1) * .5)
.makebeta <- function(family){ 
  if(!"list"%in%class(family)){
    if(length(family) == 1L)
      family <- list(family)
    else
      stop("zzz")
  }
  betas <- lapply(family, function(f){
    if(f == "gaussian"){
      return(c(-2,.1,-.1,.2))
    }else if(f == "Gamma"){
      return(c(0, -0.1, 0.1, -0.2))
    }else if(f == "genpois"){
      return(c(1, .05, -.05, .1))
    }else if(f %in% c("poisson", "negbin")){
      return(c(2,-.1,.1,.2))
    }else{
      return(c(1,-1,1,-1))
    }
  })
  do.call(rbind, betas)
}

.makesigma <- function(family){
  if(!"list"%in%class(family)){
    if(length(family) == 1L)
      family <- list(family)
    else
      stop("zzz")
  }
  lapply(family, function(f){
    if(f == "gaussian")
      return(0.16)
    else if(f == "negbin")
      return(1)
    else if(f == "genpois")
      return(-0.2)
    else if(f == "Gamma")
      return(2)
    else
      return(0)
  })
}

.quietly <- function(fun){
  capture.output(a <- fun(), file = nullfile())
  return(a)
}

.startmessage <- function(){
  cat("\f")
  message("GMVJM - Thesis plots - Quadrature nodes investigation")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  cat("\n")
}

source('stockSimfn.R')
source("datagen1.R")
source('../theme_csda.R')

.startmessage()
