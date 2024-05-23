cat("\f")

library(gmvjoint)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist) # Don't know if I actually need to load this in?
options(
	prompt = glue::glue_col("{green NormalJustification> }"),
	continue = "  "
)

.valid.families <- c("gaussian", "poisson", "binomial",
                     "genpois", "negbin", "gamma")
.check.family <- function(x) tolower(x) %in% .valid.families
toSentence <- function(x) tools::toTitleCase(x)
save.dir <- ifelse(grepl("macOS", utils::sessionInfo()$running),
                   './',
                   '/data/c0061461/THESIS/')
family.dir.name <- function(f = "Unspecified"){
	if(f == "Unspecified") cat("\n\n!!!!!\nMake sure f is specified\n!!!!!\n\n")
  stopifnot(.check.family(f))
	f <- tolower(f)
	ff <- switch(f,
	             gaussian = toSentence("gaussian"),
	             poisson = toSentence("poisson"),
	             negbin = toSentence("negativebinomial"),
	             binomial = toSentence("binomial"),
	             genpois = toSentence("generalisedpoisson"),
	             gamma = toSentence("gamma")
  )	
	paste0(ff,"-Normal-Justification")		
}

save.dir.file.path <- function(x, LHS = save.dir) gsub("\\/\\/","/",file.path(LHS, x))

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
    }else if(f == "Gamma"){
      return(diag(c(.20, .05)))
    }else if(f %in% c("poisson", "negbin")){
      return(diag(c(.50, .09)))
    }else if(f == "genpois"){
      return(matrix(0.30, 1, 1))
    }else{
      return(matrix(2, 1, 1))
    }
  })
  D <- gmvjoint:::bDiag(Ds)
  # Populate covariance between random intercepts.
  off.inds <- which(D==.25|D==.20|D==.50|D==2|D==0.30, arr.ind = T)[,1]
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
      return(-0.3)
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

.plot.both <- function(X,Y){
  stopifnot(inherits(X, "Sample") & inherits(Y, "Sample"))
  par(mfrow=c(1,2))
  plot(X)
  plot(Y)
  par(mfrow=c(1,1))
}

.startmessage <- function(){
  cat("\f")
  message("Normal justification - Thesis plots")
  cat("Current existing directories:\n")
  dirs.to.show <- gsub("\\/\\/","/",
                       dir(save.dir, pattern = "Justification$",
                           full.names = T)
                      )
  invisible(sapply(dirs.to.show, fs::dir_tree))
  cat("\n")
}

source('stockSimfn.R')
source("Createbivcontours.R")
source("Sample.R")
source('../theme_csda.R')
source("coverage.R")
cat("\nLoading c++ file...\n")
sourceCpp("mh.cpp")

.startmessage()


