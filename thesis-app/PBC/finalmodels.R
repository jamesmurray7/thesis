# ##########################
# Comparing JMb2 to gmv fits
# ##########################
source(".Rprofile")
library(xtable)
jmb.dir <- save.dir.file.path("Joint/JMb2")
(jmb.load <- dir(jmb.dir, full.names = T))
gmv.dir <- save.dir.file.path("Joint/Multivs")
(gmv.load <- dir(gmv.dir, full.names = T)[grepl("^biv|^triv", dir(gmv.dir))])


# TABULATIONS -------------------------------------------------------------
# We want to create the following table for trivariate and subseq. bivariate 
# joint models fit to PBC.

create.table <- function(ww = c("triv", "biv")){
  # Loading (local assignment)
  if(ww == "triv"){
    assign("J", get(load(save.dir.file.path("Triv.RData", jmb.dir))))
    assign("G", get(load(save.dir.file.path("triv.RData", gmv.dir))))
  }else{
    assign("J", get(load(save.dir.file.path("Biv.RData", jmb.dir))))
    assign("G", get(load(save.dir.file.path("biv.RData", gmv.dir))))
  }
  
  # xtable S3 method for `joint` gets the table we want to merge onto
  xG <- xtable::xtable(G)
  # JMbayes2:::summary.jm provides mean (SD) [95% CrI] for {beta, sigma, gamma, zeta}
  SJ <- summary(J)

  # But we need to work out these for vech(D) separately.
  jmb.D.mean <- J$statistics$Mean$D
  jmb.D.SD <- J$statistics$SD$D
  jmb.D.lb <- J$statistics$CI_low$D
  jmb.D.ub <- J$statistics$CI_upp$D
  
  .to3dp <- function(x) format(round(x, 3), nsmall = 3, justify = "right")
  vech <- gmvjoint:::vech
  xt <- xtable::xtable(G, vcov = T)
  chunks <- lapply(1:G$ModelInfo$K, function(i){
    x <- G$ModelInfo$inds$R$b[[i]]
    # Start with D
    Dx <- SJ$D[x,x]; vDx <- vech(Dx)
    inds <- which(vech(SJ$D) %in% vDx)
    xSD <- jmb.D.SD[inds]
    xLB <- jmb.D.lb[inds]
    xUB <- jmb.D.ub[inds]
    Parameterjmb <- paste0("D_{", i,",",apply(which(lower.tri(Dx, T), arr.ind = T) - 1, 1, paste, collapse=''),"}")
    MSD <- paste0(.to3dp(vDx), " (", .to3dp(xSD), ")")
    CrI <- paste0('[', .to3dp(xLB),', ',.to3dp(xUB),']')
    
    JMbD <- cbind(Parameterjmb, MSD, CrI)
    # Now beta/sigma(^2)/gamma (in that order)
    gammas <- SJ$Survival[i+3, c("Mean", "StDev", "2.5%", "97.5%"), drop=F]
    lookup <- paste0("Outcome",i)
    Out <- SJ[[lookup]][,c("Mean", "StDev", "2.5%", "97.5%")]
    if(G$ModelInfo$family[[i]]=="gaussian"){ # Simply square to get sigma^2
      Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")] <- Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")]^2
    }else{ # Utilise geometric mean to transform SD
      Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")] <- log(Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")])
      aa <- do.call(rbind, J$mcmc$sigmas)
      exp(log(J$statistics$Mean$sigmas[i])) * (sd(log(aa[,i]))/sqrt(nrow(aa)-1))
      Out[grepl("sigma", row.names(Out)), c("StDev")] <- exp(log(J$statistics$Mean$sigmas[i])) * (sd(log(aa[,i]))/sqrt(nrow(aa)-1))
      Out <- rbind(Out, dummy = c(Mean = 0, StDev = 0, `2.5%` = 0, `97.5%` = 0))
    }
    this <- rbind(Out, gammas)
    MSD <- paste0(.to3dp(this$Mean), ' (', .to3dp(this$StDev), ')')
    CrI <- paste0('[',.to3dp(this$`2.5%`), ', ', .to3dp(this$`97.5%`), ']')
    Parameterjmb <- row.names(this)
    JMbrest <- cbind(Parameterjmb, MSD, CrI)
    out <- cbind(xt$RespChunks[[i]], rbind(JMbD, JMbrest))
    if(i == G$ModelInfo$K){
      zeta <- SJ$Survival[1:length(G$coeffs$zeta),c("Mean", "StDev", "2.5%", "97.5%"),drop=F]
      MSD <- paste0(.to3dp(zeta$Mean), ' (', .to3dp(zeta$StDev), ')')
      CrI <- paste0('[',.to3dp(zeta$`2.5%`), ', ', .to3dp(zeta$`97.5%`), ']')
      out <- rbind(out, cbind(xt$zeta, Parameterjmb=row.names(zeta), MSD, CrI))
    }
    out
  })
  
  to.xt2 <- do.call(rbind, chunks)
  to.xt2 <- to.xt2[,-4]
  to.xt2[,1] <- paste0("$\\", to.xt2[,1], "$")
  xt.to.xt <- xtable(to.xt2)
  print(xt.to.xt, 
        include.rownames = FALSE,
        sanitize.text.function = identity)
}

create.table("triv")
create.table("biv")
# FIGURE(S) ---------------------------------------------------------------
# Want to make a plot for survival (?) part only (?: maybe do everything)







