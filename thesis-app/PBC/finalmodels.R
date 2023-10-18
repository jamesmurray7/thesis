# ##########################
# Comparing JMb2 to gmv fits
# ##########################
source(".Rprofile")
library(xtable)
jmb.dir <- save.dir.file.path("Joint/JMb2")
(jmb.load <- dir(jmb.dir, full.names = T))
gmv.dir <- save.dir.file.path("Joint/Multivs")
(gmv.load <- dir(gmv.dir, full.names = T)[grepl("^biv|^triv", dir(gmv.dir))])
# misc fns
.to3dp <- function(x) format(round(x, 3), nsmall = 3, justify = "right")
vech <- gmvjoint:::vech
qz <- qnorm(.975)


tab.from.joineRML <- function(G, M){
  # Check
  stopifnot(inherits(G, "joint") && inherits(M, "mjoint"))
  # summary.mjoint returns Estimate (SE) [95% CI] for {beta, sigma^2, gamma, zeta}
  SM <- summary(M)
  VarM <- vcov(M)
  # Work out vech(D)
  jML.D <- SM$D
  jML.vD <- vech(SM$D)
  # Work out its standard error from hessian
  jML.vD.SE <- sqrt(diag(VarM))[1:G$ModelInfo$Pcounts$vD]
  jML.vD.lb <- jML.vD - qz * jML.vD.SE
  jML.vD.ub <- jML.vD + qz * jML.vD.SE
  
  xt <- xtable::xtable(G, vcov = T)
  ML.long <- SM$coefs.long; rn <- row.names(ML.long)
  chunks <- lapply(1:G$ModelInfo$K, function(i){
    x <- G$ModelInfo$inds$R$b[[i]]
    # Start with D
    Dx <- SM$D[x,x]; vDx <- vech(Dx)
    inds <- which(jML.vD %in% vDx)
    xSE <- jML.vD.SE[inds]
    xLB <- jML.vD.lb[inds]
    xUB <- jML.vD.ub[inds]
    ParameterjML <- paste0("D_{", i,",",apply(which(lower.tri(Dx, T), arr.ind = T) - 1, 1, paste, collapse=''),"}")
    MSE <- paste0(.to3dp(vDx), " (", .to3dp(xSE), ")")
    CI <- paste0('[', .to3dp(xLB),', ',.to3dp(xUB),']')
    
    jMLD <- cbind(ParameterjML, MSE, CI)
    # Now beta / sigma2 / gamma 
    gammas <- SM$coefs.surv[(i+length(G$coeffs$zeta)), 1:2, drop=F]
    rows <- which(grepl(paste0("\\_", i, "$"), rn))
    # This is beta only...
    Out <- ML.long[rows,1:2,drop=F]
    # Extract sigma^2 by itself
    sigma2 <- cbind(Value = M$coefficients$sigma2[i],
                    `Std. Err` = sqrt(VarM[grepl(paste0("sigma2_", i), row.names(VarM)),
                                      grepl(paste0("sigma2_", i), colnames(VarM))]))
    this <- rbind(Out, sigma2, gammas)
    CI <- cbind(`2.5%` = this[,1] - qz * this[,2],
                `97.5%` = this[,1] + qz * this[,2])
    MSE <- paste0(.to3dp(this[,1]), " (", .to3dp(this[,2]), ")")
    CI <- paste0("[", .to3dp(CI[,1]), ", ", .to3dp(CI[,2]), "]")
    parameterjML <- row.names(this)
    Out <- cbind(parameterjML, MSE, CI)
    Out <- cbind(xt$RespChunks[[i]], rbind(jMLD, Out))
    if(i == G$ModelInfo$K){
      zeta <- SM$coefs.surv[1:length(G$coeffs$zeta),1:2,drop=F]
      MSE <- paste0(.to3dp(zeta[,1]), ' (', .to3dp(zeta[,2]), ')')
      CI <- paste0("[", .to3dp(zeta[,1] - qz * zeta[,2]), ", ", 
                        .to3dp(zeta[,1] + qz * zeta[,2]), ']')
      Out <- rbind(Out, cbind(xt$zeta, cbind(parameterjML = row.names(zeta), MSE, CI)))
    }
    Out
  })
  chunks
}

tab.from.JMbayes2 <- function(G, J){
  # Check
  stopifnot(inherits(G, "joint") && inherits(J, "jm"))
  
  # JMbayes2:::summary.jm provides mean (SD) [95% CrI] for {beta, sigma, gamma, zeta}
  SJ <- summary(J)
  
  # But we need to work out these for vech(D) separately.
  jmb.D.mean <- J$statistics$Mean$D
  jmb.D.SD <- J$statistics$SD$D
  jmb.D.lb <- J$statistics$CI_low$D
  jmb.D.ub <- J$statistics$CI_upp$D
  
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
    gammas <- SJ$Survival[i+length(G$coeffs$zeta), c("Mean", "StDev", "2.5%", "97.5%"), drop=F]
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
  chunks
}

# TABULATIONS -------------------------------------------------------------
# We want to create the following table for trivariate and subseq. bivariate 
# joint models fit to PBC.
create.table <- function(ww = c("triv", "biv")){
  # Loading (local assignment)
  if(ww == "triv"){
    # `Jmbayes2`
    assign("J", get(load(save.dir.file.path("Triv.RData", jmb.dir))))
    # `Gmvjoint`
    assign("G", get(load(save.dir.file.path("triv.RData", gmv.dir))))
  }else{
    # `Jmbayes2`
    assign("J", get(load(save.dir.file.path("Biv.RData", jmb.dir))))
    # `Gmvjoint`
    assign("G", get(load(save.dir.file.path("biv.RData", gmv.dir))))
    # `joinerMl`
    assign("M", get(load(save.dir.file.path("bivjML.RData", gmv.dir))))
  }

  to.xt2 <- do.call(rbind, tab.from.JMbayes2(G, J))
  if(ww != "triv"){
    chunks.jML <- tab.from.joineRML(G, M)
    to.xt2 <- cbind(to.xt2, do.call(rbind, chunks.jML)[,-c(1:4)])
  }
    
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







