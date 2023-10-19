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

.xx <- function(x){
  substr(x,1,1)<- toupper(substr(x,1,1))
  x
}


tab.from.joineRML <- function(G, M){
  # Check
  stopifnot(inherits(G, "joint") && inherits(M, "mjoint"))
  # summary.mjoint returns Estimate (SE) [95% CI] for {beta, sigma^2, gamma, zeta}
  SM <- summary(M)
  VarM <- joineRML:::vcov.mjoint(M)
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
create.table <- function(ww = c("triv", "biv"), reduced = FALSE){
  # Loading (local assignment)
  suff <- if(reduced) 2 else NULL
  if(ww == "triv"){
    # `Jmbayes2`
    assign("J", get(load(save.dir.file.path("Triv.RData", jmb.dir))))
    # `Gmvjoint`
    assign("G", get(load(save.dir.file.path("triv.RData", gmv.dir))))
  }else{
    file.name <- paste0("biv", suff, ".RData")
    # `Jmbayes2`
    assign("J", get(load(save.dir.file.path(.xx(file.name), jmb.dir))))
    # `Gmvjoint`
    assign("G", get(load(save.dir.file.path(file.name, gmv.dir))))
    # `joinerMl`
    file.name <- gsub("\\.RData", "jML.RData", file.name)
    assign("M", get(load(save.dir.file.path(file.name, gmv.dir))))
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

create.table("triv", FALSE)
# Before removing sex from surv. model
create.table("biv", FALSE)
# After removing sex from surv. model
create.table("biv", TRUE)


# FIGURE(S) ---------------------------------------------------------------
# Create _one_ for the final bivariate joint model...
load(save.dir.file.path("Biv2.RData", jmb.dir)) # `J3`
load(save.dir.file.path("biv2.RData", gmv.dir)) # `joint.biv2`
load(save.dir.file.path("biv2jML.RData", gmv.dir)) # `joint.biv2.joineRML`

# This gets everything needed, just need to extract relevant stuff!
one <- do.call(rbind, tab.from.JMbayes2(joint.biv2, J3))
two <- do.call(rbind, tab.from.joineRML(joint.biv2, joint.biv2.joineRML))

all3 <- as.data.frame(rbind(one[,1:3],                              # Approx. EM
      cbind(one[1:nrow(one), 1], one[,5:6]),  # JMbayes2
      cbind(two[1:nrow(two), 1], two[,5:6])   # joineRML
))

all3$Method <- c(rep("Approximate EM", nrow(one)), 
                 rep("JMbayes2", nrow(one)),
                 rep("joineRML", nrow(two)))

all3$Method <- factor(all3$Method, unique(all3$Method))

# Let's get wrangling and plotting...
library(tidyverse)

# Fixing parameter names
fix.names <- function(name.vec){
  # D // D_{k,ef} -> D[k,ef]
  # beta // beta_{kp} -> beta[kp]
  # sigma // sigma^2_k -> sigma[2]^2
  # gamma // gamma_k -> gamma[k]
  # zeta // zeta_p -> zeta[p]
  sapply(name.vec, function(x){  
    if(grepl("^D\\_\\{", x)){
      bit.in.braces <- str_extract(x, "\\{(.*?)\\}") %>% gsub("\\{|\\}", "", .)
      before.com <- str_extract(bit.in.braces, "^.*?\\,") %>% gsub("\\,", "", .)
      after.com <- str_extract(bit.in.braces, "\\,.*?$") %>% gsub("\\,", "", .)
      param <- str_extract(x, "^.*?\\_") %>% gsub("\\_", "", .)
      return(paste0(param, "[", before.com, '*","*', after.com, ']'))
    }else if(grepl("^beta\\_", x)){
      bit.in.braces <- str_extract(x, "\\{(.*?)\\}") %>% gsub("\\{|\\}", "", .)
      param <- str_extract(x, "^.*?\\_") %>% gsub("\\_", "", .)
      return(paste0(param, "[~~", bit.in.braces, ']'))
    }else if(grepl("^sigma", x)){
      return(paste0("sigma[~~", str_extract(x, "\\d$"), ']^{~~2}'))
    }else{
      return(paste0(str_extract(x, "^.*?\\_") %>% gsub("\\_", "", .),
                    "[~~", str_extract(x, "\\d$"), "]"))
    }
  })
}
all3$Parameter2 <- fix.names(all3$Parameter)
all3$Parameter2 <- gsub('00\\]', '"00"]',all3$Parameter2) # Weird behaviour of expression, turn to string...

# Extracting point estimate
all3$Estimate <- str_extract(all3$MSE %>% trimws, "^.*?\\s") %>% as.numeric

# Separating lower and upper bands
all3b <- separate(all3, CI, c("lb", "ub"), "\\,+", remove = F) %>% 
  mutate_at("lb", ~ gsub("^\\[", "", .x) %>% as.numeric) %>% 
  mutate_at("ub", ~ gsub("\\]$", "", .x) %>% as.numeric) %>% 
  mutate_at("Parameter2", fct_inorder)

# Groupings?
all3b <- all3b %>% 
  mutate(
    grouping = case_when(
      grepl("^D\\[", Parameter2) ~ "Covariance matrix",
      grepl("^beta", Parameter2) ~ "Fixed effects",
      grepl("^sigma", Parameter2) ~ "Residual variance",
      T ~ "Survival parameters"
    )
  )

all3b %>% 
  ggplot(aes(x = fct_rev(Parameter2), y =  Estimate, colour = Method)) + 
  geom_pointrange(aes(ymin = lb, ymax = ub),
                  position = position_dodge(width = 0.33),
                  size = .25, lwd = .25) +
  scale_x_discrete(labels = function(l) parse(text=l)) + 
  theme_csda() +
  coord_flip()


all3b %>% 
  ggplot(aes(x = fct_rev(Parameter2), y =  Estimate, colour = Method)) + 
  geom_pointrange(aes(ymin = lb, ymax = ub),
                  position = position_dodge(width = 0.33),
                  size = .1, lwd = .25) +
  scale_x_discrete("", labels = function(l) parse(text = l)) + 
  facet_wrap(~grouping, scales = 'free') + 
  # scale_colour_brewer(palette = "YlOrRd") + 
  # Shift YlOrRd because hard to see!
  scale_colour_manual(values = c("#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")) +
  theme_csda() +
  theme(
    strip.text = element_text(vjust = 1, size = 7),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=6)
  ) + 
  coord_flip()

ggsave("FinalModelAllMethods.png", width = 140, height = 110, units = "mm")  
