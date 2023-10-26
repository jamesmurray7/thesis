rm(list=ls())
source(".Rprofile")
library(ggplot2)
load(save.dir.file.path(x = "Joint/Multivs/biv2.RData"))
PBCredx$histologic2 <- ifelse(PBCredx$histologic%in%c("3","4"), 1, 0)

# Longit. residuals vs. fitted --------------------------------------------
R.long <- residuals(joint.biv2,'l','p')

alpha.names <- joint.biv2$ModelInfo$Resps
out <- setNames(vector("list", length(alpha.names)), alpha.names)
for(f in alpha.names){
  resp <- ifelse(f == "serBilir", "Serum bilirubin", tools::toTitleCase(f))
  
  # Make ggplots for better automatic scaling...
  R <- R.long[[f]]
  x <- attr(R, "fitted")
  out[[f]] <- data.frame(response = resp, fitted = x, pr = R, row.names = NULL)
}
out.all <- do.call(rbind, out)
out.all$response <- factor(out.all$response, c("Serum bilirubin", "Albumin"))

ggplot(out.all, aes(x = fitted, y = pr)) + 
  geom_hline(aes(yintercept = 0), lty = 5, col = 'red', lwd = 0.25) + 
  geom_point(size=.005, alpha = .2) + 
  facet_wrap(~ response, nrow = 1, scales = 'free') + 
  labs(y = "Pearson residuals", x = "Fitted") + 
  theme_csda()+ 
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    axis.title = element_text(size=5),
    strip.text = element_text(size=7,vjust=1)
  )

ggsave("FinalPearsonResiduals.png", width = 140, height = 60, units = "mm")


# Cox-Snell residuals -----------------------------------------------------
R.surv <- residuals(joint.biv2, 's')

png("FinalCoxSnellExamples.png", width = 140, height = 120, units = "mm", res = 5e2)
sf <- survfit(Surv(R.surv, unlist(joint.biv2$dmats$ph$Delta)) ~ 1)
plot(sf, mark.time = F, lwd = 0.75,
     xlab = "Cox-Snell residuals", ylab = "Survival probability")
curve(exp(-x), from = 0, to = max(sf$time), add = T, lwd = 1, col = 'steelblue')
dev.off()


# Conditional distribution of random effects ------------------------------
CR <- cond.ranefs(joint.biv2, tune = 1.1, burnin = 2000, N = 10000)
# This actually isn't interesting at all!
# png("FinalCR.png", width = 140, height = 140, units = "mm", res = 5e2)
# plot(CR)
# dev.off()
# Much better to show e.g. one subject
hist(CR$acceptance) # Looks good.
plot(CR, 6)


# LRT with joint.biv ------------------------------------------------------
load(save.dir.file.path(x = "Joint/Multivs/biv.RData"))

a <- anova(joint.biv2, joint.biv) # No evidence at 5% level to include sex in model.
                             # But there is some at 10% level.
a$LRT
a$p

# Example dynamic predictions ---------------------------------------------
# Load `joineRML` for comparison's sake?
load(save.dir.file.path("Joint/Multivs/biv2jML.RData"))
# And `JMbayes2` ?
load(save.dir.file.path("Joint/JMb2/Biv2.RData"))

# We want to make these ggplot-ified
ggdynpred <- function(d, actual.fail){
  closest.u.to.fail <- which.min(abs(d$pi$u-actual.fail))
  closest <- d$pi[closest.u.to.fail,,drop=F]
  if(actual.fail > max(d$pi$u)) col.x = "#4c98fe" else col.x = .nice.orange # for plotting those who dont fail in window
  ggplot(as.data.frame(d$pi), aes(x = u, y = mean)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = .nice.orange, alpha = .22, colour = 'grey90', lwd = .5) + 
    geom_step(lwd=.33) + 
    # geom_step(aes(y = lower), lty = 5, lwd = .25) + 
    # geom_step(aes(y = upper), lty = 5, lwd = .25) + 
    geom_point(data = closest, aes(x = u, y = mean), colour = col.x, pch = 4) + 
    expand_limits(y = c(0,1)) + 
    labs(x = expression(u),
         y = expression(Pr*"("*T[i]^"*">=u*"|"*T[i]^"*"*">"*t*","~D[i](t)*","~D*";"~Omega~" )")) + 
    theme_csda() + 
    theme(
      axis.title.y = element_text(size = 6),
      axis.text = element_text(size = 6)
    )
}

# Find example profiles
q.fup <- quantile(with(PBCredx, tapply(time, id, function(i) length(unique(i)))))
# Almost immediate failure ------
unname(which(sapply(unique(PBCredx$id), function(i){
  idat <- PBCredx[PBCredx$id==i,,drop=F]
  setNames(idat$status[1] == 1L & length(idat$time) <= q.fup[2], i)
})))

# id 1 --> Almost immediate failure
d1 <- dynPred(PBCredx, 1, joint.biv2)
g1 <- ggdynpred(d1,PBCredx[PBCredx$id == 1, 'survtime'][1])
g1 + coord_cartesian(xlim = c(d1$pi$u[1]-1e-2,5))
ggsave("./immediatefail.png", width = 140, height = 60, units = "mm")

# d1.jML <- joineRML::dynSurv(joint.biv2.joineRML, PBCredx[PBCredx$id == 1,]) # Doesnt work!

# Survive _almost_ full follow-up ----
to.check <- unname(which(sapply(unique(PBCredx$id), function(i){
  idat <- PBCredx[PBCredx$id==i,,drop=F]
  setNames(idat$status[1] == 1L & length(idat$time) > q.fup[4], i)
})))

joint.biv2$dmats$surv$ft[joint.biv2$dmats$surv$ft>q.fup[4]] # 16 possible times --> pick out prof with all 16 upcoming.
with(PBCredx[PBCredx$id%in%to.check,], tapply(time,id,tail,1))

# id: 21 from this --> 
d2 <- dynPred(PBCredx, 21, joint.biv2)
g2 <- ggdynpred(d2, PBCredx[PBCredx$id == 21, 'survtime'][1])
g2
ggsave("./laterfail.png", width = 140, height = 60, units = "mm")

# Censored during middle quantile of follow-up ----
to.check <- unname(which(sapply(unique(PBCredx$id), function(i){
  idat <- PBCredx[PBCredx$id == i,,drop=F]
  nr <- nrow(idat)
  setNames(idat$status[1] == 0L & nr > q.fup[2] & nr <= q.fup[4],
           i)
})))

sort(abs(4-with(PBCredx[PBCredx$id %in% to.check,],tapply(survtime,id,unique)))) # early censor -> 312
sort(abs(7-with(PBCredx[PBCredx$id %in% to.check,],tapply(survtime,id,unique)))) # early censor -> 312


PBCredx[PBCredx$id==312,]
d3 <- dynPred(PBCredx,312,joint.biv2)
d3b <- dynPred(PBCredx,5,joint.biv2)
d3c <- dynPred(PBCredx,7,joint.biv2)

g3 <- ggdynpred(d3, PBCredx[PBCredx$id==312,'survtime'][1])
g3b <- ggdynpred(d3b, PBCredx[PBCredx$id==5,'survtime'][1]) + coord_cartesian(xlim = c(d3b$pi$u[1]-1e-2,10))
g3c <- ggdynpred(d3c, PBCredx[PBCredx$id==7,'survtime'][1])

ggpubr::ggarrange(g3b, g3c, g3, nrow = 1) # 5 (transplanted) -> 7 (alive, late dropout) -> 312 (alive, early dropout).
ggsave("./censors.png", width = 140, height = 60, units = 'mm')

# Wrapper for JMbayes2 predictions
predJMb2 <- function(id, plot = FALSE, object = NULL){
  i.dat <- PBCredx[PBCredx$id == id, , drop = F]
  t0 <- i.dat$time[length(i.dat$time)]
  i.dat$status <- 0
  i.dat$survtime <- t0
  u <- joint.biv2$dmats$surv$ft[joint.biv2$dmats$surv$ft>t0]
  E <- JMbayes2:::predict.jm(J3, newdata = i.dat, times = u, process = 'event',
                             type_pred = 'link', return_newdata = F, n_mcmc = 200L, n_samples = 1000L)
  out <- data.frame(u = E$times,
                    mean = 1 - E$pred,
                    lower = 1 - E$upp,
                    upper = 1 - E$low,
                    id = E$id)
  if(plot){
    if(is.null(object)) stop("Provide correct dynpreds thing from AEM")
    d.mine <- object$pi[,-3]
    both <- rbind(out[-1,], d.mine)
    both$Method <- c(rep("JMbayes2", nrow(both)/2), rep("Approximate EM", nrow(both)/2)) # half n half
    return(
      ggplot(both, aes(x = u, y = mean)) + 
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = .nice.orange, alpha = .22, colour = 'grey90', lwd = .5) + 
        geom_step(lwd=.33) + 
        # geom_step(aes(y = lower), lty = 5, lwd = .25) + 
        # geom_step(aes(y = upper), lty = 5, lwd = .25) + 
        #geom_point(data = closest, aes(x = u, y = mean), colour = .nice.orange, pch = 4) + 
        expand_limits(y = c(0,1)) + 
        labs(x = expression(u),
             y = expression(Pr*"("*T[i]^"*">=u*"|"*T[i]^"*"*">"*t*","~D[i](t)*","~D*";"~Omega~" )")) + 
        theme_csda() + 
        theme(
          axis.title.y = element_text(size = 6),
          axis.text = element_text(size = 6)
        ) + 
        facet_wrap(~Method)
    )
  }else{
    return(out)
  }
}
d1.JMb <- predJMb2(1)
predJMb2(1, T, d1)
predJMb2(7, T, d3c)



# ROC ---------------------------------------------------------------------

those.who.fail <- PBCredx[PBCredx$status == 1L, c("id", "survtime")]
num.who.fail <- length(unique(those.who.fail$id))
times <- (with(those.who.fail, tapply(survtime,id,unique)))
quantile(times[times>=2],)
sum(times>2&times<3.5)/length(times)
sum(times>3.5&times<7)/length(times)
sum(times>7&times<15)/length(times)

windows <- list("w1"=c(2,3.5),"w2"=c(3.5,7),"w3"=c(7,14)) # These are windows, not {T, delta} notation!!!!

# Number of failures
num.fails <- sapply(windows, function(w){
  t <- w[1]; h <- w[2]
  sum(t < times & times <= h)
})

rates <- lapply(windows, function(w){
  t <- w[1]; h <- w[2]
  tt <- times[t < times & times <= h]
  diff(sort(tt))
})

ROCs <- lapply(windows, function(w){
  Tstart <- w[1]
  delta <- w[2]-Tstart
  ROC(joint.biv2, PBCredx, Tstart, delta, control = list(nsim = 0))
})

ROCs.MC <- lapply(windows, function(w){
  Tstart <- w[1]
  delta <- w[2]-Tstart
  ROC(joint.biv2, PBCredx, Tstart, delta, control = list(nsim = 100))
})

sapply(ROCs, function(x) c(x$PE, x$AUC)) # both PE and AUC get worse together, which is good!

save(ROCs, file = save.dir.file.path("ROC_windows.RData"))
save(ROCs.MC, file = save.dir.file.path("ROC_MC__windows.RData"))

sapply(ROCs.MC, function(x) c(x$PE, x$AUC)) # both PE and AUC get worse together, which is good!

# Three-in-one ROC plot
png("./FinalROCs.png", width = 140, height = 60, units = 'mm', res = 5e2, pointsize = 7)
par(mfrow=c(1,3))
lapply(ROCs, plot, legend = FALSE)
par(mfrow=c(1,1))
dev.off()

## can we ggplot this?
roc2ggplot <- function(r){ # cant be bothered + step behaving weird?
  m <- r$metrics
  ggplot(m, aes(x = FPR, y = TPR)) + 
    geom_abline(slope = 1, intercept = 0, lty = 5) + 
    geom_step() + 
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = expression("ROC curve for "*w[1]))
}

# Predictive error with albumin removed -----------------------------------

serBilir.only <- joint(
  long.formulas = list(serBilir ~  ns(time, knots = c(1, 4)) + histologic2 + 
                         sex + (1 + ns(time, knots = c(1, 4)) | id)
                       ),
  surv.formula = Surv(survtime, status) ~ age + histologic2,
  data = PBCredx, family = list("gaussian"), control = list(tol.rel = 5e-3)
)

ROCs.sbo <- lapply(windows, function(w){
  Tstart <- w[1]
  delta <- w[2]-Tstart
  ROC(serBilir.only, PBCredx, Tstart, delta, control = list(nsim = 0))
})


# DynPred examples --------------------------------------------------------
# Window 1
View(PBCredx[PBCredx$survtime > windows[[1]][1] & PBCredx$survtime < windows[[1]][2],])
# False negative
a <- dynPred(PBCredx, 14, joint.biv2, u = ROCs[[1]]$candidate.u)
false.neg <- ggdynpred(a, PBCredx[PBCredx$id==14,'survtime'][1]) + geom_hline(yintercept=0.91, colour = 'grey80', lty = 5)
# True positive
true.pos <- ggdynpred(dynPred(PBCredx,69,joint.biv2,ROCs[[1]]$candidate.u), 
                      PBCredx[PBCredx$id==69,'survtime'][1])+geom_hline(yintercept=.91, colour = 'grey80', lty = 5)
# False positive
View(PBCredx[PBCredx$survtime>windows[[1]][2],])
false.pos <- ggdynpred(dynPred(PBCredx,55,joint.biv2,ROCs[[1]]$candidate.u), 
                       PBCredx[PBCredx$id==55,'survtime'][1])+geom_hline(yintercept=.91, colour = 'grey80', lty = 5) 
# True negative
true.neg <- ggdynpred(dynPred(PBCredx,2,joint.biv2,ROCs[[1]]$candidate.u), 
                      PBCredx[PBCredx$id==2,'survtime'][1])+geom_hline(yintercept=.91, colour = 'grey80', lty = 5)

ggpubr::ggarrange(true.pos, true.neg, false.pos, false.neg, nrow = 2, ncol = 2)
ggsave("./four_example_dynpreds.png", 
       width = 140, height = 90, units = "mm")
