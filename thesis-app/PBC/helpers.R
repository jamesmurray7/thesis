library(ggplot2)
plotMultiJM <- function(x, sigma.limit = 50){
  stopifnot(inherits(x, "joint"))
  S <- summary(x)
  # Survival sub-model plot ----
  survs <- as.data.frame(S$Surv)[,c("Estimate", "2.5%", "97.5%")]
  survs <- cbind(Parameter = row.names(survs), survs)
  row.names(survs) <- NULL
  # Separate-out {gamma, zeta}
  zetas <- survs[grepl("^zeta\\_", survs$Parameter),]
  gammas <- survs[grepl("^gamma\\_", survs$Parameter),]
  # Clean-up parameter names
  gammas$ParameterLabels <- paste0("gamma[", 1:nrow(gammas), "]")
  gammas$ParameterLabels <- forcats::fct_inorder(gammas$ParameterLabels)
  zetas$ParameterLabels <- paste0("zeta[", 1:nrow(zetas), "]")
  zetas$ParameterLabels <- forcats::fct_inorder(zetas$ParameterLabels)
  # Make survival sub-model plot
  survs <- rbind(gammas, zetas)
  survs$isgamma <- grepl("^gamma\\_", survs$Parameter)
  PlotSurv <- ggplot(data = survs, aes(x = forcats::fct_rev(ParameterLabels), y = Estimate, colour = isgamma)) + 
    geom_hline(aes(yintercept=0), colour = "grey", lty = 5) + 
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) + 
    geom_point() + 
    scale_x_discrete("", labels = function(x) parse(text=x)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
    scale_colour_manual(values = c("darkblue", "red3")) +
    theme_csda() + 
    labs(title = "Survival sub-model")+
    theme(
      legend.position = 'none',
      axis.text.y = element_text(size = 9),
      title = element_text(size = 10)
    ) + 
    coord_flip()
  
  # Longitudinal plots ----
  K <- x$ModelInfo$K
  store <- vector("list", K)
  for(k in 1:K){
    long <- as.data.frame(S$Longits[[k]])[, c("Estimate", "2.5%", "97.5%")]
    long <- cbind(Parameter = row.names(long), long)
    row.names(long) <- NULL
    # beta/sigma split
    betas <- long[grepl(x$ModelInfo$Resps[k], long$Parameter),]
    sigmas <- long[grepl("^sigma|^phi|^shape", long$Parameter), ]
    betas$ParameterLabels <- paste0("beta[",k,'*","*',(1:nrow(betas))-1,']')
    betas$ParameterLabels <- forcats::fct_inorder(betas$ParameterLabels)
    # Work out type of sigma
    if(grepl("sigma", sigmas$Parameter)){
      sigmas$ParameterLabels <- paste0("sigma[",k,"]^2")
    }else{
      sigmas$ParameterLabels <- paste0("phi[",k,'*","*',(1:nrow(sigmas))-1,"]")  
    }
    sigmas$ParameterLabels <- forcats::fct_inorder(sigmas$ParameterLabels)
    sigmas <- sigmas[sigmas$Estimate<sigma.limit,] # Stop daft estimates skewering plot!
    store[[k]] <- rbind(betas,sigmas)
  }
  long <- do.call(rbind, store)
  # Draw horizontal lines between responses
  breaks <- sapply(store, nrow)
  cbreaks <- outer(rev(cumsum(breaks)), rev(cumsum(breaks)), '-')[1,] # Need to reverse.
  cbreaks <- cbreaks[cbreaks>0]
  long$isbeta <- grepl("beta",long$ParameterLabels)
  PlotLong <- ggplot(data = long, aes(x = forcats::fct_rev(ParameterLabels), y = Estimate, colour = isbeta)) + 
    geom_hline(aes(yintercept=0), colour = "grey", lty = 5) + 
    geom_vline(xintercept=c(cbreaks) + 0.5, colour = "black", lty = 3) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) + 
    geom_point() + 
    scale_x_discrete("", labels = function(x) parse(text=x)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(7)) + 
    scale_colour_manual(values = c("darkblue", "red3")) +
    theme_csda() + 
    labs(title = "Longitudinal sub-models")+
    theme(
      legend.position = 'none',
      axis.text.y = element_text(size = 9),
      title = element_text(size = 10)
    ) + 
    coord_flip()
  
  xx <- readline("Press Enter for survival plot >>> ")
  print(PlotSurv)
  xx <- readline("Press Enter for longitudinal plot >>> ")
  print(PlotLong)
  return(list(Surv = PlotSurv, 
              Long = PlotLong))
}
