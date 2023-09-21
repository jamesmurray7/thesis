source('../theme_csda.R')
theta <- c(-1, 0.0)



getOUT <- function(n, family, D){ # Wrapper for simulation + Sampling
  d <- .sim(family, n, D)
  btrue <- d[[2]]; D <- d[[3]]; data <- d[[1]]
  OUT <- Sample(data, btrue, family, 1:n, D)
}

plotOut <- function(OUT, save.dir = '/data/c0061461/THESIS/Gaussian-Normal-Justification/', num.to.show = 5L, id = NULL){
  fn <- paste0(save.dir, OUT$family, "_RE_mi_breakdown",id,".png")
  pp <- quantile(OUT$Acc, probs = c(0.025, 0.500, 0.975))
  cat(sprintf("Median [2.5%%, 97.5%%] acceptance rate for %s: %.2f [%.2f, %.2f]\n\n",
              OUT$family, pp[2], pp[1], pp[3]))
  df <- OUT$df
  
  # Sort into groups / make labels
  mi <- OUT$mi
  grp <- ifelse(mi <= 5, "[1,5]", ifelse(mi <= 10, "[6,10]", ifelse(mi <= 15,"[11,15]", "zzz")))
  grp <- factor(grp, 
                levels = c("[1,5]", "[6,10]", "[11,15]"),
                labels = c('"["*1*","*5*"]"', '"["*6*","*10*"]"', '"["*11*","*15*"]"'))
  
  mi.tab <- table(mi)
  mi.lab <- paste0(names(mi.tab[mi]), "\n(n=", mi.tab[mi],')')
  grp.tab <- table(grp)
  grp.lab <- paste0("m[i] %in% ", names(grp.tab[grp]))
  
  mini <- data.frame(m = mi, mi.grp = grp, mi.lab = mi.lab, mi.grp.lab = grp.lab)
  mini <- mini[!duplicated(mini),]
  mini <- mini[order(mini$m),]
  mini$mi.grp <- forcats::fct_inorder(mini$mi.grp)
  mini$mi.lab <- forcats::fct_inorder(mini$mi.lab)
  mini$mi.grp.lab <- forcats::fct_inorder(mini$mi.grp.lab)
  
  df <- left_join(df, mini, 'm')
  random.ids.to.keep <- unlist(unname(
    with(df, tapply(id, mi.grp.lab, function(x) sample(unique(x), size = num.to.show, replace = FALSE)))
  ))
  df <- df[df$id %in% random.ids.to.keep, ]
  
  # Plot (lots of) densities
  if(OUT$family != "binomial"){
    b0plot <- df %>% filter(var=='b[0]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_wrap(~mi.grp.lab, scales = 'free', labeller = label_parsed, nrow = 3, ncol = 1,
                 strip.position = 'left')+
      labs(y = '', title = bquote(b[0]), x = '') +
      theme_csda() + 
      theme(
        strip.placement = 'outside',
        plot.title.position = 'panel',
        plot.title = element_text(hjust = 0.5, face = 'plain', vjust = 0)
      )
    b1plot <- df %>% filter(var=='b[1]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_wrap(~mi.grp.lab, scales = 'free', labeller = label_parsed, nrow = 3, ncol = 1,
                 strip.position = 'left')+
      labs(y = '', x = '', title = bquote(b[1])) + 
      theme_csda()+
      theme(
        strip.text = element_blank(),
        plot.title.position = 'panel',
        plot.title = element_text(hjust = 0.5, face = 'plain', vjust = 0)
      )
    png(fn, width = 190, height = 120, units = 'mm', pointsize = 9, res = 1000)
    gridExtra::grid.arrange(b0plot, b1plot, nrow=1, ncol=2)
    dev.off()
  }else{
    b0plot <- df %>% filter(var=='b[0]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_wrap(~mi.grp.lab, scales = 'free', labeller = label_parsed, nrow = 3, ncol = 1)+
      labs(y = bquote(b[0]), x = '') + 
      theme_csda()
    png(fn, width = 190, height = 120, units = 'mm', pointsize = 9, res = 1000)
    print(b0plot)
    dev.off()
  }
  
  cat("Done for", OUT$family, "showing", num.to.show, "posteriors,\n")
  cat("saved in", fn, ".\n\n")
}

getPlot <- function(n, family, num.to.show = 5L, D = NULL, id = NULL){
  OUT <- getOUT(n, family, D)
  plotOut(OUT, id = id)
  rm(OUT)
  on.exit(gc())
}


getPlot(60, "poisson")
theta <- c(-1, 0.1)
getPlot(60, "binomial")


# Gaussian ----------------------------------------------------------------
# getPlot(60, "gaussian")
getPlot(60, "gaussian", D = diag(c(.25*3, .09*3),2), id = "diagD3x")
getPlot(100, "gaussian", D = diag(c(.25*10, .09*10),2), id = "diagD10x")
getPlot(100, "gaussian", D = matrix(c(2.5, .125, .125, 0.9), 2, 2), id = "diagD10xCorr")
getPlot(100, "gaussian", D = matrix(c(5, 0, 0, 1), 2, 2), id = "large")
