PBC2 <- PBC %>% 
  select(id, drug, survtime, status, time, albumin, serBilir, platelets, alkaline) %>% 
  na.omit(.) %>% 
  mutate_at("serBilir", log) %>% 
  tidyr::pivot_longer(serBilir:alkaline, names_to = 'biomarker', values_to = "value") %>% 
  mutate_at("drug", ~ifelse(.x == 1, "Treatment", "Placebo")) %>% 
  mutate_at("drug", forcats::fct_inorder) %>% 
  mutate_at("status", ~ifelse(.x == 1, "Died", "Survived")) %>% 
  mutate(
  biom = case_when(
    biomarker == "alkaline" ~ "Alkaline phosphotase",
    biomarker == "platelets" ~ "Platelet count",
    biomarker == "serBilir" ~ "log(Serum bilirubin)",
    T ~ "AAA"
  ),
  fam = case_when(
    biomarker == "alkaline" ~ "Negative binomial",
    biomarker == "platelets" ~ "Poisson",
    biomarker == "serBilir" ~ "Gaussian",
    T ~ "AAA"
  )) %>% 
  mutate(
    biom = factor(biom, c("log(Serum bilirubin)", "Platelet count", "Alkaline phosphotase"))
  )

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE, digits = 2)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  l <- gsub("10\\^\\+", "10^", l)
  l <- ifelse(grepl("10\\^00", l), "0", l)
  ind1 <- which(grepl("1\\.0", l))
  if(length(ind1))
    l[1:ind1] <- gsub("\\.0", "", l[1:ind1])
  #print(l)
  # return this as an expression
  parse(text=l)
}

plotter <- function(resp){
  .title <- switch(resp,
                serBilir = "Gaussian",
                albumin = "Gaussian",
                platelets = "Poisson",
                alkaline = "Negative binomial")
  if(resp == "alkaline") PBC2 <- PBC2 %>% dplyr::filter(value < 8000)
  P <- PBC2 %>% 
    dplyr::filter(biomarker == resp) %>% 
    ggplot(aes(x = time, y = value, group = id)) + 
    geom_line(lwd = 0.25, alpha = 0.25) + 
    geom_smooth(aes(group = status, colour = status), se = FALSE, formula = y~x, method = "loess") + 
    facet_wrap(~biom, nrow = 1, scales = "free_y", strip.position = "left") + 
    expand_limits(x = c(0,15.5)) + 
    scale_x_continuous(breaks = seq(0,15,5), labels = seq(0,15,5) %>% as.character) +
    theme_csda() + 
    labs(y = NULL, colour = NULL, x = NULL, title = .title) + 
    scale_colour_manual(values = c(nclred, "#fd8d3c")) + 
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 21, colour = nclblue, vjust = 1),
      legend.position = "none",
      plot.title = element_text(size = 21, colour = nclblue, vjust = 0.25),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 15)
    )
  
  if(resp == "alkaline"){
    P <- P + 
      scale_y_continuous(labels = fancy_scientific,
                         breaks = scales::pretty_breaks(n = 5))
  }
  if(resp == "serBilir") P <- P + theme(legend.position = c(0.7,0.1), legend.key = element_blank())
  P
}

P1 <- plotter("serBilir")  
P2 <- plotter("platelets")
P3 <- plotter("alkaline")

P <- gridExtra::grid.arrange(P1,P2,P3, ncol = 3, nrow = 1,
                             bottom = grid::textGrob(label = "Follow-up time (years)",
                                                     gp = grid::gpar(fontsize = 19)))


png("PBCapp.png", width = linewidth, height = 130, units = "mm", res = 1500)
grid::grid.draw(P)
dev.off()
