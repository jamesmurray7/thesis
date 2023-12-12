out.list <- setNames(vector("list", 6), .valid.families)
for(f in .valid.families){
  cat(f,'\n')
  assign("f.dat", get(load(paste0(save.dir.file.path(family.dir.name(f)), "/psimi.RData"))))
  # Go through and get new plots outputted !!!!!!
  out <- paste0(save.dir.file.path(family.dir.name(f)), "/ellipse-",f,"2.png")
  # print(to.eval)
  png(file = out, width = 140, height = 110, units = "mm", res = 500)
  plot(f.dat[[2]])
  dev.off()
  out.list[[f]] <- f.dat[[1]]
  cat(sprintf("----%s done; saved in %s----\n", f, out))
}

fams <- as.data.frame(do.call(rbind, out.list))

# This plot ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
fams %>%
  mutate(famname2 = purrr::map_chr(family, family.dir.name)) %>% 
  mutate_at("famname2", ~gsub("\\-Normal\\-Justification$", "", .x)) %>%
  mutate_at("famname2", ~gsub("isedpoi", "ised Poi", .x)) %>% 
  mutate_at("famname2", ~gsub("ivebin", "ive bin", .x)) %>% 
  filter(famname2!='Binomial') %>% 
  ggplot(aes(x = as.factor(mi), y = psi, fill = as.factor(surv))) + 
  geom_boxplot(outlier.alpha = 0.33, lwd = 0.25, 
               fatten = 2, outlier.size = 0.50) + 
  facet_wrap(~famname2, scales = 'free_y') + 
  scale_fill_manual(values = c("#FFEDA0", "#F03B20")) + 
  labs(x = expression(m[i]), y = expression(psi["  "*i])) + 
  theme_csda() + 
  theme(
    legend.position = "none",
    strip.text = element_text(vjust = 1)
  )
# Save...
ggsave(save.dir.file.path("psi-mi-families2.png"), width = 140, height = 90, units = "mm")
