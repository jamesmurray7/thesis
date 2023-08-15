library(gmvjoint)
library(tidyverse)
source('../theme_csda.R')
data("PBC", package = "gmvjoint")

dat <- PBC %>% 
  select(id, survtime, status, drug, age, sex, histologic, # Baseline
         time,
         serBilir:prothrombin,                             # Longit. biomarkers (cont/count)
         ascites:spiders) %>%                              # Longit. biomarkers (binary)
  select(-serChol) %>% 
  mutate(
    id = as.numeric(as.character(id)),
    drug = factor(drug, 0:1),
    serBilir = log(serBilir),
    AST = log(SGOT)
  ) %>% 
  group_by(id) %>% 
  mutate(
    rr = row_number()
  ) %>% 
  ungroup %>% 
  mutate(
    hist_bl = ifelse(rr == 1, histologic, NA)
  ) %>% select(-rr, -histologic, -SGOT) %>% 
  fill(hist_bl) %>% 
  mutate_at("hist_bl", ~factor(.x, 1:4)) %>% as.data.frame


# Longitudinal trajectories for the six markers ---------------------------
#                                    (continuous)
dat.long <- dat %>% 
  select(-c(ascites, hepatomegaly, spiders)) %>% 
  pivot_longer(cols = c(serBilir:prothrombin, AST),
               names_to = "biomarker", values_to = "measurement") %>% 
  mutate(
    biom.lab = case_when(
      biomarker == "serBilir" ~ "log(Serum~bilirubin)",
      biomarker == "albumin" ~ "Serum~albumin",
      biomarker == "alkaline" ~ "Alkaline~phosphotase",
      biomarker == "platelets" ~ "Platelet~count",
      biomarker == "prothrombin" ~ "Prothrombin~time",
      biomarker == "AST" ~ "log(AST)",
      T ~ "999"
    )
  ) %>% 
  mutate_at("biom.lab", forcats::fct_inorder) %>% 
  mutate_at("status", ~factor(.x, 0:1, c("Survived", "Died"))) 

# Make sure font matches rest of thesis later 15/8/24 =+=+=+=+=+=+=
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  if(max(l, na.rm = T) < 1.5e3) return(l) # SO only alkaline is affected
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

dat.long %>% 
  filter(measurement <= 8e3) %>% # Presentation purposes on Alkaline
  ggplot(aes(x = time, y = measurement, group = id)) + 
  geom_line(alpha = .266, lwd = .125) + 
  facet_wrap(~biom.lab, scales = 'free_y', labeller = label_parsed,
             strip.position = "left") + 
  labs(y = NULL, x = "Follow-up time (years)", colour = "") + 
  geom_smooth(aes(group = status, colour = status), 
              formula = y~x, method = "loess", se = FALSE) + 
  scale_color_manual(values = rev(c("#e31836","#fd8d3c"))) + 
  theme_csda("") + 
  expand_limits(x = c(0,15)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(6),
                     labels = fancy_scientific) +
  theme(
    strip.placement = "outside",
    legend.position = c(.94,.975),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.text = element_text(size = 6),
    legend.text = element_text(size = 5),
    axis.text = element_text(size = 4.5),
    axis.title = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    legend.spacing = unit(2, "mm")
  )

png("./PBC_Cha5.png", width = 140, height = 90, units = "mm", res = 1e3)
last_plot()
dev.off()

# Tile/heatmap for binary ones? -------------------------------------------

# What do we want to show?
# For each of the three binary markers, the proportion of total sample 
# with presence at given time, coloured somehow by failure?
bin.long <- dat %>% 
  select(id, survtime, status, time, ascites, hepatomegaly, spiders) %>% 
  pivot_longer(cols = c(ascites, hepatomegaly, spiders),
               names_to = "biomarker", values_to = "measurement") %>% 
  filter(!is.na(measurement))

cuts <- unique(quantile(sort(bin.long$time), seq(0.0, 1.0, length.out = 20)))
bin.long$time2 <- cut(bin.long$time, cuts, include.lowest = T)

bin.long %>% 
  # group_by(time, biomarker, status) %>% 
  group_by(biomarker, time2) %>% 
  # slice(1:250) %>% 
  mutate(
    n = n()
  ) %>% 
  mutate(
    present.pc = 100 * sum(measurement)/n
  ) %>% ungroup %>% 
  mutate_at("biomarker", stringr::str_to_sentence) %>% 
  ggplot(aes(x = time2, y = biomarker))+ 
  geom_tile(aes(fill = present.pc), colour = "black") + 
  theme_csda() + 
  labs(
    y = NULL, x = NULL,# x = "Follow-up time (years)",
    fill = "% Present in time window"
  ) + 
  scale_fill_continuous(low = "purple4", high = "#fd8d3c") + 
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 8, angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 6, angle = 360-45, vjust=-3),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(angle = 270, vjust = 0.25, size=7.5),
    legend.text = element_text(angle=270, hjust = 0.4, size=6.5),
    legend.key.width = unit(3, "mm")
  )


png("./binheatmap.png", width = 140, height = 90, units = "mm", res = 1e3)
last_plot()
dev.off()
