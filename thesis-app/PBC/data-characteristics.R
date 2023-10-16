rm(list=ls())
library(tidyverse)
library(xtable)
data("pbc2", package = "joineRML")

# Copy from ../../thesis-app/PBC/.Rprofile
# Create some extra datasets 
PBCredx <- na.omit(pbc2[,c("id", "years", "status2", "drug", "sex",
                          "age", "histologic", "year", "hepatomegaly", "spiders",
                          "alkaline", "SGOT", "platelets", "prothrombin", 
                          "serBilir", "albumin")])

BL <- PBCredx[!duplicated(PBCredx$id), c("id", "histologic")]
BL$histologic <- factor(BL$histologic, 1:4)

PBCredx <- merge(dplyr::select(PBCredx, -histologic), 
                 BL, "id")
surv.data <- PBCredx[!duplicated(PBCredx[,'id']), ]

PBCredx <- PBCredx %>% 
  select(id, years, drug, sex, status2, age, histologic) %>% 
  group_by(id) %>% 
  mutate(r = n()) %>% 
  ungroup %>% 
  distinct

# All but histologic & sex
conts <- PBCredx %>% 
  group_by(status2) %>% 
  summarise(.groups = "keep",
    n = n(),
    med.fup = median(r),
    p25.fup = quantile(r, .25),
    p75.fup = quantile(r, .75),
    mn = mean(age), sd = sd(age),
  ) %>% t

with(surv.data, t.test(age ~ status2))
with(PBCredx, wilcox.test(r ~ status2))

PBCredx %>% count(status2, sex) %>% group_by(status2) %>% mutate(p = prop.table(n) * 100)
surv.data %>% count(status2, histologic) %>% group_by(status2) %>% mutate(p = prop.table(n)*100)
surv.data %>% count(status2, drug) %>% group_by(status2) %>% mutate(p = prop.table(n)*100)

chisq.test(table(surv.data$sex,surv.data$status2))
chisq.test(table(surv.data$hist,surv.data$status2))
chisq.test(table(surv.data$drug,surv.data$status2))
