rm(list=ls())
source(".Rprofile")

in.dir <- save.dir.file.path("fits")
(files <- dir(in.dir, pattern = "\\.RData"))
nms <- sapply(files, function(f){
  r <- as.numeric(gsub("\\.RData|^gh", "", f))
  r <- 
    if(r == 97)
      "rho(10)" # I didn't know you could do this in R!
    else if(r == 98)
      "rho(25)"
    else if(r == 99)
      "rho(33)"
    else
      paste0("rho = ", r)
  r
})

.loader <- function(file){
  assign("out", get(load(save.dir.file.path(file, in.dir))))
}

parsed <- setNames(lapply(files, function(f){
  fits <- .loader(f)
  this.parsed <- lapply(fits, extract.from.joint)
  Omega <- sapply(this.parsed, function(y) y$Omega)
  SE <- sapply(this.parsed, function(y) y$SE)
  lb <- sapply(this.parsed, function(y) y$lb)
  ub <- sapply(this.parsed, function(y) y$ub)
  et <- unname(sapply(this.parsed, function(y) y$elapsed))
  tt <- unname(sapply(this.parsed, function(y) y$total))
  gh <- unname(sapply(this.parsed, function(y) y$gh.nodes))
  cat(sprintf("\n%s done!\n", f))
  list(Omega = Omega, SE = SE, lb = lb, ub = ub, elapsed = et, total = tt,
       gh = gh)
}), nms)

TM <- t(create.targetmat(N = 200L)) # One for all N sets

step <- lapply(seq_along(parsed), function(i){
  this.nm <- nms[i];
  this.parsed <- parsed[[i]]
  Omega <- this.parsed$Omega; rn <- row.names(Omega)
  to.keep <- c(grep("^Y\\.3", rn), grep("^sigma", rn), grep("gamma\\_1", rn):length(rn))
  Omega2 <- Omega[to.keep,]
  ubs <- this.parsed$ub[to.keep,]; lbs <- this.parsed$lb[to.keep,]
  # Shorten target mat (which has different number of entries, annoyingly so cant use `to.keep`).
  TM2 <- TM[grep("\\\\beta\\_\\{30\\}", row.names(TM)):length(row.names(TM)),]
  
  # Summary measures -- unsure what we need?
  Bias <- Omega2-TM2
  rt.Bias <- sqrt(abs(Bias))
  MSE <- (TM2-Omega2)^2
  row.names(MSE) <- row.names(Bias)
  CP <- rowSums(lbs<TM2&ubs>TM2)/200
  
  # Summaries on these (again unsure what to plot)
  mn.Bias <- rowMeans(Bias)
  md.Bias <- apply(Bias, 1, median)
  q.Bias <- apply(Bias, 1, quantile, c(.025,.975))
  mn.abs.Bias <- rowMeans(abs(Bias))
  md.abs.Bias <- apply(abs(Bias), 1, median)
  q.abs.Bias <- apply(abs(Bias), 1, quantile, c(.025,.975))
  
  df <- data.frame(id = this.nm,
                   Parameter = names(mn.Bias),
                   # Raw bias summaries
                   mnBias = mn.Bias, md.Bias = md.Bias,
                   p025Bias = q.Bias[1,], p975Bias = q.Bias[2,],
                   # Absolute bias summaries
                   mnabsBias = mn.abs.Bias, mdabsBias = md.abs.Bias,
                   p025absBias = q.abs.Bias[1,], p975absBias = q.abs.Bias[2,],
                   row.names = NULL, stringsAsFactors = FALSE
                   )
  
  return(list(nm = this.nm, df = df,
         Bias = Bias, rt.Bias = rt.Bias, MSE = MSE, CP = CP))

})


# Idea 1: caterpillar plot ------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
# Helper functions
t.adf <- function(x) as.data.frame(t(x))
to.long <- function(x) pivot_longer(x, `Y.3_(Intercept)`:`zeta_bin`,
                                    names_to = "Parameter", values_to = "value")
parsed.steps <- lapply(seq_along(step), function(i){
  this <- step[[i]]
  this.nm <- unname(this$nm)
  this.Biases <- t.adf(this$Bias)
  this.Biases$measure <- "Bias"; this.Biases$i <- 1:nrow(this.Biases)
  this.rt.Bias <- t.adf(this$rt.Bias)
  this.rt.Bias$measure <- "sqrtBias"; this.rt.Bias$i <- 1:nrow(this.rt.Bias)
  this.MSE <- t.adf(this$MSE)
  this.MSE$measure <- "MSE"; this.MSE$i <- 1:nrow(this.MSE)
  
  out <- rbind(to.long(this.Biases), 
               to.long(this.rt.Bias), 
               to.long(this.MSE))
  out$id <- this.nm
  out
})

parsed.steps <- do.call(rbind, parsed.steps)

parsed.steps2 <- parsed.steps %>% 
  mutate(
    param.lab = case_when(
      grepl("gamma\\_\\d", Parameter) ~ paste0("gamma[", stringr::str_extract(Parameter, "\\d"), "]"),
      grepl("zeta", Parameter) ~ "zeta",
      grepl("sigma", Parameter) ~ "sigma[epsilon]^2",
      Parameter == "Y.3_(Intercept)" ~ "beta[30]",
      Parameter == "Y.3_time" ~ "beta[31]",
      Parameter == "Y.3_cont" ~ "beta[32]",
      Parameter == "Y.3_bin" ~ "beta[33]",
    ),
    id = gsub("\\=", "==", id)
  )


parsed.steps2 %>% 
  filter(measure == "Bias", Parameter == "sigma^2_1") %>% 
  ggplot(aes(x = i, y = abs(value), group = id, colour = id)) + 
  geom_line(lwd=.2)

parsed.steps2 %>% 
  filter(measure == "Bias") %>% 
  ggplot(aes(x = id, y = (value), group = id, colour = id)) + 
  # geom_line(lwd=.2)
  geom_boxplot() + 
  facet_wrap(~param.lab, labeller = label_parsed, scales = "free_y")


# Idea 2: Median and quantiles of |Bias| ----------------------------------

summ.tab <- do.call(rbind, lapply(step, '[[', "df"))
summ.tab$gh <- ifelse(grepl("\\=", summ.tab$id), as.numeric(stringr::str_extract(summ.tab$id, "\\d?\\d")),
                      100+as.numeric(stringr::str_extract(summ.tab$id, "\\d?\\d")))
summ.tab$id <- gsub("\\=","==",summ.tab$id)
# Sort out x-axis
summ.tab$x <- ifelse(summ.tab$gh>100, paste0("rho(",summ.tab$gh-100,")"),as.character(summ.tab$gh))
summ.tab2 <- arrange(summ.tab, gh)

fmtter <- function(l, digits = 5){
  l <- scales::scientific(l, digits = digits)
  l <- gsub("e","%*%10",l)
  l <- gsub("-","^{-", l)
  # print(l)
  l <- sapply(l, function(a){
    if(!is.na(a)) paste0(a, "}") else a
  })
  # print(l)
  parse(text=l)
}

summ.tab2 %>% 
  mutate(
    param.lab = case_when(
      grepl("gamma\\_\\d", Parameter) ~ paste0("gamma[", stringr::str_extract(Parameter, "\\d"), "]"),
      grepl("zeta", Parameter) ~ "zeta",
      grepl("sigma", Parameter) ~ "sigma[epsilon]^2",
      Parameter == "Y.3_(Intercept)" ~ "beta[30]",
      Parameter == "Y.3_time" ~ "beta[31]",
      Parameter == "Y.3_cont" ~ "beta[32]",
      Parameter == "Y.3_bin" ~ "beta[33]",
    )
  ) %>% 
  ggplot(aes(x = forcats::fct_inorder(x), y = mdabsBias, group=Parameter)) + 
  geom_line(lwd=.5) + geom_point(size=.85)+
  # geom_line(aes(y=p025absBias), lty = 5) +
  # geom_line(aes(y=p975absBias), lty = 5) + 
  facet_wrap(~param.lab, scales = "free", labeller = label_parsed) + 
  scale_x_discrete(labels = function(x){
    sapply(x, function(i){
      print(i)
      if(!grepl("\\(", i)){ 
        return(as.character(i))
      }else{
        print(i)
        P <- stringr::str_extract(i,"\\d?\\d")
        print(P)
        a <- bquote(rho*"("*.(P)*")")
        print(a)
        return(a)
      }
    })
  }, name=bquote(rho)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(5), labels = function(l) fmtter(l)) +
  labs(y="Median absolute bias") + 
  theme_csda() + 
  theme(
    # axis.ticks.x=element_blank(),
    axis.text.y = element_text(size=4.5),
    axis.text.x = element_text(size=4.05),#, angle = 360-45)
    axis.title = element_text(size=6),
    strip.text = element_text(size=7)
  )

ggsave(save.dir.file.path("GHinvestigation.png", in.dir),
       width = 140, height = 90, units = "mm")

