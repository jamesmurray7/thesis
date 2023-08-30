rm(list=ls())
source('.Rprofile')

# Getting best models from parseLongModels --------------------------------
ff <- readLines("./parseLongModels.R")
header.inds <- grep("->\\s+\\=+", ff)
model.inds <- grep("->>", ff)

# Get into nice readable output format -->
invisible(unname(sapply(diag(outer(header.inds, model.inds, FUN = function(x,y) paste(x,y))), function(x){
  header <- as.numeric(stringr::str_extract(x,"\\d?\\d?\\d?\\s"))
  model <- as.numeric(stringr::str_extract(x,"\\s\\d?\\d?\\d?")  )
  c(header, model)
  cat(ff[header], "\n", ff[model], "\n")
})))

# capture this locally (by rerunning).
capture.output(invisible(unname(sapply(diag(outer(header.inds, model.inds, FUN = function(x,y) paste(x,y))), function(x){
  header <- as.numeric(stringr::str_extract(x,"\\d?\\d?\\d?\\s"))
  model <- as.numeric(stringr::str_extract(x,"\\s\\d?\\d?\\d?")  )
  c(header, model)
  cat(ff[header], "\n", ff[model], "\n")
}))), file = "./printout.txt")


# Fitting with gmvjoint ---------------------------------------------------
head(PBCredx)
library(gmvjoint)
library(splines)

# albumin -> ====
alb <- joint(
  long.formulas = list(albumin ~ age + histologic + time + I(time^2) + 
                         (1 + time + I(time^2)|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian"),
  control = list(tol.rel = 5e-3)
)

plot(resid(alb))
cc <- cond.ranefs(alb); plot(cc) # Very well-fitting (NB as we expect)

# serBilir -> ====
sbil <- joint(
  long.formulas = list(serBilir ~ sex + histologic + ns(time, knots = c(1, 4)) + 
                         (1 + ns(time, knots = c(1, 4))|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian"),
  control = list(tol.rel = 5e-3)
)

plot(resid(sbil))
cc <- cond.ranefs(sbil); plot(cc) # Decent, though intercept is perhaps better suited by skew-normal?

# SGOT -> ====
sgot <- joint(
  long.formulas = list(SGOT ~ age + histologic + drug * (time + I(time^2)) + 
                         (1 + time + I(time^2)|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("gaussian"),
  control = list(tol.rel = 5e-3)
)

plot(resid(sgot))
cc <- cond.ranefs(sgot); plot(cc) # pretty good

# Spiders -> ====
# Seems unstable? :/
# Proceed with just an intercept-only one.
# Might be due to low prevalence at later follow-up times maybe.
spiders <- joint(
  long.formulas = list(spiders ~ age + histologic + time + 
                         (1|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("binomial"),
  control = list(tol.rel = 5e-3, verbose = T)
)

cc <- cond.ranefs(spiders); plot(cc) # slightly left skew but passable; large variance.

# Hepatomegaly -> ====
hepa <- joint(
  long.formulas = list(hepatomegaly ~ sex + histologic + time + I(time^2) +  
                         (1 + time + I(time^2)|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  data = PBCredx,
  family = list("binomial"),
  control = list(tol.rel = 5e-3, verbose = T)
)

plot(resid(hepa, type = 'pearson'))
cc <- cond.ranefs(hepa); plot(cc) # Very good (again).

# Prothrombin -> ====
proth <- joint(
  long.formulas = list(prothrombin ~ age + histologic + drug * (time + I(time^2)) + 
                         (1 + time + I(time^2)|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  disp.formulas = list(~time+age+drug+sex+histologic),
  data = PBCredx,
  family = list("Gamma"),
  control = list(tol.rel = 5e-3, verbose = T)
)

plot(resid(proth, type = 'pearson'))
cc <- cond.ranefs(proth); plot(cc) # Decent.

# Platelets -> ====
plat <- joint(
  long.formulas = list(platelets ~ sex + histologic + drug * ns(time, knots = c(1, 4)) + 
                         (1 + ns(time, knots = c(1, 4))|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  disp.formulas = list(~time + age + sex + histologic),
  data = PBCredx,
  family = list("genpois"),
  control = list(tol.rel = 5e-3, verbose = T)
)

plot(resid(plat, type = 'pearson')) # A little splayed at edges but acceptable.
cc <- cond.ranefs(plat); plot(cc) # Decent.

# Alkaline -> ====
alka <- joint(
  long.formulas = list(alkaline ~ age + histologic + ns(time, knots = c(1, 4)) + 
                         (1 + ns(time, knots = c(1, 4))|id)),
  surv.formula = Surv(survtime, status) ~ age + histologic,
  disp.formulas = list(~time+drug+age+histologic),
  data = PBCredx,
  family = list("negbin"),
  control = list(tol.rel = 5e-3, verbose = T)
)

plot(resid(alka, type = 'pearson')) # A little splayed at edges but acceptable.
cc <- cond.ranefs(alka); plot(cc)   # A little off!

# Loop over univariate joint models and print the surviva sub-model summaries.
jms <- ls()[sapply(ls(), function(i) eval(parse(text=paste0("class(",i,")"))))=="joint"]

# So all significant by themselves.
for(j in jms){
  cat(j, "\n")
  print(eval(parse(text = paste0(
    "summary(", j, ")$Surv"
  ))))
  cat(cli::rule())
  cat("\n")
}

